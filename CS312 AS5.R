#install packages
install.packages("optmatch")
install.packages("RItools", dependencies=T)

library(ggplot2)

library(optmatch)
library(RItools)

#read data
rm(list = ls())

library(readxl)
df <- read_excel("Airway RME matching data.xlsx")

#change column names for df
colnames(df) <- c("id", "gender", "age", "cbct_interval", "sagittal_class", "tongue_posture", "group")

#replication step 1 - optmatch

#calculate propensity score
psm <- glm(group ~ gender + age + cbct_interval + sagittal_class + tongue_posture,
           family = binomial(), data = df)

df$pscore <- psm$fitted.values
df$pscore

summary(psm)

#matching using optmatch
ps.pm = pairmatch(psm, data = df) 
summary(ps.pm)
xBalance(group ~ gender + age + cbct_interval + sagittal_class + tongue_posture, 
         ps.pm, data = df)

#compare replicated matched data id with those in the paper 
matched_indices <- which(!is.na(as.numeric(ps.pm)))
matched_indices

df2 <- read_excel("Airway RME dimensions .xlsx")
sort(df2$id)

#show the difference
df[17, ]
df[27, ]

# Plot matched and unmatched observations
df$matched <- ifelse(df$id %in% matched_indices, 1, 0)
df$matched <- as.factor(df$matched)

#histogram of treated, matched control and raw control
raw_treated <- df[df$matched == 1 & df$group == 1, ]
hist(raw_treated$pscore)

matched_treated <- df[df$matched == 1 & df$group == 1, ]
hist(matched_treated$pscore)

raw_control <- df[df$group == 0, ]
hist(raw_control$pscore)

matched_control <- df[df$matched == 1 & df$group == 0, ]
hist(matched_control$pscore)

#jitter plot of matched and unmatched data
ggplot(df, aes(x = pscore, y = group, color = matched)) +
  geom_jitter(width = 0.05, height = 0.15, alpha = 0.6) +
  labs(x = "Propensity Score", y = "Group") +
  theme_minimal() +
  theme(axis.title.y = element_blank()) +
  scale_color_manual(values = c("1" = "blue", "0" = "red")) +
  scale_y_continuous(breaks = c(0, 1), labels = c("control", "treatment"))

#matching using optmatch with replacement
ps.pm1 = pairmatch(psm, replace = TRUE, data = df)
summary(ps.pm1)
xBalance(group ~ gender + age + cbct_interval + sagittal_class + tongue_posture, 
         ps.pm1, data = df)

#xBalance without propensity score
xBalance(group ~ gender + age + cbct_interval + sagittal_class + tongue_posture, 
         data = df)


#extension - step 2 genetic matching

#clean the matched data
df2_clean <- na.omit(df2)
colnames(df2_clean)
colnames(df2_clean) <- c("id", "gender", "age", "cbct_interval", "sagittal_class", 
                         "tongue_posture", "group", "distance", "weights", "subclass",
                         "Nasopharyngeal_airway_Volumn_T0", 
                         "Retropalatal_airway_Volumn_T0", "Retropalatal_airway_MCA_T0",
                         "Retroglossal_airway_Volumn_T0", "Retroglossal_airway_MCA_T0",
                         "Nasopharyngeal_airway_Volumn_T1",                         
                         "Retropalatal_airway_Volumn_T1", "Retropalatal_airway_MCA_T1",
                         "Retroglossal_airway_Volumn_T1", "Retroglossal_airway_MCA_T1")
str(df2_clean)
df2_clean$Retropalatal_airway_MCA_T0 <- as.numeric(df2_clean$Retropalatal_airway_MCA_T0)
df2_clean$Retropalatal_airway_MCA_T1 <- as.numeric(df2_clean$Retropalatal_airway_MCA_T1)
df2_clean$Retroglossal_airway_MCA_T0 <- as.numeric(df2_clean$Retroglossal_airway_MCA_T0)
df2_clean$Retroglossal_airway_MCA_T1 <- as.numeric(df2_clean$Retroglossal_airway_MCA_T1)

#genetic matching
library(Matching)
set.seed(123)

gen_match <- GenMatch(Tr = df2_clean$group,
                      X = df2_clean[, c("gender", "age", "cbct_interval", "sagittal_class", "tongue_posture", "Retropalatal_airway_MCA_T0")],
                      wait.generations=10,
                      M=1, pop.size = 1000)

matched <- Match(Y = df2_clean$Retropalatal_airway_MCA_T1, Tr = df2_clean$group, 
                      X = df2_clean[, c("gender", "age", "cbct_interval", "sagittal_class", "tongue_posture", "Retropalatal_airway_MCA_T0")],
                      Weight.matrix = gen_match)
summary(matched)

balance_check <- MatchBalance(group ~ gender + age + cbct_interval + sagittal_class + tongue_posture + Retropalatal_airway_MCA_T0, 
                              data = df2_clean, match.out = matched)

gen_match2 <- GenMatch(Tr = df2_clean$group,
                      X = df2_clean[, c("gender", "age", "cbct_interval", "sagittal_class", "tongue_posture", "Retroglossal_airway_MCA_T0")],
                      exact = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE),
                      M=1, pop.size = 1000)

matched2 <- Match(Y = df2_clean$Retroglossal_airway_MCA_T1, Tr = df2_clean$group, 
                 X = df2_clean[, c("gender", "age", "cbct_interval", "sagittal_class", "tongue_posture", "Retroglossal_airway_MCA_T0")],
                 Weight.matrix = gen_match2)
summary(matched2)

MatchBalance(group ~ gender + age + cbct_interval + sagittal_class + tongue_posture + Retroglossal_airway_MCA_T0, 
                              data = df2_clean, match.out = matched2)

#For Figure 3 and 4

install.packages("MatchIt")
library(MatchIt)

#load dataset
ds2 <- read_excel("Airway RME matching data.xlsx")

#column name change
colnames(ds2) <- c("id", "gender", "age", "cbct_interval", "sagittal_class", "tongue_posture", "group")

#optimal matching without replacement
matchingdata <- matchit(group ~ gender + age + cbct_interval + sagittal_class + tongue_posture, data = ds2, method = "optimal", ratio = 1)
summary(matchingdata)

#hist plot
plot(matchingdata, type ="hist")
#jitter plot
plot(matchingdata, type ="jitter")


