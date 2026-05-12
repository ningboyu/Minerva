# Thesis analysis script
# Author: Boyu Ning
# 
# Note:
# This script was used for the thesis analysis.
# The underlying article-level dataset is not included in this public repository
# because it contains private platform and business-related data.
#
# To reproduce the analysis, place the clean analytic dataset at:
# data/master_data_20260425_final.xlsx

# -------------------------------
# 0. Packages
# -------------------------------

library(readxl)
library(dplyr)
library(tidyr)
library(MatchIt)
library(ggplot2)

# -------------------------------
# 1. Load data
# -------------------------------

df <- read_excel("data/master_data_20260425_final.xlsx", sheet = "Master Data")

# -------------------------------
# 2. Data cleaning and variable construction
# -------------------------------

df_clean <- df %>%
  mutate(
    # Stable article identifier for matching and pair-level analysis
    row_id = row_number(),
    
    # Treatment variable
    # 0 = human-only, 1 = AI-assisted
    treatment = as.numeric(treatment_group),
    treatment_group = factor(
      treatment_group,
      levels = c(0, 1),
      labels = c("human", "human_ai")
    ),
    
    # Matching variables
    theme = factor(theme),
    cta_binary = factor(cta_binary),
    
    # Numeric variables
    word_count = as.numeric(word_count),
    delivered_users = as.numeric(delivered_users),
    reads = as.numeric(reads),
    shares = as.numeric(shares),
    likes = as.numeric(likes),
    new_follows = as.numeric(new_follows),
    comments = as.numeric(comments),
    
    # One completion-rate value is recorded as "na" in the source data.
    # It is treated as missing before numeric conversion.
    completion_rate = na_if(as.character(completion_rate), "na"),
    completion_rate = as.numeric(completion_rate)
  ) %>%
  mutate(
    # Pre-click metric
    view_rate = if_else(delivered_users > 0, reads / delivered_users, NA_real_),
    
    # Post-click engagement metric
    total_engagement = rowSums(
      across(c(shares, likes, comments, new_follows)),
      na.rm = TRUE
    ),
    engagement_rate = if_else(reads > 0, total_engagement / reads, NA_real_)
  )


# -------------------------------
# Table 3: Theme distribution by workflow condition before matching
# -------------------------------

table3_theme <- df_clean %>%
  count(treatment_group, theme) %>%
  group_by(treatment_group) %>%
  mutate(
    pct = n / sum(n),
    n_pct = paste0(n, " (", round(100 * pct, 1), "%)")
  ) %>%
  ungroup() %>%
  select(theme, treatment_group, n_pct) %>%
  pivot_wider(
    names_from = treatment_group,
    values_from = n_pct,
    values_fill = "0 (0.0%)"
  ) %>%
  rename(
    Theme = theme,
    `Human-only n (%)` = human,
    `AI-assisted n (%)` = human_ai
  )

View(table3_theme)

# -------------------------------
# Table 4: Paper-style descriptive statistics before matching
# Numeric values only
# -------------------------------

table4 <- df_clean %>%
  mutate(
    cta_presence = as.numeric(as.character(cta_binary))
  ) %>%
  select(
    treatment_group,
    word_count,
    cta_presence,
    view_rate,
    engagement_rate,
    completion_rate
  ) %>%
  pivot_longer(
    cols = -treatment_group,
    names_to = "variable",
    values_to = "value"
  ) %>%
  group_by(treatment_group, variable) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Variable = recode(
      variable,
      "word_count" = "Word Count",
      "cta_presence" = "CTA Presence",
      "view_rate" = "View Rate",
      "engagement_rate" = "Engagement Rate",
      "completion_rate" = "Completion Rate"
    )
  ) %>%
  select(-variable) %>%
  pivot_wider(
    names_from = treatment_group,
    values_from = c(mean, sd)
  ) %>%
  transmute(
    Variable,
    `Human-only Mean` = mean_human,
    `Human-only SD` = sd_human,
    `AI-assisted Mean` = mean_human_ai,
    `AI-assisted SD` = sd_human_ai
  ) %>%
  mutate(
    Variable = factor(
      Variable,
      levels = c(
        "Word Count",
        "CTA Presence",
        "View Rate",
        "Engagement Rate",
        "Completion Rate"
      )
    )
  ) %>%
  arrange(Variable) %>%
  mutate(Variable = as.character(Variable))

table4_view <- table4 %>%
  mutate(
    `Human-only Mean` = round(`Human-only Mean`, 4),
    `Human-only SD` = round(`Human-only SD`, 4),
    `AI-assisted Mean` = round(`AI-assisted Mean`, 4),
    `AI-assisted SD` = round(`AI-assisted SD`, 4)
  )

View(table4_view)

# -------------------------------
# 3. Matching: Nearest Neighbor Matching
# -------------------------------

set.seed(123)

m.out <- matchit(
  treatment ~ word_count + theme + cta_binary,
  data = df_clean,
  method = "nearest",
  exact = ~ theme + cta_binary,
  replace = TRUE
)

# Extract matched dataset
matched_data <- match.data(m.out)

# -------------------------------
# 4. Matching balance summary
# -------------------------------

match_matrix <- m.out$match.matrix

matched_pairs <- tibble(
  treated_row = as.integer(rownames(match_matrix)),
  control_row = as.integer(match_matrix[, 1])
) %>%
  left_join(
    df_clean %>%
      select(
        row_id,
        ai_theme = theme,
        ai_cta = cta_binary,
        ai_word_count = word_count,
        ai_engagement_rate = engagement_rate,
        ai_completion_rate = completion_rate,
        ai_view_rate = view_rate
      ),
    by = c("treated_row" = "row_id")
  ) %>%
  left_join(
    df_clean %>%
      select(
        row_id,
        human_theme = theme,
        human_cta = cta_binary,
        human_word_count = word_count,
        human_engagement_rate = engagement_rate,
        human_completion_rate = completion_rate,
        human_view_rate = view_rate
      ),
    by = c("control_row" = "row_id")
  ) %>%
  mutate(
    abs_word_count_diff = abs(ai_word_count - human_word_count),
    diff_engagement_rate = ai_engagement_rate - human_engagement_rate,
    diff_completion_rate = ai_completion_rate - human_completion_rate,
    diff_view_rate = ai_view_rate - human_view_rate
  )

balance_summary <- matched_pairs %>%
  summarise(
    n_pairs = n(),
    unique_controls = n_distinct(control_row),
    max_control_reuse = max(as.numeric(table(control_row))),
    exact_theme_match_rate = mean(ai_theme == human_theme),
    exact_cta_match_rate = mean(ai_cta == human_cta),
    mean_abs_word_count_diff = mean(abs_word_count_diff, na.rm = TRUE),
    median_abs_word_count_diff = median(abs_word_count_diff, na.rm = TRUE),
    max_abs_word_count_diff = max(abs_word_count_diff, na.rm = TRUE)
  )

balance_summary %>%
  as.data.frame()

# -------------------------------
# 5. Treatment effect estimate: Engagement Rate
# -------------------------------

estimate_effect <- function(diff_vector) {
  diff_vector <- diff_vector[!is.na(diff_vector)]
  test <- t.test(diff_vector, mu = 0)
  
  tibble(
    n_pairs = length(diff_vector),
    estimate = mean(diff_vector),
    se = sd(diff_vector) / sqrt(length(diff_vector)),
    t_statistic = unname(test$statistic),
    p_value = test$p.value,
    ci_low = test$conf.int[1],
    ci_high = test$conf.int[2]
  )
}

engagement_effect <- estimate_effect(matched_pairs$diff_engagement_rate)

engagement_effect


# -------------------------------
# 6. Supporting outcomes: Completion Rate and View Rate
# -------------------------------

completion_effect <- estimate_effect(matched_pairs$diff_completion_rate)
view_effect <- estimate_effect(matched_pairs$diff_view_rate)

supporting_effects <- bind_rows(
  engagement_effect %>% mutate(outcome = "Engagement Rate"),
  completion_effect %>% mutate(outcome = "Completion Rate"),
  view_effect %>% mutate(outcome = "View Rate")
) %>%
  select(outcome, n_pairs, estimate, se, t_statistic, p_value, ci_low, ci_high)

supporting_effects

# -------------------------------
# 7. Full-sample differences vs matched differences
# -------------------------------

full_means <- df_clean %>%
  group_by(treatment_group) %>%
  summarise(
    mean_engagement_rate = mean(engagement_rate, na.rm = TRUE),
    mean_completion_rate = mean(completion_rate, na.rm = TRUE),
    mean_view_rate = mean(view_rate, na.rm = TRUE),
    .groups = "drop"
  )

full_differences <- tibble(
  outcome = c("Engagement Rate", "Completion Rate", "View Rate"),
  full_difference = c(
    full_means$mean_engagement_rate[full_means$treatment_group == "human_ai"] -
      full_means$mean_engagement_rate[full_means$treatment_group == "human"],
    full_means$mean_completion_rate[full_means$treatment_group == "human_ai"] -
      full_means$mean_completion_rate[full_means$treatment_group == "human"],
    full_means$mean_view_rate[full_means$treatment_group == "human_ai"] -
      full_means$mean_view_rate[full_means$treatment_group == "human"]
  )
)

matched_differences <- supporting_effects %>%
  select(outcome, matched_difference = estimate)

full_vs_matched <- full_differences %>%
  left_join(matched_differences, by = "outcome") %>%
  mutate(
    change_after_matching = matched_difference - full_difference
  )

full_vs_matched

# -------------------------------
# 8. Pair-level performance dispersion
# -------------------------------

cv_pair_summary <- tibble(
  group = c("Matched human-only", "AI-assisted"),
  n = c(
    sum(!is.na(matched_pairs$human_engagement_rate)),
    sum(!is.na(matched_pairs$ai_engagement_rate))
  ),
  mean_engagement_rate = c(
    mean(matched_pairs$human_engagement_rate, na.rm = TRUE),
    mean(matched_pairs$ai_engagement_rate, na.rm = TRUE)
  ),
  sd_engagement_rate = c(
    sd(matched_pairs$human_engagement_rate, na.rm = TRUE),
    sd(matched_pairs$ai_engagement_rate, na.rm = TRUE)
  ),
  mean_completion_rate = c(
    mean(matched_pairs$human_completion_rate, na.rm = TRUE),
    mean(matched_pairs$ai_completion_rate, na.rm = TRUE)
  ),
  sd_completion_rate = c(
    sd(matched_pairs$human_completion_rate, na.rm = TRUE),
    sd(matched_pairs$ai_completion_rate, na.rm = TRUE)
  )
) %>%
  mutate(
    cv_engagement_rate = sd_engagement_rate / mean_engagement_rate,
    cv_completion_rate = sd_completion_rate / mean_completion_rate
  )

cv_pair_summary %>%
  as.data.frame()

# -------------------------------
# 9. Bootstrap CI for CV difference
# -------------------------------

set.seed(123)

cv <- function(x) {
  x <- x[!is.na(x)]
  sd(x) / mean(x)
}

bootstrap_cv_difference <- function(ai, human, n_boot = 10000) {
  ai <- ai[!is.na(ai)]
  human <- human[!is.na(human)]
  
  boot_diffs <- replicate(n_boot, {
    ai_sample <- sample(ai, size = length(ai), replace = TRUE)
    human_sample <- sample(human, size = length(human), replace = TRUE)
    
    cv(ai_sample) - cv(human_sample)
  })
  
  tibble(
    cv_ai = cv(ai),
    cv_human = cv(human),
    cv_difference = cv(ai) - cv(human),
    ci_low = quantile(boot_diffs, 0.025),
    ci_high = quantile(boot_diffs, 0.975),
    n_boot = n_boot
  )
}

cv_bootstrap_engagement <- bootstrap_cv_difference(
  ai = matched_pairs$ai_engagement_rate,
  human = matched_pairs$human_engagement_rate,
  n_boot = 10000
)

cv_bootstrap_engagement

# -------------------------------
# 10. Robustness check: matching without replacement
# -------------------------------

set.seed(123)

m.out_no_replace <- matchit(
  treatment ~ word_count + theme + cta_binary,
  data = df_clean,
  method = "nearest",
  exact = ~ theme + cta_binary,
  replace = FALSE
)

match_matrix_no_replace <- m.out_no_replace$match.matrix

matched_pairs_no_replace <- tibble(
  treated_row = as.integer(rownames(match_matrix_no_replace)),
  control_row = as.integer(match_matrix_no_replace[, 1])
) %>%
  filter(!is.na(control_row)) %>%
  left_join(
    df_clean %>%
      select(
        row_id,
        ai_theme = theme,
        ai_cta = cta_binary,
        ai_word_count = word_count,
        ai_engagement_rate = engagement_rate,
        ai_completion_rate = completion_rate,
        ai_view_rate = view_rate
      ),
    by = c("treated_row" = "row_id")
  ) %>%
  left_join(
    df_clean %>%
      select(
        row_id,
        human_theme = theme,
        human_cta = cta_binary,
        human_word_count = word_count,
        human_engagement_rate = engagement_rate,
        human_completion_rate = completion_rate,
        human_view_rate = view_rate
      ),
    by = c("control_row" = "row_id")
  ) %>%
  mutate(
    abs_word_count_diff = abs(ai_word_count - human_word_count),
    diff_engagement_rate = ai_engagement_rate - human_engagement_rate,
    diff_completion_rate = ai_completion_rate - human_completion_rate,
    diff_view_rate = ai_view_rate - human_view_rate
  )

balance_summary_no_replace <- matched_pairs_no_replace %>%
  summarise(
    n_pairs = n(),
    unique_controls = n_distinct(control_row),
    max_control_reuse = max(as.numeric(table(control_row))),
    exact_theme_match_rate = mean(ai_theme == human_theme),
    exact_cta_match_rate = mean(ai_cta == human_cta),
    mean_abs_word_count_diff = mean(abs_word_count_diff, na.rm = TRUE),
    median_abs_word_count_diff = median(abs_word_count_diff, na.rm = TRUE),
    max_abs_word_count_diff = max(abs_word_count_diff, na.rm = TRUE)
  )

engagement_effect_no_replace <- estimate_effect(
  matched_pairs_no_replace$diff_engagement_rate
)

balance_summary_no_replace%>%
  as.data.frame()
engagement_effect_no_replace%>%
  as.data.frame()

engagement_robustness_check <- bind_rows(
  engagement_effect %>%
    mutate(specification = "Main matching, with replacement"),
  engagement_effect_no_replace %>%
    mutate(specification = "Robustness check, without replacement")
) %>%
  select(
    specification,
    n_pairs,
    estimate,
    se,
    t_statistic,
    p_value,
    ci_low,
    ci_high
  )

engagement_robustness_check