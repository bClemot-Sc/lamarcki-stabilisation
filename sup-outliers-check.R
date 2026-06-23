### Head stabilisation in the diurnal dung beetle Kheper lamarcki
### Response outlier check for head stabilisation
### Bastien CLEMOT, Lund Vision Group, Biology Department, Lund University


# Packages and Data -------------------------------------------------------

# install.packages(c(
#   "tidyverse", "svglite"
# ))

library(tidyverse)
library(svglite)

# load processed data
df <- read.csv("data_stabilisation_processed.csv", fileEncoding = "UTF-8")


# 1 - Prepare response variables ------------------------------------------

# response variables to check
response_vars <- c(
  "head_angle_corrected",
  "thorax_angle_corrected",
  "abdomen_angle_corrected"
)

# keep only response variables present in the data
response_vars <- response_vars[response_vars %in% names(df)]

if (length(response_vars) == 0) {
  stop("No corrected response-angle columns found in the data.")
}

if (!("uid" %in% names(df))) {
  stop("Column uid is needed to save outlier flags back to the processed data.")
}

# comment column to print with flagged rows
comment_vars <- if ("comments" %in% names(df)) "comments" else character(0)

# long format
response_long <- df %>%
  mutate(
    view = recode(view, "left" = "sideL", "right" = "sideR"),
    view = factor(view, levels = c("front", "sideL", "sideR")),
    view_group = case_when(
      view == "front" ~ "front",
      view %in% c("sideL", "sideR") ~ "side",
      TRUE ~ NA_character_
    ),
    view_group = factor(view_group, levels = c("front", "side"))
  ) %>%
  select(
    any_of(c(
      "uid",
      "beetle_id",
      "condition",
      "antennae",
      "view",
      "view_group",
      "axis_entry",
      comment_vars,
      response_vars
    ))
  ) %>%
  pivot_longer(
    cols = all_of(response_vars),
    names_to = "response",
    values_to = "angle"
  ) %>%
  filter(!is.na(angle))


# 2 - Flag response outliers ----------------------------------------------

# values are flagged within each response variable
response_outliers <- response_long %>%
  group_by(response) %>%
  mutate(
    angle_mean = mean(angle, na.rm = TRUE),
    angle_sd = sd(angle, na.rm = TRUE),
    angle_z = (angle - angle_mean) / angle_sd,
    q1 = quantile(angle, 0.25, na.rm = TRUE),
    q3 = quantile(angle, 0.75, na.rm = TRUE),
    iqr = q3 - q1,
    lower_iqr = q1 - 1.5 * iqr,
    upper_iqr = q3 + 1.5 * iqr,
    flagged_z3 = abs(angle_z) > 3,
    flagged_iqr = angle < lower_iqr | angle > upper_iqr,
    flagged = flagged_z3 | flagged_iqr,
    flag_reason = case_when(
      flagged_z3 & flagged_iqr ~ "z3_and_iqr",
      flagged_z3 ~ "z3",
      flagged_iqr ~ "iqr",
      TRUE ~ "not_flagged"
    )
  ) %>%
  ungroup()

# summary printed to console
outlier_summary <- response_outliers %>%
  group_by(response) %>%
  summarise(
    n = n(),
    mean_angle = round(mean(angle, na.rm = TRUE), 2),
    sd_angle = round(sd(angle, na.rm = TRUE), 2),
    min_angle = round(min(angle, na.rm = TRUE), 2),
    max_angle = round(max(angle, na.rm = TRUE), 2),
    n_flagged_z3 = sum(flagged_z3, na.rm = TRUE),
    n_flagged_iqr = sum(flagged_iqr, na.rm = TRUE),
    n_flagged = sum(flagged, na.rm = TRUE),
    .groups = "drop"
  )

print(outlier_summary)

# flagged rows printed to console for manual inspection
flagged_outliers <- response_outliers %>%
  filter(flagged) %>%
  select(
    any_of(c(
      "uid",
      "beetle_id",
      "condition",
      "antennae",
      "view",
      "view_group",
      "axis_entry",
      comment_vars
    )),
    response,
    angle,
    angle_z,
    flag_reason,
    flagged_z3,
    flagged_iqr
  ) %>%
  arrange(response, desc(abs(angle_z)))

print(flagged_outliers, n = Inf, width = Inf)


# 3 - Exploratory plots ----------------------------------------------------

# create figure folder if needed
dir.create("raw_figures", showWarnings = FALSE)

# response boxplots
p_response_boxplot <- ggplot(
  response_long,
  aes(
    x = response,
    y = angle
  )
) +
  geom_boxplot(outlier.alpha = 0.6) +
  labs(
    title = "Response angle distributions",
    x = "Response variable",
    y = "Angle"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = "raw_figures/sup-outlier-response-boxplots.svg",
  plot = p_response_boxplot,
  width = 9,
  height = 5,
  device = svglite,
  bg = "white"
)

# response boxplots by view
p_response_by_view <- ggplot(
  response_long,
  aes(
    x = view,
    y = angle
  )
) +
  geom_boxplot(outlier.alpha = 0.6) +
  facet_wrap(~ response, scales = "free_y") +
  labs(
    title = "Response angle distributions by view",
    x = "View",
    y = "Angle"
  ) +
  theme_classic()

ggsave(
  filename = "raw_figures/super-outlier-view.svg",
  plot = p_response_by_view,
  width = 10,
  height = 7,
  device = svglite,
  bg = "white"
)

# flagged values highlighted
p_response_flagged <- ggplot(
  response_outliers,
  aes(
    x = response,
    y = angle,
    colour = flag_reason
  )
) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  scale_colour_manual(
    values = c(
      "not_flagged" = "grey70",
      "z3" = "#0072B2",
      "iqr" = "#D55E00",
      "z3_and_iqr" = "#CC79A7"
    ),
    breaks = c("z3", "iqr", "z3_and_iqr", "not_flagged")
  ) +
  labs(
    title = "Flagged response-angle values",
    subtitle = "Flagged using absolute z-score > 3 or 1.5 x IQR rule",
    x = "Response variable",
    y = "Angle",
    colour = "Flag reason"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = "raw_figures/sup-outlier-flagged.svg",
  plot = p_response_flagged,
  width = 9,
  height = 5,
  device = svglite,
  bg = "white"
)

print(p_response_boxplot)
print(p_response_by_view)
print(p_response_flagged)

# 4 - Save outlier flags ---------------------------------------------------

# response-specific outlier flags for later use
outlier_flags <- response_outliers %>%
  mutate(
    response_short = recode(
      response,
      "head_angle_corrected" = "head",
      "thorax_angle_corrected" = "thorax",
      "abdomen_angle_corrected" = "abdomen"
    )
  ) %>%
  select(
    uid,
    response_short,
    flagged_z3,
    flagged_iqr
  ) %>%
  pivot_wider(
    names_from = response_short,
    values_from = c(flagged_z3, flagged_iqr),
    names_glue = "outlier_{str_remove(.value, 'flagged_')}_{response_short}",
    values_fill = FALSE
  )

# remove old outlier flags if the script is run more than once
old_outlier_flags <- names(df) %>%
  str_subset("^outlier_(z3|iqr|flagged_z3|flagged_iqr|any)")

df <- df %>%
  select(-any_of(old_outlier_flags)) %>%
  left_join(outlier_flags, by = "uid") %>%
  mutate(
    across(starts_with("outlier_"), ~ replace_na(.x, FALSE))
  )

# overwrite processed data with added outlier flags
write.csv(
  df,
  "data_stabilisation_processed.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
