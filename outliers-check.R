# Head stabilisation in Kheper lamarcki - outlier flagging
# Flags head, thorax and abdomen angles with two rules, within each response: |z| > 3 and 1.5 x IQR
# The flags are saved in the processed data, main-analysis.R decides which values are removed

library(tidyverse)
library(svglite)

data_file <- "data_stabilisation_processed.csv"

dir.create("raw_figures", showWarnings = FALSE)

# 1. Data ------------------------------------------------------------------

df <- read.csv(data_file, fileEncoding = "UTF-8")

# One row per angle (head, thorax, abdomen)
response_long <- df %>%
  mutate(view = factor(recode(view, "left" = "sideL", "right" = "sideR"),
                       levels = c("front", "sideL", "sideR"))) %>%
  rename(head = head_angle_corrected, thorax = thorax_angle_corrected,
         abdomen = abdomen_angle_corrected) %>%
  select(uid, beetle_id, condition, antennae, view, axis_entry, comments, head, thorax, abdomen) %>%
  pivot_longer(c(head, thorax, abdomen), names_to = "response", values_to = "angle") %>%
  filter(!is.na(angle))

# 2. Flag outliers ---------------------------------------------------------

# Values are flagged within each response variable
outliers <- response_long %>%
  group_by(response) %>%
  mutate(
    angle_z = (angle - mean(angle)) / sd(angle),
    iqr_low = quantile(angle, 0.25) - 1.5 * IQR(angle),
    iqr_high = quantile(angle, 0.75) + 1.5 * IQR(angle),
    flagged_z3 = abs(angle_z) > 3,
    flagged_iqr = angle < iqr_low | angle > iqr_high,
    flagged = flagged_z3 | flagged_iqr,
    flag_reason = case_when(
      flagged_z3 & flagged_iqr ~ "z3_and_iqr",
      flagged_z3 ~ "z3",
      flagged_iqr ~ "iqr",
      TRUE ~ "not_flagged"
    )
  ) %>%
  ungroup()

# Summary and flagged rows for manual inspection
outliers %>%
  group_by(response) %>%
  summarise(n = n(), mean_angle = mean(angle), sd_angle = sd(angle),
            min_angle = min(angle), max_angle = max(angle),
            n_flagged_z3 = sum(flagged_z3), n_flagged_iqr = sum(flagged_iqr),
            n_flagged = sum(flagged)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  print()

outliers %>%
  filter(flagged) %>%
  select(uid, beetle_id, condition, antennae, view, axis_entry, comments,
         response, angle, angle_z, flag_reason) %>%
  arrange(response, desc(abs(angle_z))) %>%
  print(n = Inf, width = Inf)

# 3. Plots -----------------------------------------------------------------

save_svg <- function(plot, name, width = 9, height = 5) {
  ggsave(paste0("raw_figures/", name, ".svg"), plot,
         width = width, height = height, device = svglite, bg = "white")
}

# Flagged values highlighted
p_flagged <- ggplot(outliers, aes(x = response, y = angle, colour = flag_reason)) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  scale_colour_manual(
    values = c(not_flagged = "grey70", z3 = "#0072B2", iqr = "#D55E00", z3_and_iqr = "#CC79A7"),
    breaks = c("z3", "iqr", "z3_and_iqr", "not_flagged")
  ) +
  labs(title = "Flagged response-angle values",
       subtitle = "Flagged using absolute z-score > 3 or 1.5 x IQR rule",
       x = "Response variable", y = "Angle", colour = "Flag reason") +
  theme_classic()
save_svg(p_flagged, "sup-outlier-flagged")

# 4. Save outlier flags ----------------------------------------------------

# One column per rule and response (outlier_z3_head, outlier_iqr_head, ...)
outlier_flags <- outliers %>%
  select(uid, response, flagged_z3, flagged_iqr) %>%
  pivot_wider(names_from = response, values_from = c(flagged_z3, flagged_iqr),
              names_glue = "outlier_{str_remove(.value, 'flagged_')}_{response}",
              values_fill = FALSE)

# Remove old flags so the script can be run more than once, then add the new ones to the processed data
df <- df %>%
  select(-any_of(str_subset(names(df), "^outlier_(z3|iqr|flagged_z3|flagged_iqr|any)"))) %>%
  left_join(outlier_flags, by = "uid") %>%
  mutate(across(starts_with("outlier_"), ~ replace_na(.x, FALSE)))

write.csv(df, data_file, row.names = FALSE, fileEncoding = "UTF-8")