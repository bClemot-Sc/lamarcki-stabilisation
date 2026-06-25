### Head stabilisation in the diurnal dung beetle Kheper lamarcki
### Main analysis with cleaned outlier values and external-factor checks
### Bastien CLEMOT, Lund Vision Group, Biology Department, Lund University


# Packages and Data -------------------------------------------------------

# install.packages(c(
#   "tidyverse", "lme4", "lmerTest", "emmeans", "svglite"
# ))

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(svglite)

# load processed data
df <- read.csv("data_stabilisation_processed.csv", fileEncoding = "UTF-8")

dir.create("raw_figures", showWarnings = FALSE)


# 0 - Prepare analysis data -----------------------------------------------

analysis_data <- df %>%
  mutate(
    condition = factor(condition, levels = c("flat0", "tilted10")),
    antennae = factor(antennae, levels = c("both", "right", "left", "none")),
    view = recode(view, "left" = "sideL", "right" = "sideR"),
    view = factor(view, levels = c("front", "sideL", "sideR")),
    rotation_plane = case_when(
      view == "front" ~ "roll",
      view %in% c("sideL", "sideR") ~ "pitch"
    ),
    rotation_plane = factor(rotation_plane, levels = c("roll", "pitch")),
    date = factor(date),
    axis_entry = factor(axis_entry, levels = c("butt", "head")),
    has_comment = !is.na(comments) & str_trim(comments) != "",
    outlier_remove_head = (outlier_z3_head & outlier_iqr_head) |
      ((outlier_z3_head | outlier_iqr_head) & has_comment),
    outlier_remove_thorax = (outlier_z3_thorax & outlier_iqr_thorax) |
      ((outlier_z3_thorax | outlier_iqr_thorax) & has_comment),
    outlier_remove_abdomen = (outlier_z3_abdomen & outlier_iqr_abdomen) |
      ((outlier_z3_abdomen | outlier_iqr_abdomen) & has_comment),
    head_angle_clean = if_else(outlier_remove_head, NA_real_, head_angle_corrected),
    thorax_angle_clean = if_else(outlier_remove_thorax, NA_real_, thorax_angle_corrected),
    abdomen_angle_clean = if_else(outlier_remove_abdomen, NA_real_, abdomen_angle_corrected)
  )

cat("\n===== CLEANING SUMMARY =====\n")
cat("Rows in processed data:", nrow(analysis_data), "\n")
cat("Head values set to NA:", sum(analysis_data$outlier_remove_head, na.rm = TRUE), "\n")
cat("Thorax values set to NA:", sum(analysis_data$outlier_remove_thorax, na.rm = TRUE), "\n")
cat("Abdomen values set to NA:", sum(analysis_data$outlier_remove_abdomen, na.rm = TRUE), "\n\n")

cat("===== EXTERNAL FACTORS USED IN ADJUSTED MODELS =====\n")
cat("Rows missing date:", sum(is.na(analysis_data$date)), "\n")
cat("Rows missing time_decimal_c:", sum(is.na(analysis_data$time_decimal_c)), "\n")
cat("Rows missing axis_entry:", sum(is.na(analysis_data$axis_entry)), "\n")
cat("Rows missing progress_c:", sum(is.na(analysis_data$progress_c)), "\n")

# bootstrap confidence interval across beetles
boot_ci <- function(x, stat_fun, n_boot = 5000) {
  boot_values <- replicate(
    n_boot,
    stat_fun(sample(x, size = length(x), replace = TRUE))
  )
  unname(quantile(boot_values, c(0.025, 0.975), na.rm = TRUE))
}

# summarize values with bootstrap confidence intervals
boot_summary <- function(x, n_boot = 5000) {
  mean_ci <- boot_ci(x, mean, n_boot)
  median_ci <- boot_ci(x, median, n_boot)
  
  tibble(
    statistic = c("mean", "median"),
    estimate = c(mean(x), median(x)),
    ci_low = c(mean_ci[1], median_ci[1]),
    ci_high = c(mean_ci[2], median_ci[2])
  )
}

# format p-values for printed tables
format_p <- function(p) {
  formatC(p, format = "e", digits = 2)
}

# external factors used in adjusted models
external_terms <- "date + time_decimal_c + axis_entry * (progress_c + progress_c2)"

# combine raw and Bonferroni-adjusted contrast p-values
combine_contrasts <- function(contrast_object, group_cols = character()) {
  raw <- as.data.frame(summary(contrast_object, infer = c(TRUE, TRUE), adjust = "none")) %>%
    select(
      all_of(group_cols),
      contrast,
      estimate,
      SE,
      df,
      lower.CL,
      upper.CL,
      t.ratio,
      p_raw = p.value
    )
  
  bonferroni <- as.data.frame(summary(contrast_object, infer = c(FALSE, TRUE), adjust = "bonferroni")) %>%
    select(
      all_of(group_cols),
      contrast,
      p_bonferroni = p.value
    )
  
  raw %>%
    left_join(bonferroni, by = c(group_cols, "contrast"))
}

# print contrast tables
print_contrasts <- function(data) {
  print(
    as_tibble(data %>%
                mutate(
                  across(
                    c(estimate, SE, df, lower.CL, upper.CL, t.ratio),
                    ~ round(.x, 3)
                  ),
                  p_raw = format_p(p_raw),
                  p_bonferroni = format_p(p_bonferroni)
                )),
    n = Inf
  )
}


# 1 - Body segment movement during flat-horizon dances --------------------

# flat horizon with both antennae
movement_data <- analysis_data %>%
  filter(
    condition == "flat0",
    antennae == "both"
  ) %>%
  select(
    uid,
    beetle_id,
    rotation_plane,
    view,
    head_angle_clean,
    thorax_angle_clean,
    abdomen_angle_clean
  ) %>%
  pivot_longer(
    cols = c(
      head_angle_clean,
      thorax_angle_clean,
      abdomen_angle_clean
    ),
    names_to = "body_segment",
    values_to = "angle"
  ) %>%
  mutate(
    body_segment = recode(
      body_segment,
      "head_angle_clean" = "head",
      "thorax_angle_clean" = "thorax",
      "abdomen_angle_clean" = "abdomen"
    ),
    body_segment = factor(body_segment, levels = c("head", "thorax", "abdomen"))
  ) %>%
  filter(!is.na(angle))

# center each body segment around each beetle's own mean
movement_centered <- movement_data %>%
  group_by(rotation_plane, body_segment, beetle_id) %>%
  mutate(
    beetle_mean = mean(angle),
    angle_centered = angle - beetle_mean
  ) %>%
  ungroup()

# within-beetle movement range
beetle_movement <- movement_centered %>%
  group_by(rotation_plane, body_segment, beetle_id) %>%
  summarise(
    n_obs = n(),
    min_angle = min(angle_centered),
    max_angle = max(angle_centered),
    movement_range = max_angle - min_angle,
    .groups = "drop"
  )

# sample sizes
movement_n <- movement_centered %>%
  group_by(rotation_plane, body_segment) %>%
  summarise(
    n_obs = n(),
    n_beetles = n_distinct(beetle_id),
    .groups = "drop"
  )

# bootstrap summaries across beetles
movement_summary <- beetle_movement %>%
  group_by(rotation_plane, body_segment) %>%
  group_modify(~ boot_summary(.x$movement_range)) %>%
  ungroup() %>%
  left_join(movement_n, by = c("rotation_plane", "body_segment")) %>%
  select(
    rotation_plane,
    body_segment,
    statistic,
    estimate,
    ci_low,
    ci_high,
    n_obs,
    n_beetles
  )

cat("\n===== BODY SEGMENT MOVEMENT RANGE — BOOTSTRAP SUMMARY =====\n")
cat("Movement range is estimated within beetles after centering each beetle around its own mean.\n")
cat("Confidence intervals are bootstrap 95% CIs across beetles.\n\n")
print(
  as_tibble(movement_summary %>%
              mutate(
                estimate = round(estimate, 2),
                ci_low = round(ci_low, 2),
                ci_high = round(ci_high, 2)
              )),
  n = Inf
)

# compare movement ranges between body segments
movement_range_contrasts <- map_dfr(
  levels(beetle_movement$rotation_plane),
  function(plane_label) {
    plane_data <- beetle_movement %>%
      filter(rotation_plane == plane_label) %>%
      mutate(body_segment = droplevels(body_segment))
    
    fit <- lmer(
      movement_range ~ body_segment + (1 | beetle_id),
      data = plane_data,
      REML = TRUE
    )
    
    emm <- emmeans(fit, ~ body_segment)
    segment_contrasts <- contrast(
      emm,
      method = "pairwise",
      adjust = "none"
    )
    
    combine_contrasts(segment_contrasts) %>%
      mutate(
        rotation_plane = plane_label,
        n_obs = nobs(fit),
        n_beetles = n_distinct(fit@frame$beetle_id)
      )
  }
) %>%
  select(
    rotation_plane,
    contrast,
    estimate,
    SE,
    df,
    lower.CL,
    upper.CL,
    t.ratio,
    p_raw,
    p_bonferroni,
    n_obs,
    n_beetles
  )

cat("\n===== BODY SEGMENT MOVEMENT RANGE — MODEL CONTRASTS =====\n")
cat("Movement ranges are compared at the beetle level within each rotation plane.\n")
cat("Negative estimates mean that the first body segment in the contrast has a smaller movement range than the second.\n")
cat("Both raw and Bonferroni-adjusted p-values are shown.\n\n")
print(
  as_tibble(movement_range_contrasts %>%
              mutate(
                across(
                  c(estimate, SE, df, lower.CL, upper.CL, t.ratio),
                  ~ round(.x, 3)
                ),
                p_raw = format_p(p_raw),
                p_bonferroni = format_p(p_bonferroni)
              )),
  n = Inf
)

# centered angle distribution
p_movement_centered <- ggplot(
  movement_centered,
  aes(
    x = body_segment,
    y = angle_centered
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_violin(fill = "grey85", colour = NA, alpha = 0.8, trim = TRUE) +
  geom_jitter(width = 0.12, alpha = 0.20, size = 0.8) +
  geom_point(
    data = beetle_movement,
    aes(y = min_angle),
    colour = "#0072B2",
    alpha = 0.8,
    size = 1.8,
    position = position_jitter(width = 0.10, height = 0)
  ) +
  geom_point(
    data = beetle_movement,
    aes(y = max_angle),
    colour = "#D55E00",
    alpha = 0.8,
    size = 1.8,
    position = position_jitter(width = 0.10, height = 0)
  ) +
  facet_wrap(~ rotation_plane, scales = "free_x") +
  labs(
    title = "Centered body-segment angles during flat-horizon dances",
    subtitle = "Grey = frame-level values; blue/orange = each beetle's minimum and maximum centered angle",
    x = "Body segment",
    y = "Centered angle (deg)"
  ) +
  theme_classic()

# beetle-level movement range
movement_summary_plot <- movement_summary %>%
  mutate(
    statistic = factor(statistic, levels = c("mean", "median"))
  )

p_movement_range <- ggplot() +
  geom_jitter(
    data = beetle_movement,
    aes(
      x = body_segment,
      y = movement_range
    ),
    width = 0.12,
    alpha = 0.45,
    size = 1.8
  ) +
  geom_pointrange(
    data = movement_summary_plot,
    aes(
      x = body_segment,
      y = estimate,
      ymin = ci_low,
      ymax = ci_high,
      shape = statistic
    ),
    position = position_dodge(width = 0.35),
    size = 0.7,
    linewidth = 0.8
  ) +
  facet_wrap(~ rotation_plane, scales = "free_x") +
  labs(
    title = "Within-beetle body-segment movement range",
    subtitle = "Points are individual beetles; intervals are bootstrap 95% CIs across beetles",
    x = "Body segment",
    y = "Movement range (deg)",
    shape = "Summary"
  ) +
  theme_classic()

# save plots
ggsave(
  filename = "raw_figures/main-1-centered-segment-angles.svg",
  plot = p_movement_centered,
  width = 9,
  height = 5,
  device = svglite,
  bg = "white"
)

ggsave(
  filename = "raw_figures/main-1-segment-movement-range.svg",
  plot = p_movement_range,
  width = 9,
  height = 5,
  device = svglite,
  bg = "white"
)

print(p_movement_centered)
print(p_movement_range)

# clean section 1 objects
rm(
  movement_data,
  movement_centered,
  beetle_movement,
  movement_n,
  movement_summary,
  movement_range_contrasts,
  movement_summary_plot,
  p_movement_centered,
  p_movement_range
)


# 2 - Neck movement during flat-horizon dances ----------------------------

# flat horizon with both antennae
neck_data <- analysis_data %>%
  filter(
    condition == "flat0",
    antennae == "both",
    !is.na(head_angle_clean),
    !is.na(thorax_angle_clean)
  ) %>%
  mutate(
    neck_angle = head_angle_clean - thorax_angle_clean
  ) %>%
  select(
    uid,
    beetle_id,
    rotation_plane,
    view,
    neck_angle
  )

# center neck angle around each beetle's own mean
neck_centered <- neck_data %>%
  group_by(rotation_plane, beetle_id) %>%
  mutate(
    beetle_mean = mean(neck_angle),
    neck_centered = neck_angle - beetle_mean
  ) %>%
  ungroup()

# within-beetle neck movement range
beetle_neck <- neck_centered %>%
  group_by(rotation_plane, beetle_id) %>%
  summarise(
    n_obs = n(),
    min_neck = min(neck_centered),
    max_neck = max(neck_centered),
    neck_range = max_neck - min_neck,
    .groups = "drop"
  )

# sample sizes
neck_n <- neck_centered %>%
  group_by(rotation_plane) %>%
  summarise(
    n_obs = n(),
    n_beetles = n_distinct(beetle_id),
    .groups = "drop"
  )

# bootstrap summaries across beetles
neck_summary <- beetle_neck %>%
  group_by(rotation_plane) %>%
  group_modify(~ boot_summary(.x$neck_range)) %>%
  ungroup() %>%
  left_join(neck_n, by = "rotation_plane") %>%
  select(
    rotation_plane,
    statistic,
    estimate,
    ci_low,
    ci_high,
    n_obs,
    n_beetles
  )

cat("\n===== NECK MOVEMENT RANGE — BOOTSTRAP SUMMARY =====\n")
cat("Neck angle is head angle minus thorax angle.\n")
cat("Outlier head and thorax values were set to NA before computing neck angle.\n")
cat("Movement range is estimated within beetles after centering each beetle around its own mean.\n")
cat("Confidence intervals are bootstrap 95% CIs across beetles.\n\n")
print(
  as_tibble(neck_summary %>%
              mutate(
                estimate = round(estimate, 2),
                ci_low = round(ci_low, 2),
                ci_high = round(ci_high, 2)
              )),
  n = Inf
)

# compare neck movement range to body-segment movement ranges
segment_neck_body_ranges <- analysis_data %>%
  filter(
    condition == "flat0",
    antennae == "both"
  ) %>%
  select(
    beetle_id,
    rotation_plane,
    head_angle_clean,
    thorax_angle_clean,
    abdomen_angle_clean
  ) %>%
  pivot_longer(
    cols = c(
      head_angle_clean,
      thorax_angle_clean,
      abdomen_angle_clean
    ),
    names_to = "range_type",
    values_to = "angle"
  ) %>%
  mutate(
    range_type = recode(
      range_type,
      "head_angle_clean" = "head",
      "thorax_angle_clean" = "thorax",
      "abdomen_angle_clean" = "abdomen"
    ),
    range_type = factor(range_type, levels = c("head", "thorax", "abdomen", "neck"))
  ) %>%
  filter(!is.na(angle)) %>%
  group_by(rotation_plane, range_type, beetle_id) %>%
  mutate(angle_centered = angle - mean(angle)) %>%
  summarise(
    movement_range = max(angle_centered) - min(angle_centered),
    .groups = "drop"
  )

segment_neck_ranges <- bind_rows(
  segment_neck_body_ranges,
  beetle_neck %>%
    transmute(
      rotation_plane,
      beetle_id,
      range_type = factor("neck", levels = c("head", "thorax", "abdomen", "neck")),
      movement_range = neck_range
    )
)

segment_neck_range_contrasts <- map_dfr(
  levels(segment_neck_ranges$rotation_plane),
  function(plane_label) {
    plane_data <- segment_neck_ranges %>%
      filter(rotation_plane == plane_label) %>%
      mutate(range_type = droplevels(range_type))
    
    range_levels <- levels(plane_data$range_type)
    neck_position <- which(range_levels == "neck")
    segment_positions <- which(range_levels != "neck")
    
    neck_contrasts <- map(
      segment_positions,
      function(segment_position) {
        contrast_vector <- rep(0, length(range_levels))
        contrast_vector[neck_position] <- 1
        contrast_vector[segment_position] <- -1
        contrast_vector
      }
    )
    names(neck_contrasts) <- paste("neck -", range_levels[segment_positions])
    
    fit <- lmer(
      movement_range ~ range_type + (1 | beetle_id),
      data = plane_data,
      REML = TRUE
    )
    
    emm <- emmeans(fit, ~ range_type)
    range_contrasts <- contrast(
      emm,
      method = neck_contrasts,
      adjust = "none"
    )
    
    combine_contrasts(range_contrasts) %>%
      mutate(
        rotation_plane = plane_label,
        n_obs = nobs(fit),
        n_beetles = n_distinct(fit@frame$beetle_id)
      )
  }
) %>%
  select(
    rotation_plane,
    contrast,
    estimate,
    SE,
    df,
    lower.CL,
    upper.CL,
    t.ratio,
    p_raw,
    p_bonferroni,
    n_obs,
    n_beetles
  )

cat("\n===== NECK VERSUS BODY-SEGMENT MOVEMENT RANGE — MODEL CONTRASTS =====\n")
cat("Movement ranges are compared at the beetle level within each rotation plane.\n")
cat("Positive estimates mean that neck movement range is larger than the body-segment range in the contrast.\n")
cat("Both raw and Bonferroni-adjusted p-values are shown.\n\n")
print(
  as_tibble(segment_neck_range_contrasts %>%
              mutate(
                across(
                  c(estimate, SE, df, lower.CL, upper.CL, t.ratio),
                  ~ round(.x, 3)
                ),
                p_raw = format_p(p_raw),
                p_bonferroni = format_p(p_bonferroni)
              )),
  n = Inf
)

# centered neck angle distribution
p_neck_centered <- ggplot(
  neck_centered,
  aes(
    x = rotation_plane,
    y = neck_centered
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_violin(fill = "grey85", colour = NA, alpha = 0.8, trim = TRUE) +
  geom_jitter(width = 0.12, alpha = 0.20, size = 0.8) +
  geom_point(
    data = beetle_neck,
    aes(y = min_neck),
    colour = "#0072B2",
    alpha = 0.8,
    size = 1.8,
    position = position_jitter(width = 0.10, height = 0)
  ) +
  geom_point(
    data = beetle_neck,
    aes(y = max_neck),
    colour = "#D55E00",
    alpha = 0.8,
    size = 1.8,
    position = position_jitter(width = 0.10, height = 0)
  ) +
  labs(
    title = "Centered neck angles during flat-horizon dances",
    subtitle = "Grey = frame-level values; blue/orange = each beetle's minimum and maximum centered neck angle",
    x = "Rotation plane",
    y = "Centered neck angle (deg)"
  ) +
  theme_classic()

# beetle-level neck movement range
neck_summary_plot <- neck_summary %>%
  mutate(
    statistic = factor(statistic, levels = c("mean", "median"))
  )

p_neck_range <- ggplot() +
  geom_jitter(
    data = beetle_neck,
    aes(
      x = rotation_plane,
      y = neck_range
    ),
    width = 0.12,
    alpha = 0.45,
    size = 1.8
  ) +
  geom_pointrange(
    data = neck_summary_plot,
    aes(
      x = rotation_plane,
      y = estimate,
      ymin = ci_low,
      ymax = ci_high,
      shape = statistic
    ),
    position = position_dodge(width = 0.35),
    size = 0.7,
    linewidth = 0.8
  ) +
  labs(
    title = "Within-beetle neck movement range",
    subtitle = "Points are individual beetles; intervals are bootstrap 95% CIs across beetles",
    x = "Rotation plane",
    y = "Neck movement range (deg)",
    shape = "Summary"
  ) +
  theme_classic()

# save plots
ggsave(
  filename = "raw_figures/main-2-centered-neck-angles.svg",
  plot = p_neck_centered,
  width = 7,
  height = 5,
  device = svglite,
  bg = "white"
)

ggsave(
  filename = "raw_figures/main-2-neck-movement-range.svg",
  plot = p_neck_range,
  width = 7,
  height = 5,
  device = svglite,
  bg = "white"
)

print(p_neck_centered)
print(p_neck_range)

# clean section 2 objects
rm(
  neck_data,
  neck_centered,
  beetle_neck,
  neck_n,
  neck_summary,
  segment_neck_body_ranges,
  segment_neck_ranges,
  segment_neck_range_contrasts,
  neck_summary_plot,
  p_neck_centered,
  p_neck_range
)


# 3 - Active head stabilisation by neck movement -------------------------

# flat horizon with both antennae
compensation_data <- analysis_data %>%
  filter(
    condition == "flat0",
    antennae == "both",
    !is.na(head_angle_clean),
    !is.na(thorax_angle_clean)
  ) %>%
  mutate(
    neck_angle = head_angle_clean - thorax_angle_clean
  ) %>%
  group_by(rotation_plane, beetle_id) %>%
  mutate(
    thorax_dev = thorax_angle_clean - mean(thorax_angle_clean),
    head_dev = head_angle_clean - mean(head_angle_clean),
    neck_dev = neck_angle - mean(neck_angle)
  ) %>%
  ungroup()

# extract slope estimate and compensation value
slope_summary <- function(data, response, model_label, model_version, adjusted = FALSE) {
  model_formula <- if (adjusted) {
    as.formula(paste0(response, " ~ thorax_dev + ", external_terms, " + (1 | beetle_id)"))
  } else {
    as.formula(paste0(response, " ~ thorax_dev + (1 | beetle_id)"))
  }
  
  fit <- lmer(
    model_formula,
    data = data,
    REML = TRUE
  )
  coefs <- coef(summary(fit))
  ci <- confint(fit, parm = "thorax_dev", method = "Wald")
  slope <- coefs["thorax_dev", "Estimate"]
  
  tibble(
    model_version = model_version,
    model = model_label,
    slope = slope,
    se = coefs["thorax_dev", "Std. Error"],
    df = coefs["thorax_dev", "df"],
    t_value = coefs["thorax_dev", "t value"],
    p_value = coefs["thorax_dev", "Pr(>|t|)"],
    ci_low = ci[1],
    ci_high = ci[2],
    compensation_percent = case_when(
      response == "head_dev" ~ (1 - slope) * 100,
      response == "neck_dev" ~ abs(slope) * 100
    ),
    n_obs = nobs(fit),
    n_beetles = n_distinct(fit@frame$beetle_id)
  )
}

# fit head and neck compensation models by rotation_plane
compensation_summary <- map_dfr(
  levels(compensation_data$rotation_plane),
  function(plane_label) {
    plane_data <- compensation_data %>%
      filter(rotation_plane == plane_label)
    
    bind_rows(
      slope_summary(
        data = plane_data,
        response = "head_dev",
        model_label = "head_vs_thorax",
        model_version = "unadjusted",
        adjusted = FALSE
      ),
      slope_summary(
        data = plane_data,
        response = "neck_dev",
        model_label = "neck_vs_thorax",
        model_version = "unadjusted",
        adjusted = FALSE
      ),
      slope_summary(
        data = plane_data,
        response = "head_dev",
        model_label = "head_vs_thorax",
        model_version = "external_adjusted",
        adjusted = TRUE
      ),
      slope_summary(
        data = plane_data,
        response = "neck_dev",
        model_label = "neck_vs_thorax",
        model_version = "external_adjusted",
        adjusted = TRUE
      )
    ) %>%
      mutate(rotation_plane = plane_label, .before = model_version)
  }
)

cat("\n===== ACTIVE HEAD STABILISATION — SLOPE SUMMARY =====\n")
cat("All variables are centered within each beetle and rotation plane.\n")
cat("Unadjusted and external-factor-adjusted models are both shown.\n")
cat("For head_vs_thorax, slope = 1 means no stabilisation and slope = 0 means perfect stabilisation.\n")
cat("For neck_vs_thorax, slope = 0 means no neck response and slope = -1 means perfect compensation.\n\n")
print(
  as_tibble(compensation_summary %>%
              mutate(
                across(
                  c(slope, se, df, t_value, ci_low, ci_high, compensation_percent),
                  ~ round(.x, 3)
                ),
                p_value = format_p(p_value)
              )),
  n = Inf
)

# model-slope lines for plots
compensation_plot_lines <- compensation_summary %>%
  filter(model_version == "external_adjusted") %>%
  select(rotation_plane, model, slope, ci_low, ci_high) %>%
  left_join(
    compensation_data %>%
      group_by(rotation_plane) %>%
      summarise(
        x_min = min(thorax_dev),
        x_max = max(thorax_dev),
        .groups = "drop"
      ),
    by = "rotation_plane"
  ) %>%
  rowwise() %>%
  mutate(thorax_dev = list(seq(x_min, x_max, length.out = 100))) %>%
  unnest(thorax_dev) %>%
  ungroup() %>%
  mutate(
    fit = slope * thorax_dev,
    ci_a = ci_low * thorax_dev,
    ci_b = ci_high * thorax_dev,
    ci_lower = pmin(ci_a, ci_b),
    ci_upper = pmax(ci_a, ci_b)
  )

head_plot_lines <- compensation_plot_lines %>%
  filter(model == "head_vs_thorax")

neck_plot_lines <- compensation_plot_lines %>%
  filter(model == "neck_vs_thorax")

# head stabilisation plot
p_head_stabilisation <- ggplot(
  compensation_data,
  aes(
    x = thorax_dev,
    y = head_dev
  )
) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#D55E00") +
  geom_point(alpha = 0.25, size = 1.0) +
  geom_ribbon(
    data = head_plot_lines,
    aes(
      x = thorax_dev,
      ymin = ci_lower,
      ymax = ci_upper
    ),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.35
  ) +
  geom_line(
    data = head_plot_lines,
    aes(
      x = thorax_dev,
      y = fit
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.9
  ) +
  facet_wrap(~ rotation_plane) +
  coord_equal() +
  labs(
    title = "Does the head follow thorax movement?",
    subtitle = "Variables are centered within each beetle; solid line = adjusted mixed-model slope; dashed line = no stabilisation",
    x = "Thorax deviation from beetle mean (deg)",
    y = "Head deviation from beetle mean (deg)"
  ) +
  theme_classic()

# neck compensation plot
p_neck_compensation <- ggplot(
  compensation_data,
  aes(
    x = thorax_dev,
    y = neck_dev
  )
) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed", colour = "#D55E00") +
  geom_point(alpha = 0.25, size = 1.0) +
  geom_ribbon(
    data = neck_plot_lines,
    aes(
      x = thorax_dev,
      ymin = ci_lower,
      ymax = ci_upper
    ),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.35
  ) +
  geom_line(
    data = neck_plot_lines,
    aes(
      x = thorax_dev,
      y = fit
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.9
  ) +
  facet_wrap(~ rotation_plane) +
  coord_equal() +
  labs(
    title = "Does the neck compensate thorax movement?",
    subtitle = "Variables are centered within each beetle; solid line = adjusted mixed-model slope; dashed line = perfect compensation",
    x = "Thorax deviation from beetle mean (deg)",
    y = "Neck deviation from beetle mean (deg)"
  ) +
  theme_classic()

# save plots
ggsave(
  filename = "raw_figures/main-3-head-stabilisation.svg",
  plot = p_head_stabilisation,
  width = 8,
  height = 5,
  device = svglite,
  bg = "white"
)

ggsave(
  filename = "raw_figures/main-3-neck-compensation.svg",
  plot = p_neck_compensation,
  width = 8,
  height = 5,
  device = svglite,
  bg = "white"
)

print(p_head_stabilisation)
print(p_neck_compensation)

# clean section 3 objects
rm(
  compensation_data,
  compensation_summary,
  compensation_plot_lines,
  head_plot_lines,
  neck_plot_lines,
  p_head_stabilisation,
  p_neck_compensation,
  slope_summary
)


# 4 - Horizon effect on head angle ---------------------------------------

# both antennae only, raw corrected head angle
horizon_data <- analysis_data %>%
  filter(
    antennae == "both",
    !is.na(head_angle_clean)
  ) %>%
  mutate(
    horizon_context = case_when(
      view == "front" ~ "roll",
      view %in% c("sideL", "sideR") ~ "pitch"
    ),
    horizon_group = case_when(
      condition == "flat0" & view == "front" ~ "flat0_roll",
      condition == "tilted10" & view == "front" ~ "tilted10_roll",
      condition == "flat0" & view %in% c("sideL", "sideR") ~ "flat0_pitch",
      condition == "tilted10" & view == "sideL" ~ "tilted10_pitchL",
      condition == "tilted10" & view == "sideR" ~ "tilted10_pitchR"
    ),
    horizon_context = factor(horizon_context, levels = c("roll", "pitch")),
    horizon_group = factor(
      horizon_group,
      levels = c(
        "flat0_roll",
        "tilted10_roll",
        "flat0_pitch",
        "tilted10_pitchL",
        "tilted10_pitchR"
      )
    )
  ) %>%
  filter(!is.na(horizon_group))

# group descriptive statistics based on beetle means
horizon_beetle_means <- horizon_data %>%
  group_by(horizon_context, horizon_group, beetle_id) %>%
  summarise(
    beetle_mean = mean(head_angle_clean),
    .groups = "drop"
  )

horizon_n <- horizon_data %>%
  group_by(horizon_context, horizon_group) %>%
  summarise(
    n_obs = n(),
    n_beetles = n_distinct(beetle_id),
    .groups = "drop"
  )

horizon_descriptive <- horizon_beetle_means %>%
  group_by(horizon_context, horizon_group) %>%
  group_modify(~ boot_summary(.x$beetle_mean)) %>%
  ungroup() %>%
  left_join(horizon_n, by = c("horizon_context", "horizon_group")) %>%
  select(
    horizon_context,
    horizon_group,
    statistic,
    estimate,
    ci_low,
    ci_high,
    n_obs,
    n_beetles
  )

cat("\n===== HORIZON EFFECT — DESCRIPTIVE HEAD ANGLES =====\n")
cat("Values are beetle-level mean head angles with bootstrap 95% CIs across beetles.\n\n")
print(
  as_tibble(horizon_descriptive %>%
              mutate(
                estimate = round(estimate, 2),
                ci_low = round(ci_low, 2),
                ci_high = round(ci_high, 2)
              )),
  n = Inf
)

# roll model: flat horizon vs tilted horizon in front view
horizon_roll <- horizon_data %>%
  filter(horizon_context == "roll") %>%
  mutate(
    horizon_group = factor(
      horizon_group,
      levels = c("flat0_roll", "tilted10_roll")
    )
  )

mod_horizon_roll <- lmer(
  head_angle_clean ~ horizon_group + (1 | beetle_id),
  data = horizon_roll,
  REML = TRUE
)

mod_horizon_roll_adjusted <- lmer(
  as.formula(paste0("head_angle_clean ~ horizon_group + ", external_terms, " + (1 | beetle_id)")),
  data = horizon_roll,
  REML = TRUE
)

emm_horizon_roll <- emmeans(mod_horizon_roll, ~ horizon_group)
emm_horizon_roll_adjusted <- emmeans(mod_horizon_roll_adjusted, ~ horizon_group)

contrast_horizon_roll <- contrast(
  emm_horizon_roll,
  method = list("tilted10_roll - flat0_roll" = c(-1, 1)),
  adjust = "none"
)

contrast_horizon_roll_adjusted <- contrast(
  emm_horizon_roll_adjusted,
  method = list("tilted10_roll - flat0_roll" = c(-1, 1)),
  adjust = "none"
)

# pitch model: flat sides merged, tilted sides kept separate
horizon_pitch <- horizon_data %>%
  filter(horizon_context == "pitch") %>%
  mutate(
    horizon_group = factor(
      horizon_group,
      levels = c("flat0_pitch", "tilted10_pitchL", "tilted10_pitchR")
    )
  )

mod_horizon_pitch <- lmer(
  head_angle_clean ~ horizon_group + (1 | beetle_id),
  data = horizon_pitch,
  REML = TRUE
)

mod_horizon_pitch_adjusted <- lmer(
  as.formula(paste0("head_angle_clean ~ horizon_group + ", external_terms, " + (1 | beetle_id)")),
  data = horizon_pitch,
  REML = TRUE
)

emm_horizon_pitch <- emmeans(mod_horizon_pitch, ~ horizon_group)
emm_horizon_pitch_adjusted <- emmeans(mod_horizon_pitch_adjusted, ~ horizon_group)

contrast_horizon_pitch <- contrast(
  emm_horizon_pitch,
  method = list(
    "tilted10_pitchL - flat0_pitch" = c(-1, 1, 0),
    "tilted10_pitchR - flat0_pitch" = c(-1, 0, 1)
  ),
  adjust = "none"
)

contrast_horizon_pitch_adjusted <- contrast(
  emm_horizon_pitch_adjusted,
  method = list(
    "tilted10_pitchL - flat0_pitch" = c(-1, 1, 0),
    "tilted10_pitchR - flat0_pitch" = c(-1, 0, 1)
  ),
  adjust = "none"
)

horizon_emmeans <- bind_rows(
  as.data.frame(confint(emm_horizon_roll_adjusted)) %>%
    mutate(horizon_context = "roll"),
  as.data.frame(confint(emm_horizon_pitch_adjusted)) %>%
    mutate(horizon_context = "pitch")
) %>%
  rename(
    ci_low = lower.CL,
    ci_high = upper.CL
  ) %>%
  mutate(
    horizon_context = factor(horizon_context, levels = c("roll", "pitch")),
    horizon_group = factor(
      horizon_group,
      levels = levels(horizon_data$horizon_group)
    )
  )

horizon_contrasts <- bind_rows(
  combine_contrasts(contrast_horizon_roll) %>%
    mutate(horizon_context = "roll", model_version = "unadjusted"),
  combine_contrasts(contrast_horizon_pitch) %>%
    mutate(horizon_context = "pitch", model_version = "unadjusted"),
  combine_contrasts(contrast_horizon_roll_adjusted) %>%
    mutate(horizon_context = "roll", model_version = "external_adjusted"),
  combine_contrasts(contrast_horizon_pitch_adjusted) %>%
    mutate(horizon_context = "pitch", model_version = "external_adjusted")
) %>%
  select(
    model_version,
    horizon_context,
    contrast,
    estimate,
    SE,
    df,
    lower.CL,
    upper.CL,
    t.ratio,
    p_raw,
    p_bonferroni
  )

cat("\n===== HORIZON EFFECT — MODEL CONTRASTS =====\n")
cat("Positive estimates mean that the tilted horizon group has a higher corrected head angle than the flat reference.\n")
cat("Both unadjusted and external-factor-adjusted models are shown.\n")
cat("Both raw and Bonferroni-adjusted p-values are shown.\n\n")
print_contrasts(horizon_contrasts)

# plot labels
group_labels <- c(
  "flat0_roll" = "Flat",
  "tilted10_roll" = "Tilted",
  "flat0_pitch" = "Flat",
  "tilted10_pitchL" = "Tilted L",
  "tilted10_pitchR" = "Tilted R"
)

# horizon effect plot
p_horizon_head <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_violin(
    data = horizon_data,
    aes(
      x = horizon_group,
      y = head_angle_clean
    ),
    fill = "grey85",
    colour = NA,
    alpha = 0.8,
    trim = TRUE
  ) +
  geom_jitter(
    data = horizon_data,
    aes(
      x = horizon_group,
      y = head_angle_clean
    ),
    width = 0.12,
    alpha = 0.18,
    size = 0.8
  ) +
  geom_jitter(
    data = horizon_beetle_means,
    aes(
      x = horizon_group,
      y = beetle_mean
    ),
    width = 0.12,
    alpha = 0.85,
    size = 2.2,
    colour = "#0072B2"
  ) +
  geom_pointrange(
    data = horizon_emmeans,
    aes(
      x = horizon_group,
      y = emmean,
      ymin = ci_low,
      ymax = ci_high
    ),
    colour = "grey25",
    size = 0.75,
    linewidth = 0.8
  ) +
  facet_wrap(~ horizon_context, scales = "free_x") +
  scale_x_discrete(labels = group_labels) +
  labs(
    title = "Effect of horizon tilt on head angle",
    subtitle = "Grey = frames; blue = beetle means; dark grey = adjusted model EMM ± 95% CI",
    x = "Horizon condition",
    y = "Corrected head angle (deg)"
  ) +
  theme_classic()

# save plot
ggsave(
  filename = "raw_figures/main-4-horizon-head-angle.svg",
  plot = p_horizon_head,
  width = 9,
  height = 5,
  device = svglite,
  bg = "white"
)

print(p_horizon_head)

# clean section 4 objects
rm(
  horizon_data,
  horizon_beetle_means,
  horizon_n,
  horizon_descriptive,
  horizon_roll,
  horizon_pitch,
  mod_horizon_roll,
  mod_horizon_roll_adjusted,
  mod_horizon_pitch,
  mod_horizon_pitch_adjusted,
  emm_horizon_roll,
  emm_horizon_roll_adjusted,
  emm_horizon_pitch,
  emm_horizon_pitch_adjusted,
  contrast_horizon_roll,
  contrast_horizon_roll_adjusted,
  contrast_horizon_pitch,
  contrast_horizon_pitch_adjusted,
  horizon_emmeans,
  horizon_contrasts,
  group_labels,
  p_horizon_head
)


# 5 - Antennae effect on head angle --------------------------------------

# flat horizon only, raw corrected head angle
antennae_data <- analysis_data %>%
  filter(
    condition == "flat0",
    !is.na(head_angle_clean)
  ) %>%
  mutate(
    antennae_sep = factor(antennae, levels = c("both", "right", "left", "none")),
    antennae_grp = case_when(
      antennae == "both" ~ "both",
      antennae %in% c("right", "left") ~ "one",
      antennae == "none" ~ "none"
    ),
    antennae_grp = factor(antennae_grp, levels = c("both", "one", "none"))
  )

# descriptive statistics with left and right antennae separated
antennae_sep_beetle_means <- antennae_data %>%
  group_by(rotation_plane, antennae_sep, beetle_id) %>%
  summarise(
    beetle_mean = mean(head_angle_clean),
    .groups = "drop"
  )

antennae_sep_n <- antennae_data %>%
  group_by(rotation_plane, antennae_sep) %>%
  summarise(
    n_obs = n(),
    n_beetles = n_distinct(beetle_id),
    .groups = "drop"
  )

antennae_sep_descriptive <- antennae_sep_beetle_means %>%
  group_by(rotation_plane, antennae_sep) %>%
  group_modify(~ boot_summary(.x$beetle_mean)) %>%
  ungroup() %>%
  left_join(antennae_sep_n, by = c("rotation_plane", "antennae_sep")) %>%
  select(
    rotation_plane,
    antennae_sep,
    statistic,
    estimate,
    ci_low,
    ci_high,
    n_obs,
    n_beetles
  )

cat("\n===== ANTENNAE EFFECT — DESCRIPTIVE HEAD ANGLES, LEFT/RIGHT SEPARATED =====\n")
cat("Values are beetle-level mean head angles with bootstrap 95% CIs across beetles.\n\n")
print(
  as_tibble(antennae_sep_descriptive %>%
              mutate(
                estimate = round(estimate, 2),
                ci_low = round(ci_low, 2),
                ci_high = round(ci_high, 2)
              )),
  n = Inf
)

# descriptive statistics with left and right antennae merged
antennae_grp_beetle_means <- antennae_data %>%
  group_by(rotation_plane, antennae_grp, beetle_id) %>%
  summarise(
    beetle_mean = mean(head_angle_clean),
    .groups = "drop"
  )

antennae_grp_n <- antennae_data %>%
  group_by(rotation_plane, antennae_grp) %>%
  summarise(
    n_obs = n(),
    n_beetles = n_distinct(beetle_id),
    .groups = "drop"
  )

antennae_grp_descriptive <- antennae_grp_beetle_means %>%
  group_by(rotation_plane, antennae_grp) %>%
  group_modify(~ boot_summary(.x$beetle_mean)) %>%
  ungroup() %>%
  left_join(antennae_grp_n, by = c("rotation_plane", "antennae_grp")) %>%
  select(
    rotation_plane,
    antennae_grp,
    statistic,
    estimate,
    ci_low,
    ci_high,
    n_obs,
    n_beetles
  )

cat("\n===== ANTENNAE EFFECT — DESCRIPTIVE HEAD ANGLES, SINGLE ANTENNA MERGED =====\n")
cat("Values are beetle-level mean head angles with bootstrap 95% CIs across beetles.\n\n")
print(
  as_tibble(antennae_grp_descriptive %>%
              mutate(
                estimate = round(estimate, 2),
                ci_low = round(ci_low, 2),
                ci_high = round(ci_high, 2)
              )),
  n = Inf
)

# model helper for antennae contrasts
antennae_model_summary <- function(data, group_var, adjusted = FALSE) {
  model_formula <- if (adjusted) {
    as.formula(paste0("head_angle_clean ~ ", group_var, " + ", external_terms, " + (1 | beetle_id)"))
  } else {
    as.formula(paste0("head_angle_clean ~ ", group_var, " + (1 | beetle_id)"))
  }
  
  fit <- lmer(
    model_formula,
    data = data,
    REML = TRUE
  )
  emm <- emmeans(fit, as.formula(paste0("~ ", group_var)))
  contrasts <- contrast(
    emm,
    method = "pairwise",
    adjust = "none"
  )
  list(
    fit = fit,
    emm = emm,
    contrasts = contrasts
  )
}

# models with left and right antennae separated
antennae_sep_results <- map(
  levels(antennae_data$rotation_plane),
  function(plane_label) {
    plane_data <- antennae_data %>%
      filter(rotation_plane == plane_label)
    list(
      unadjusted = antennae_model_summary(
        data = plane_data,
        group_var = "antennae_sep",
        adjusted = FALSE
      ),
      external_adjusted = antennae_model_summary(
        data = plane_data,
        group_var = "antennae_sep",
        adjusted = TRUE
      )
    )
  }
)
names(antennae_sep_results) <- levels(antennae_data$rotation_plane)

antennae_sep_emmeans <- map_dfr(
  names(antennae_sep_results),
  function(plane_label) {
    as.data.frame(confint(antennae_sep_results[[plane_label]]$external_adjusted$emm)) %>%
      mutate(rotation_plane = plane_label)
  }
) %>%
  rename(
    ci_low = lower.CL,
    ci_high = upper.CL
  ) %>%
  mutate(
    rotation_plane = factor(rotation_plane, levels = levels(antennae_data$rotation_plane)),
    antennae_sep = factor(antennae_sep, levels = levels(antennae_data$antennae_sep))
  )

antennae_sep_contrasts <- map_dfr(
  names(antennae_sep_results),
  function(plane_label) {
    bind_rows(
      combine_contrasts(antennae_sep_results[[plane_label]]$unadjusted$contrasts) %>%
        mutate(model_version = "unadjusted"),
      combine_contrasts(antennae_sep_results[[plane_label]]$external_adjusted$contrasts) %>%
        mutate(model_version = "external_adjusted")
    ) %>%
      mutate(rotation_plane = plane_label)
  }
) %>%
  select(
    model_version,
    rotation_plane,
    contrast,
    estimate,
    SE,
    df,
    lower.CL,
    upper.CL,
    t.ratio,
    p_raw,
    p_bonferroni
  )

cat("\n===== ANTENNAE EFFECT — MODEL CONTRASTS, LEFT/RIGHT SEPARATED =====\n")
cat("Contrasts are all pairwise comparisons within each rotation plane.\n")
cat("Both unadjusted and external-factor-adjusted models are shown.\n")
cat("Both raw and Bonferroni-adjusted p-values are shown.\n\n")
print_contrasts(antennae_sep_contrasts)

# models with left and right antennae merged
antennae_grp_results <- map(
  levels(antennae_data$rotation_plane),
  function(plane_label) {
    plane_data <- antennae_data %>%
      filter(rotation_plane == plane_label)
    list(
      unadjusted = antennae_model_summary(
        data = plane_data,
        group_var = "antennae_grp",
        adjusted = FALSE
      ),
      external_adjusted = antennae_model_summary(
        data = plane_data,
        group_var = "antennae_grp",
        adjusted = TRUE
      )
    )
  }
)
names(antennae_grp_results) <- levels(antennae_data$rotation_plane)

antennae_grp_emmeans <- map_dfr(
  names(antennae_grp_results),
  function(plane_label) {
    as.data.frame(confint(antennae_grp_results[[plane_label]]$external_adjusted$emm)) %>%
      mutate(rotation_plane = plane_label)
  }
) %>%
  rename(
    ci_low = lower.CL,
    ci_high = upper.CL
  ) %>%
  mutate(
    rotation_plane = factor(rotation_plane, levels = levels(antennae_data$rotation_plane)),
    antennae_grp = factor(antennae_grp, levels = levels(antennae_data$antennae_grp))
  )

antennae_grp_contrasts <- map_dfr(
  names(antennae_grp_results),
  function(plane_label) {
    bind_rows(
      combine_contrasts(antennae_grp_results[[plane_label]]$unadjusted$contrasts) %>%
        mutate(model_version = "unadjusted"),
      combine_contrasts(antennae_grp_results[[plane_label]]$external_adjusted$contrasts) %>%
        mutate(model_version = "external_adjusted")
    ) %>%
      mutate(rotation_plane = plane_label)
  }
) %>%
  select(
    model_version,
    rotation_plane,
    contrast,
    estimate,
    SE,
    df,
    lower.CL,
    upper.CL,
    t.ratio,
    p_raw,
    p_bonferroni
  )

cat("\n===== ANTENNAE EFFECT — MODEL CONTRASTS, SINGLE ANTENNA MERGED =====\n")
cat("Contrasts are all pairwise comparisons within each rotation plane.\n")
cat("Both unadjusted and external-factor-adjusted models are shown.\n")
cat("Both raw and Bonferroni-adjusted p-values are shown.\n\n")
print_contrasts(antennae_grp_contrasts)

# plot labels
antennae_sep_labels <- c(
  "both" = "Both",
  "right" = "Right only",
  "left" = "Left only",
  "none" = "None"
)

antennae_grp_labels <- c(
  "both" = "Both",
  "one" = "One",
  "none" = "None"
)

# antennae effect plot with left and right separated
p_antennae_sep <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_violin(
    data = antennae_data,
    aes(
      x = antennae_sep,
      y = head_angle_clean
    ),
    fill = "grey85",
    colour = NA,
    alpha = 0.8,
    trim = TRUE
  ) +
  geom_jitter(
    data = antennae_data,
    aes(
      x = antennae_sep,
      y = head_angle_clean
    ),
    width = 0.12,
    alpha = 0.18,
    size = 0.8
  ) +
  geom_jitter(
    data = antennae_sep_beetle_means,
    aes(
      x = antennae_sep,
      y = beetle_mean
    ),
    width = 0.12,
    alpha = 0.85,
    size = 2.2,
    colour = "#0072B2"
  ) +
  geom_pointrange(
    data = antennae_sep_emmeans,
    aes(
      x = antennae_sep,
      y = emmean,
      ymin = ci_low,
      ymax = ci_high
    ),
    colour = "grey25",
    size = 0.75,
    linewidth = 0.8
  ) +
  facet_wrap(~ rotation_plane) +
  scale_x_discrete(labels = antennae_sep_labels) +
  labs(
    title = "Effect of antennae condition on head angle",
    subtitle = "Single-antenna treatments separated; grey = frames; blue = beetle means; dark grey = adjusted model EMM ± 95% CI",
    x = "Antennae condition",
    y = "Corrected head angle (deg)"
  ) +
  theme_classic()

# antennae effect plot with left and right merged
p_antennae_grp <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_violin(
    data = antennae_data,
    aes(
      x = antennae_grp,
      y = head_angle_clean
    ),
    fill = "grey85",
    colour = NA,
    alpha = 0.8,
    trim = TRUE
  ) +
  geom_jitter(
    data = antennae_data,
    aes(
      x = antennae_grp,
      y = head_angle_clean
    ),
    width = 0.12,
    alpha = 0.18,
    size = 0.8
  ) +
  geom_jitter(
    data = antennae_grp_beetle_means,
    aes(
      x = antennae_grp,
      y = beetle_mean
    ),
    width = 0.12,
    alpha = 0.85,
    size = 2.2,
    colour = "#0072B2"
  ) +
  geom_pointrange(
    data = antennae_grp_emmeans,
    aes(
      x = antennae_grp,
      y = emmean,
      ymin = ci_low,
      ymax = ci_high
    ),
    colour = "grey25",
    size = 0.75,
    linewidth = 0.8
  ) +
  facet_wrap(~ rotation_plane) +
  scale_x_discrete(labels = antennae_grp_labels) +
  labs(
    title = "Effect of antennae condition on head angle",
    subtitle = "Single-antenna treatments merged; grey = frames; blue = beetle means; dark grey = adjusted model EMM ± 95% CI",
    x = "Antennae condition",
    y = "Corrected head angle (deg)"
  ) +
  theme_classic()

# save plots
ggsave(
  filename = "raw_figures/main-5-antennae-head-angle-separated.svg",
  plot = p_antennae_sep,
  width = 9,
  height = 5,
  device = svglite,
  bg = "white"
)

ggsave(
  filename = "raw_figures/main-5-antennae-head-angle-merged.svg",
  plot = p_antennae_grp,
  width = 9,
  height = 5,
  device = svglite,
  bg = "white"
)

print(p_antennae_sep)
print(p_antennae_grp)

# clean section 5 objects
rm(
  antennae_data,
  antennae_sep_beetle_means,
  antennae_sep_n,
  antennae_sep_descriptive,
  antennae_grp_beetle_means,
  antennae_grp_n,
  antennae_grp_descriptive,
  antennae_model_summary,
  antennae_sep_results,
  antennae_sep_emmeans,
  antennae_sep_contrasts,
  antennae_grp_results,
  antennae_grp_emmeans,
  antennae_grp_contrasts,
  antennae_sep_labels,
  antennae_grp_labels,
  p_antennae_sep,
  p_antennae_grp
)


# 6 - Combined horizon and antennae effects on head angle -----------------

# all antennae conditions, single antennae merged
combined_data <- analysis_data %>%
  filter(!is.na(head_angle_clean)) %>%
  mutate(
    antennae_grp = case_when(
      antennae == "both" ~ "both",
      antennae %in% c("right", "left") ~ "one",
      antennae == "none" ~ "none"
    ),
    antennae_grp = factor(antennae_grp, levels = c("both", "one", "none"))
  )

# contrast print helper for section 6
format_contrast_table <- function(raw_table) {
  raw_table %>%
    mutate(
      p_bonferroni = p.adjust(p_raw, method = "bonferroni")
    ) %>%
    mutate(
      across(
        c(estimate, SE, df, lower.CL, upper.CL, t.ratio),
        ~ round(.x, 3)
      ),
      p_raw = format_p(p_raw),
      p_bonferroni = format_p(p_bonferroni)
    )
}

# roll interaction
combined_roll <- combined_data %>%
  filter(view == "front") %>%
  mutate(
    condition = factor(condition, levels = c("flat0", "tilted10"))
  )

combined_roll_descriptive <- combined_roll %>%
  group_by(condition, antennae_grp, beetle_id) %>%
  summarise(
    beetle_mean = mean(head_angle_clean),
    .groups = "drop"
  ) %>%
  group_by(condition, antennae_grp) %>%
  group_modify(~ boot_summary(.x$beetle_mean)) %>%
  ungroup()

combined_roll_n <- combined_roll %>%
  group_by(condition, antennae_grp) %>%
  summarise(
    n_obs = n(),
    n_beetles = n_distinct(beetle_id),
    .groups = "drop"
  )

combined_roll_descriptive <- combined_roll_descriptive %>%
  left_join(combined_roll_n, by = c("condition", "antennae_grp"))

cat("\n===== COMBINED HORIZON AND ANTENNAE EFFECTS — ROLL DESCRIPTIVE HEAD ANGLES =====\n")
cat("Values are beetle-level mean head angles with bootstrap 95% CIs across beetles.\n\n")
print(
  as_tibble(combined_roll_descriptive %>%
              mutate(
                estimate = round(estimate, 2),
                ci_low = round(ci_low, 2),
                ci_high = round(ci_high, 2)
              )),
  n = Inf
)

mod_combined_roll <- lmer(
  head_angle_clean ~ condition * antennae_grp + (1 | beetle_id),
  data = combined_roll,
  REML = TRUE
)

mod_combined_roll_adjusted <- lmer(
  as.formula(paste0("head_angle_clean ~ condition * antennae_grp + ", external_terms, " + (1 | beetle_id)")),
  data = combined_roll,
  REML = TRUE
)

cat("\n===== COMBINED HORIZON AND ANTENNAE EFFECTS — ROLL MODEL, UNADJUSTED =====\n")
print(anova(mod_combined_roll, type = 3))

cat("\n===== COMBINED HORIZON AND ANTENNAE EFFECTS — ROLL MODEL, EXTERNAL-FACTOR ADJUSTED =====\n")
print(anova(mod_combined_roll_adjusted, type = 3))

emm_combined_roll <- emmeans(
  mod_combined_roll,
  ~ condition | antennae_grp
)

emm_combined_roll_adjusted <- emmeans(
  mod_combined_roll_adjusted,
  ~ condition | antennae_grp
)

contrast_combined_roll <- contrast(
  emm_combined_roll,
  method = list("tilted10 - flat0" = c(-1, 1)),
  adjust = "none"
)

contrast_combined_roll_adjusted <- contrast(
  emm_combined_roll_adjusted,
  method = list("tilted10 - flat0" = c(-1, 1)),
  adjust = "none"
)

combined_roll_contrasts <- bind_rows(
  combine_contrasts(contrast_combined_roll, group_cols = "antennae_grp") %>%
    mutate(model_version = "unadjusted"),
  combine_contrasts(contrast_combined_roll_adjusted, group_cols = "antennae_grp") %>%
    mutate(model_version = "external_adjusted")
) %>%
  select(
    model_version,
    antennae_grp,
    contrast,
    estimate,
    SE,
    df,
    lower.CL,
    upper.CL,
    t.ratio,
    p_raw,
    p_bonferroni
  )

cat("\n===== COMBINED HORIZON AND ANTENNAE EFFECTS — ROLL HORIZON CONTRASTS =====\n")
cat("Contrasts show the tilted-minus-flat horizon effect within each antennae condition.\n")
cat("Both unadjusted and external-factor-adjusted models are shown.\n")
cat("Both raw and Bonferroni-adjusted p-values are shown.\n\n")
print_contrasts(combined_roll_contrasts)

# compare antennae groups within each roll horizon condition
emm_combined_roll_antennae <- emmeans(
  mod_combined_roll,
  ~ antennae_grp | condition
)

emm_combined_roll_antennae_adjusted <- emmeans(
  mod_combined_roll_adjusted,
  ~ antennae_grp | condition
)

contrast_combined_roll_antennae <- contrast(
  emm_combined_roll_antennae,
  method = "pairwise",
  adjust = "none"
)

contrast_combined_roll_antennae_adjusted <- contrast(
  emm_combined_roll_antennae_adjusted,
  method = "pairwise",
  adjust = "none"
)

combined_roll_antennae_contrasts <- bind_rows(
  combine_contrasts(contrast_combined_roll_antennae, group_cols = "condition") %>%
    mutate(model_version = "unadjusted"),
  combine_contrasts(contrast_combined_roll_antennae_adjusted, group_cols = "condition") %>%
    mutate(model_version = "external_adjusted")
) %>%
  select(
    model_version,
    condition,
    contrast,
    estimate,
    SE,
    df,
    lower.CL,
    upper.CL,
    t.ratio,
    p_raw,
    p_bonferroni
  )

cat("\n===== COMBINED HORIZON AND ANTENNAE EFFECTS — ROLL ANTENNAE CONTRASTS WITHIN HORIZON CONDITION =====\n")
cat("Contrasts compare antennae groups within each roll horizon condition.\n")
cat("Negative estimates mean that the first antennae group has a lower corrected head angle than the second.\n")
cat("Both unadjusted and external-factor-adjusted models are shown.\n")
cat("Both raw and Bonferroni-adjusted p-values are shown.\n\n")
print_contrasts(combined_roll_antennae_contrasts)

combined_roll_emmeans <- as.data.frame(confint(emmeans(mod_combined_roll_adjusted, ~ condition * antennae_grp))) %>%
  rename(
    ci_low = lower.CL,
    ci_high = upper.CL
  )

# pitch interaction
combined_pitch <- combined_data %>%
  filter(view %in% c("sideL", "sideR")) %>%
  mutate(
    horizon_pitch = case_when(
      condition == "flat0" ~ "flat0_pitch",
      condition == "tilted10" & view == "sideL" ~ "tilted10_pitchL",
      condition == "tilted10" & view == "sideR" ~ "tilted10_pitchR"
    ),
    horizon_pitch = factor(
      horizon_pitch,
      levels = c("flat0_pitch", "tilted10_pitchL", "tilted10_pitchR")
    )
  )

combined_pitch_descriptive <- combined_pitch %>%
  group_by(horizon_pitch, antennae_grp, beetle_id) %>%
  summarise(
    beetle_mean = mean(head_angle_clean),
    .groups = "drop"
  ) %>%
  group_by(horizon_pitch, antennae_grp) %>%
  group_modify(~ boot_summary(.x$beetle_mean)) %>%
  ungroup()

combined_pitch_n <- combined_pitch %>%
  group_by(horizon_pitch, antennae_grp) %>%
  summarise(
    n_obs = n(),
    n_beetles = n_distinct(beetle_id),
    .groups = "drop"
  )

combined_pitch_descriptive <- combined_pitch_descriptive %>%
  left_join(combined_pitch_n, by = c("horizon_pitch", "antennae_grp"))

cat("\n===== COMBINED HORIZON AND ANTENNAE EFFECTS — PITCH DESCRIPTIVE HEAD ANGLES =====\n")
cat("Values are beetle-level mean head angles with bootstrap 95% CIs across beetles.\n\n")
print(
  as_tibble(combined_pitch_descriptive %>%
              mutate(
                estimate = round(estimate, 2),
                ci_low = round(ci_low, 2),
                ci_high = round(ci_high, 2)
              )),
  n = Inf
)

mod_combined_pitch <- lmer(
  head_angle_clean ~ horizon_pitch * antennae_grp + (1 | beetle_id),
  data = combined_pitch,
  REML = TRUE
)

mod_combined_pitch_adjusted <- lmer(
  as.formula(paste0("head_angle_clean ~ horizon_pitch * antennae_grp + ", external_terms, " + (1 | beetle_id)")),
  data = combined_pitch,
  REML = TRUE
)

cat("\n===== COMBINED HORIZON AND ANTENNAE EFFECTS — PITCH MODEL, UNADJUSTED =====\n")
print(anova(mod_combined_pitch, type = 3))

cat("\n===== COMBINED HORIZON AND ANTENNAE EFFECTS — PITCH MODEL, EXTERNAL-FACTOR ADJUSTED =====\n")
print(anova(mod_combined_pitch_adjusted, type = 3))

emm_combined_pitch <- emmeans(
  mod_combined_pitch,
  ~ horizon_pitch | antennae_grp
)

emm_combined_pitch_adjusted <- emmeans(
  mod_combined_pitch_adjusted,
  ~ horizon_pitch | antennae_grp
)

contrast_combined_pitch <- contrast(
  emm_combined_pitch,
  method = list(
    "tilted10_pitchL - flat0_pitch" = c(-1, 1, 0),
    "tilted10_pitchR - flat0_pitch" = c(-1, 0, 1)
  ),
  adjust = "none"
)

contrast_combined_pitch_adjusted <- contrast(
  emm_combined_pitch_adjusted,
  method = list(
    "tilted10_pitchL - flat0_pitch" = c(-1, 1, 0),
    "tilted10_pitchR - flat0_pitch" = c(-1, 0, 1)
  ),
  adjust = "none"
)

combined_pitch_contrasts <- bind_rows(
  combine_contrasts(contrast_combined_pitch, group_cols = "antennae_grp") %>%
    mutate(model_version = "unadjusted"),
  combine_contrasts(contrast_combined_pitch_adjusted, group_cols = "antennae_grp") %>%
    mutate(model_version = "external_adjusted")
) %>%
  select(
    model_version,
    antennae_grp,
    contrast,
    estimate,
    SE,
    df,
    lower.CL,
    upper.CL,
    t.ratio,
    p_raw,
    p_bonferroni
  )

cat("\n===== COMBINED HORIZON AND ANTENNAE EFFECTS — PITCH HORIZON CONTRASTS =====\n")
cat("Contrasts show the tilted-minus-flat horizon effect within each antennae condition.\n")
cat("Both unadjusted and external-factor-adjusted models are shown.\n")
cat("Both raw and Bonferroni-adjusted p-values are shown.\n\n")
print_contrasts(combined_pitch_contrasts)

# compare antennae groups within each pitch horizon condition
emm_combined_pitch_antennae <- emmeans(
  mod_combined_pitch,
  ~ antennae_grp | horizon_pitch
)

emm_combined_pitch_antennae_adjusted <- emmeans(
  mod_combined_pitch_adjusted,
  ~ antennae_grp | horizon_pitch
)

contrast_combined_pitch_antennae <- contrast(
  emm_combined_pitch_antennae,
  method = "pairwise",
  adjust = "none"
)

contrast_combined_pitch_antennae_adjusted <- contrast(
  emm_combined_pitch_antennae_adjusted,
  method = "pairwise",
  adjust = "none"
)

combined_pitch_antennae_contrasts <- bind_rows(
  combine_contrasts(contrast_combined_pitch_antennae, group_cols = "horizon_pitch") %>%
    mutate(model_version = "unadjusted"),
  combine_contrasts(contrast_combined_pitch_antennae_adjusted, group_cols = "horizon_pitch") %>%
    mutate(model_version = "external_adjusted")
) %>%
  select(
    model_version,
    horizon_pitch,
    contrast,
    estimate,
    SE,
    df,
    lower.CL,
    upper.CL,
    t.ratio,
    p_raw,
    p_bonferroni
  )

cat("\n===== COMBINED HORIZON AND ANTENNAE EFFECTS — PITCH ANTENNAE CONTRASTS WITHIN HORIZON CONDITION =====\n")
cat("Contrasts compare antennae groups within each pitch horizon condition.\n")
cat("Negative estimates mean that the first antennae group has a lower corrected head angle than the second.\n")
cat("Both unadjusted and external-factor-adjusted models are shown.\n")
cat("Both raw and Bonferroni-adjusted p-values are shown.\n\n")
print_contrasts(combined_pitch_antennae_contrasts)

combined_pitch_emmeans <- as.data.frame(confint(emmeans(mod_combined_pitch_adjusted, ~ horizon_pitch * antennae_grp))) %>%
  rename(
    ci_low = lower.CL,
    ci_high = upper.CL
  )

# plots
combined_antennae_labels <- c(
  "both" = "Both",
  "one" = "One",
  "none" = "None"
)

combined_antennae_colours <- c(
  "both" = "#1B9E77",
  "one" = "#D95F02",
  "none" = "#7570B3"
)

combined_roll_condition_labels <- c(
  "flat0" = "Flat",
  "tilted10" = "Tilted"
)

combined_pitch_horizon_labels <- c(
  "flat0_pitch" = "Flat",
  "tilted10_pitchL" = "Tilted L",
  "tilted10_pitchR" = "Tilted R"
)

combined_roll_beetles <- combined_roll %>%
  group_by(condition, antennae_grp, beetle_id) %>%
  summarise(beetle_mean = mean(head_angle_clean), .groups = "drop")

combined_pitch_beetles <- combined_pitch %>%
  group_by(horizon_pitch, antennae_grp, beetle_id) %>%
  summarise(beetle_mean = mean(head_angle_clean), .groups = "drop")

p_combined_roll <- ggplot(
  combined_roll,
  aes(
    x = condition,
    y = head_angle_clean,
    fill = antennae_grp,
    colour = antennae_grp
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_violin(
    position = position_dodge(width = 0.75),
    width = 0.65,
    alpha = 0.20,
    colour = NA,
    trim = TRUE
  ) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.75),
    alpha = 0.12,
    size = 0.7,
    show.legend = FALSE
  ) +
  geom_point(
    data = combined_roll_beetles,
    aes(
      x = condition,
      y = beetle_mean,
      fill = antennae_grp
    ),
    position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.75),
    shape = 21,
    colour = "grey35",
    stroke = 0.25,
    alpha = 0.45,
    size = 2.2
  ) +
  geom_pointrange(
    data = combined_roll_emmeans,
    aes(
      x = condition,
      y = emmean,
      ymin = ci_low,
      ymax = ci_high
    ),
    position = position_dodge(width = 0.75),
    colour = "grey25",
    linewidth = 0.8,
    size = 0.75
  ) +
  scale_x_discrete(labels = combined_roll_condition_labels) +
  scale_colour_manual(values = combined_antennae_colours, labels = combined_antennae_labels) +
  scale_fill_manual(values = combined_antennae_colours, labels = combined_antennae_labels) +
  labs(
    title = "Combined horizon and antennae effects on head angle — roll",
    subtitle = "Violin = frame distribution; pale dots = beetle means; dark grey intervals = adjusted model EMM ± 95% CI",
    x = "Horizon condition",
    y = "Corrected head angle (deg)",
    colour = "Antennae",
    fill = "Antennae"
  ) +
  theme_classic()

p_combined_pitch <- ggplot(
  combined_pitch,
  aes(
    x = horizon_pitch,
    y = head_angle_clean,
    fill = antennae_grp,
    colour = antennae_grp
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_violin(
    position = position_dodge(width = 0.75),
    width = 0.65,
    alpha = 0.20,
    colour = NA,
    trim = TRUE
  ) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.75),
    alpha = 0.12,
    size = 0.7,
    show.legend = FALSE
  ) +
  geom_point(
    data = combined_pitch_beetles,
    aes(
      x = horizon_pitch,
      y = beetle_mean,
      fill = antennae_grp
    ),
    position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.75),
    shape = 21,
    colour = "grey35",
    stroke = 0.25,
    alpha = 0.45,
    size = 2.2
  ) +
  geom_pointrange(
    data = combined_pitch_emmeans,
    aes(
      x = horizon_pitch,
      y = emmean,
      ymin = ci_low,
      ymax = ci_high
    ),
    position = position_dodge(width = 0.75),
    colour = "grey25",
    linewidth = 0.8,
    size = 0.75
  ) +
  scale_x_discrete(labels = combined_pitch_horizon_labels) +
  scale_colour_manual(values = combined_antennae_colours, labels = combined_antennae_labels) +
  scale_fill_manual(values = combined_antennae_colours, labels = combined_antennae_labels) +
  labs(
    title = "Combined horizon and antennae effects on head angle — pitch",
    subtitle = "Violin = frame distribution; pale dots = beetle means; dark grey intervals = adjusted model EMM ± 95% CI",
    x = "Horizon condition",
    y = "Corrected head angle (deg)",
    colour = "Antennae",
    fill = "Antennae"
  ) +
  theme_classic()

# save plots
ggsave(
  filename = "raw_figures/main-6-combined-horizon-antennae-roll.svg",
  plot = p_combined_roll,
  width = 8,
  height = 5,
  device = svglite,
  bg = "white"
)

ggsave(
  filename = "raw_figures/main-6-combined-horizon-antennae-pitch.svg",
  plot = p_combined_pitch,
  width = 9,
  height = 5,
  device = svglite,
  bg = "white"
)

print(p_combined_roll)
print(p_combined_pitch)

# clean section 6 objects
rm(
  combined_data,
  format_contrast_table,
  combined_roll,
  combined_roll_descriptive,
  combined_roll_n,
  mod_combined_roll,
  mod_combined_roll_adjusted,
  emm_combined_roll,
  emm_combined_roll_adjusted,
  contrast_combined_roll,
  contrast_combined_roll_adjusted,
  combined_roll_contrasts,
  emm_combined_roll_antennae,
  emm_combined_roll_antennae_adjusted,
  contrast_combined_roll_antennae,
  contrast_combined_roll_antennae_adjusted,
  combined_roll_antennae_contrasts,
  combined_roll_emmeans,
  combined_pitch,
  combined_pitch_descriptive,
  combined_pitch_n,
  mod_combined_pitch,
  mod_combined_pitch_adjusted,
  emm_combined_pitch,
  emm_combined_pitch_adjusted,
  contrast_combined_pitch,
  contrast_combined_pitch_adjusted,
  combined_pitch_contrasts,
  emm_combined_pitch_antennae,
  emm_combined_pitch_antennae_adjusted,
  contrast_combined_pitch_antennae,
  contrast_combined_pitch_antennae_adjusted,
  combined_pitch_antennae_contrasts,
  combined_pitch_emmeans,
  combined_antennae_labels,
  combined_antennae_colours,
  combined_roll_condition_labels,
  combined_pitch_horizon_labels,
  combined_roll_beetles,
  combined_pitch_beetles,
  p_combined_roll,
  p_combined_pitch
)
