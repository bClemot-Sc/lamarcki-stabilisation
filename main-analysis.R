# Head stabilisation in Kheper lamarcki - main analysis

library(tidyverse)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(emmeans)
library(svglite)
library(gghalves)

# Check if dhglm and car are installed (functions directly called later)
for (pkg in c("dhglm", "car")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required. Install it with install.packages('", pkg, "') and rerun the script.")
  }
}
if (sink.number() > 0) stop("An output sink is already active. Restart R and run the script again.")

set.seed(42)

# Setting the factor contrast
options(contrasts = c("contr.sum", "contr.poly"))

SCRIPT_VERSION <- "v26-10000bootstraps"
data_file <- "data_stabilisation_processed.csv"
output_file <- "reviews/main_analysis_review_v24.txt"
n_boot <- 10000

dir.create("raw_figures", showWarnings = FALSE)
dir.create("reviews", showWarnings = FALSE)

# 1. Data ------------------------------------------------------------------

# Remove a value if both outlier rules flag it, or one rule flags it and the row has a comment
flag_outlier <- function(z3, iqr, has_comment) {
  z3 <- coalesce(z3, FALSE)
  iqr <- coalesce(iqr, FALSE)
  (z3 & iqr) | ((z3 | iqr) & has_comment)
}

# Import and prepare data
analysis_data <- read.csv(data_file, fileEncoding = "UTF-8") %>%
  mutate(
    condition = factor(condition, levels = c("flat0", "tilted10")),
    antennae = factor(antennae, levels = c("both", "right", "left", "none")),
    view = factor(recode(view, "left" = "sideL", "right" = "sideR"),
                  levels = c("front", "sideL", "sideR")),
    rotation_plane = factor(
      case_when(view == "front" ~ "roll", view %in% c("sideL", "sideR") ~ "pitch"),
      levels = c("roll", "pitch")
    ),
    has_comment = !is.na(comments) & str_trim(comments) != "",
    outlier_remove_head = flag_outlier(outlier_z3_head, outlier_iqr_head, has_comment),
    outlier_remove_thorax = flag_outlier(outlier_z3_thorax, outlier_iqr_thorax, has_comment),
    head_angle_clean = if_else(outlier_remove_head, NA_real_, head_angle_corrected),
    thorax_angle_clean = if_else(outlier_remove_thorax, NA_real_, thorax_angle_corrected),
    antennae_grp = factor(
      case_when(
        antennae == "both" ~ "both",
        antennae %in% c("left", "right") ~ "one",
        antennae == "none" ~ "none"
      ),
      levels = c("both", "one", "none")
    ),
    # pitch: sideL = low tilt, sideR = high tilt
    horizon_pitch = factor(
      case_when(
        rotation_plane == "pitch" & condition == "flat0" ~ "flat0_pitch",
        rotation_plane == "pitch" & condition == "tilted10" & view == "sideL" ~ "tilted10_pitchL",
        rotation_plane == "pitch" & condition == "tilted10" & view == "sideR" ~ "tilted10_pitchR"
      ),
      levels = c("flat0_pitch", "tilted10_pitchL", "tilted10_pitchR")
    )
  ) %>%
  filter(is.na(comments) | !str_detect(comments, regex("discard", ignore_case = TRUE)))

cat("\nDATA SUMMARY\n")
cat("Rows retained after explicit discard comments:", nrow(analysis_data), "\n")
cat("Unique beetles:", n_distinct(analysis_data$beetle_id), "\n")
cat("Head values removed as outliers:", sum(analysis_data$outlier_remove_head, na.rm = TRUE), "\n")
cat("Thorax values removed as outliers:", sum(analysis_data$outlier_remove_thorax, na.rm = TRUE), "\n")

# 2. Helpers ---------------------------------------------------------------

# Format p values
format_p <- function(p) ifelse(is.na(p), NA_character_, formatC(p, format = "e", digits = 2))
# Add Holm correct p values column
holm <- function(x) if (is.null(x) || nrow(x) == 0) x else mutate(x, p_holm = p.adjust(p_raw, method = "holm"))
# Data output functions
save_svg <- function(plot, name, width = 8, height = 5) {
  ggsave(paste0("raw_figures/", name, ".svg"), plot,
         width = width, height = height, device = svglite, bg = "white")
}
# Round values for output
print_rounded <- function(x, digits = 3) {
  print.data.frame(as.data.frame(mutate(x, across(where(is.numeric), ~ round(.x, digits)))), row.names = FALSE)
}

# Rename emmeans CI for cleaner output
ci_names <- function(s) {
  c(intersect(c("lower.CL", "asymp.LCL"), names(s))[1],
    intersect(c("upper.CL", "asymp.UCL"), names(s))[1])
}
add_ci <- function(s) {
  nm <- ci_names(s)
  mutate(s, lower.CL = .data[[nm[1]]], upper.CL = .data[[nm[2]]])
}

# Check if a glmmTMB model has converged and fitted properly
glmmTMB_ok <- function(model, label) {
  ok <- model$fit$convergence == 0 && isTRUE(model$sdr$pdHess)
  if (!ok) warning(label, " convergence/Hessian problem")
  ok
}

# Extract p value from Anova for output
anova_p <- function(tab, effect) {
  tab <- as.data.frame(tab)
  p_col <- grep("^Pr\\(", names(tab), value = TRUE)[1]
  if (is.na(p_col) || !(effect %in% rownames(tab))) return(NA_real_)
  as.numeric(tab[effect, p_col])
}

# Adjust which horizon contrats are tested based on test families and calculate contrats with emmeans
horizon_contrasts <- function(emm, horizon_var) {
  methods <- if (horizon_var == "condition") {
    list("tilted - flat" = c(-1, 1))
  } else {
    list("low - flat" = c(-1, 1, 0), "high - flat" = c(-1, 0, 1), "low - high" = c(0, 1, -1))
  }
  contrast(emm, method = methods, adjust = "none")
}

# Re-organise emmeans contrats output 
tidy_contrast <- function(x) {
  s <- as.data.frame(summary(x, infer = c(TRUE, TRUE), adjust = "none", type = "link"))
  nm <- ci_names(s)
  stat <- intersect(c("t.ratio", "z.ratio"), names(s))[1]
  if (anyNA(c(nm, stat))) stop("Could not identify CI/statistic columns in emmeans output.")
  s %>%
    mutate(lower.CL = .data[[nm[1]]], upper.CL = .data[[nm[2]]],
           test_statistic = .data[[stat]], p_raw = p.value) %>%
    select(any_of(c("contrast", "condition", "horizon_pitch", "antennae_grp")),
           estimate, SE, any_of("df"), lower.CL, upper.CL, test_statistic, p_raw)
}

# Output contrast tables in the review file
print_contrasts <- function(x) {
  if (is.null(x) || nrow(x) == 0) {
    cat("No post-hoc contrasts in this family.\n")
    return(invisible(NULL))
  }
  rounded <- c("estimate", "SE", "df", "lower.CL", "upper.CL", "test_statistic",
               "variance_ratio", "variance_ratio_low", "variance_ratio_high",
               "sd_ratio", "sd_ratio_low", "sd_ratio_high")
  x %>%
    mutate(across(any_of(rounded), ~ round(.x, 3)),
           across(any_of(c("p_raw", "p_holm")), format_p)) %>%
    as.data.frame() %>%
    print.data.frame(row.names = FALSE)
}

# 3. Q1 - baseline head-thorax slope ---------------------------------------
# Baseline = flat horizon + both antennae. Thorax is centred within beetle.

baseline_paired <- analysis_data %>%
  filter(condition == "flat0", antennae == "both",
         !is.na(head_angle_clean), !is.na(thorax_angle_clean), !is.na(rotation_plane))

# Centre data
q1_data <- baseline_paired %>%
  group_by(rotation_plane, beetle_id) %>%
  mutate(thorax_dev = thorax_angle_clean - mean(thorax_angle_clean),
         head_dev = head_angle_clean - mean(head_angle_clean)) %>%
  ungroup()

# Define model (will be tested for roll and pitch), fitted models are kept for the diagnostics
q1_models <- list()
fit_q1 <- function(plane) {
  d <- filter(q1_data, rotation_plane == plane)
  # === The model fitted here
  fit <- lmer(head_angle_clean ~ thorax_dev + (1 | beetle_id), data = d, REML = TRUE)
  # ===
  q1_models[[plane]] <<- fit
  cf <- coef(summary(fit))["thorax_dev", ]
  ci <- confint(fit, parm = "thorax_dev", method = "Wald")
  t_vs_1 <- (cf[["Estimate"]] - 1) / cf[["Std. Error"]]
  tibble(
    rotation_plane = plane,
    slope = cf[["Estimate"]], SE = cf[["Std. Error"]], df = cf[["df"]],
    ci_low = ci[1], ci_high = ci[2],
    t_vs_1 = t_vs_1, p_vs_1 = 2 * pt(-abs(t_vs_1), df = cf[["df"]]),
    n_paired_frames = nobs(fit), n_beetles = n_distinct(d$beetle_id)
  )
}

# Resample whole beetles and refit the same model
cluster_bootstrap <- function(plane, n_boot, seed) {
  set.seed(seed)
  d <- filter(q1_data, rotation_plane == plane)
  ids <- unique(d$beetle_id)
  slopes <- rep(NA_real_, n_boot)
  
  for (b in seq_len(n_boot)) {
    sampled_ids <- sample(ids, size = length(ids), replace = TRUE)
    db <- map2_dfr(sampled_ids, seq_along(sampled_ids), function(id, draw_id) {
      d %>% filter(beetle_id == id) %>% mutate(boot_beetle = factor(paste0(as.character(id), "_", draw_id)))
    })
    fit <- try(suppressWarnings(lmer(head_angle_clean ~ thorax_dev + (1 | boot_beetle), data = db, REML = TRUE)),
               silent = TRUE)
    if (!inherits(fit, "try-error")) slopes[b] <- fixef(fit)["thorax_dev"]
  }
  
  good <- slopes[is.finite(slopes)]
  if (length(good) < 0.9 * n_boot) warning("More than 10% of Q1 bootstrap fits failed for ", plane)
  
  tibble(
    rotation_plane = plane, n_boot_requested = n_boot, n_boot_success = length(good),
    bootstrap_se = sd(good),
    bootstrap_ci_low = unname(quantile(good, 0.025)),
    bootstrap_ci_high = unname(quantile(good, 0.975))
  )
}

# Table output for the model
q1_summary <- bind_rows(fit_q1("roll"), fit_q1("pitch")) %>%
  left_join(bind_rows(cluster_bootstrap("roll", n_boot, seed = 42),
                      cluster_bootstrap("pitch", n_boot, seed = 43)),
            by = "rotation_plane")

# Plot data
q1_lines <- q1_summary %>%
  left_join(q1_data %>% group_by(rotation_plane) %>%
              summarise(x_min = min(thorax_dev), x_max = max(thorax_dev), .groups = "drop"),
            by = "rotation_plane") %>%
  rowwise() %>%
  mutate(thorax_dev = list(seq(x_min, x_max, length.out = 100))) %>%
  unnest(thorax_dev) %>%
  ungroup() %>%
  mutate(fit = slope * thorax_dev,
         ribbon_low = pmin(bootstrap_ci_low * thorax_dev, bootstrap_ci_high * thorax_dev),
         ribbon_high = pmax(bootstrap_ci_low * thorax_dev, bootstrap_ci_high * thorax_dev))

q1_labels <- q1_summary %>%
  mutate(label = paste0(if_else(rotation_plane == "roll", "Roll", "Pitch"),
                        ": slope = ", sprintf("%.2f", slope)))
p_q1 <- ggplot(q1_data, aes(thorax_dev, head_dev, colour = rotation_plane, shape = rotation_plane)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey85") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.6, colour = "grey45") +
  geom_ribbon(data = q1_lines, aes(x = thorax_dev, ymin = ribbon_low, ymax = ribbon_high, fill = rotation_plane),
              inherit.aes = FALSE, alpha = 0.12, colour = NA) +
  geom_point(alpha = 0.30, size = 1.4) +
  geom_line(data = q1_lines, aes(x = thorax_dev, y = fit, colour = rotation_plane),
            inherit.aes = FALSE, linewidth = 1.1) +
  geom_text(data = q1_labels,
            aes(x = Inf, y = Inf, label = label, colour = rotation_plane),
            inherit.aes = FALSE, hjust = 1.08, vjust = c(2.0, 3.4),
            fontface = "bold", size = 3.4, show.legend = FALSE) +
  scale_colour_manual(values = c(roll = "#0072B2", pitch = "#D55E00"),
                      labels = c(roll = "Roll", pitch = "Pitch"), name = NULL) +
  scale_fill_manual(values = c(roll = "#0072B2", pitch = "#D55E00"), guide = "none") +
  scale_shape_manual(values = c(roll = 16, pitch = 17),
                     labels = c(roll = "Roll", pitch = "Pitch"), name = NULL) +
  guides(colour = guide_legend(override.aes = list(shape = c(16, 17), linewidth = 0, size = 2.5)),
         shape = "none") +
  coord_equal() +
  labs(title = "Q1: Baseline head-thorax coupling",
       subtitle = "Solid = model slope; shading = 95% bootstrap CI; dashed = 1:1",
       x = "Thorax deviation from beetle mean (°)",
       y = "Head deviation from beetle mean (°)") +
  theme_classic(base_size = 11) +
  theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.35),
        panel.grid.minor.y = element_blank()) +
  theme(legend.position = "top",
        legend.justification = "left",
        legend.margin = margin(0, 0, 2, 0))
save_svg(p_q1, "q1-baseline-head-thorax-slope")

# 4. Q2 - baseline head versus thorax SD -----------------------------------
# Paired within-beetle SDs; one-sample t-test of log(head SD / thorax SD) against 0.

# Get and prepare data, need at least 3 frames to estimate
q2_sd_by_beetle <- baseline_paired %>%
  group_by(rotation_plane, beetle_id) %>%
  summarise(n_frames = n(), head_sd = sd(head_angle_clean), thorax_sd = sd(thorax_angle_clean),
            .groups = "drop") %>%
  filter(n_frames >= 3, is.finite(head_sd), is.finite(thorax_sd), head_sd > 0, thorax_sd > 0) %>%
  mutate(log_sd_ratio = log(head_sd / thorax_sd))

# test for SD difference
fit_q2 <- function(plane) {
  d <- filter(q2_sd_by_beetle, rotation_plane == plane)
  if (nrow(d) < 2) stop("Q2 baseline SD comparison: fewer than two beetles available for ", plane)
  test <- t.test(d$log_sd_ratio, mu = 0)
  tibble(
    rotation_plane = plane,
    n_beetles = nrow(d),
    n_matched_frames = sum(d$n_frames),
    head_geometric_mean_sd = exp(mean(log(d$head_sd))),
    thorax_geometric_mean_sd = exp(mean(log(d$thorax_sd))),
    sd_ratio_head_over_thorax = exp(mean(d$log_sd_ratio)),
    sd_ratio_low = exp(test$conf.int[1]),
    sd_ratio_high = exp(test$conf.int[2]),
    t_value = unname(test$statistic),
    df = unname(test$parameter),
    p_value = test$p.value
  )
}

# Table Output
q2_summary <- bind_rows(fit_q2("roll"), fit_q2("pitch"))

# First plot for the within beetles SD
q2_sd_labels <- q2_summary %>%
  select(rotation_plane, head_geometric_mean_sd, thorax_geometric_mean_sd) %>%
  pivot_longer(c(head_geometric_mean_sd, thorax_geometric_mean_sd),
               names_to = "body_segment", values_to = "sd_value") %>%
  mutate(body_segment = recode(body_segment,
                               head_geometric_mean_sd = "Head",
                               thorax_geometric_mean_sd = "Thorax"),
         x = if_else(body_segment == "Head", 0.85, 1.15),
         plot_group = paste(rotation_plane, tolower(body_segment), sep = "_"),
         label = paste0("SD = ", sprintf("%.2f", sd_value), "°"))

q2_sd_plot_data <- q2_sd_by_beetle %>%
  select(rotation_plane, beetle_id, head_sd, thorax_sd) %>%
  pivot_longer(c(head_sd, thorax_sd), names_to = "body_segment", values_to = "movement_sd") %>%
  mutate(body_segment = recode(body_segment, head_sd = "Head", thorax_sd = "Thorax"),
         x = if_else(body_segment == "Head", 0.85, 1.15),
         plot_group = paste(rotation_plane, tolower(body_segment), sep = "_"))

p_q2_sd <- ggplot(q2_sd_plot_data, aes(y = movement_sd)) +
  geom_hline(yintercept = seq(0, 25, 5), linewidth = 0.3, colour = "grey92") +
  gghalves::geom_half_violin(data = filter(q2_sd_plot_data, body_segment == "Head"),
                             aes(x = 0.65, fill = plot_group, colour = plot_group),
                             side = "l", alpha = 0.28, trim = FALSE, width = 0.35, linewidth = 0.55) +
  gghalves::geom_half_violin(data = filter(q2_sd_plot_data, body_segment == "Thorax"),
                             aes(x = 1.35, fill = plot_group, colour = plot_group),
                             side = "r", alpha = 0.28, trim = FALSE, width = 0.35, linewidth = 0.55) +
  geom_boxplot(data = filter(q2_sd_plot_data, body_segment == "Head"),
               aes(x = 0.74, colour = plot_group), width = 0.07,
               outlier.shape = NA, fill = "white", linewidth = 0.5) +
  geom_boxplot(data = filter(q2_sd_plot_data, body_segment == "Thorax"),
               aes(x = 1.26, colour = plot_group), width = 0.07,
               outlier.shape = NA, fill = "white", linewidth = 0.5) +
  geom_line(aes(x = x, group = beetle_id, colour = rotation_plane),
            alpha = 0.24, linewidth = 0.4) +
  geom_point(aes(x = x, colour = plot_group), size = 1.9, alpha = 0.85) +
  geom_text(data = q2_sd_labels, aes(x = x, y = Inf, label = label, colour = plot_group),
            inherit.aes = FALSE, vjust = 1.4, fontface = "bold", size = 3.2, show.legend = FALSE) +
  facet_wrap(~ rotation_plane, labeller = as_labeller(c(roll = "Roll", pitch = "Pitch")),
             axes = "all_y", axis.labels = "all_y") +
  scale_colour_manual(values = c(roll_head = "#8CC8E8", roll_thorax = "#0066A6",
                                 pitch_head = "#F3A074", pitch_thorax = "#C94B00",
                                 roll = "#0066A6", pitch = "#C94B00"), guide = "none") +
  scale_fill_manual(values = c(roll_head = "#8CC8E8", roll_thorax = "#0066A6",
                               pitch_head = "#F3A074", pitch_thorax = "#C94B00"), guide = "none") +
  scale_x_continuous(breaks = c(0.85, 1.15), labels = c("Head", "Thorax"),
                     limits = c(0.45, 1.55), expand = c(0, 0)) +
  labs(title = "Q2: Baseline angular variability",
       subtitle = "Connected points = paired within-beetle SDs",
       x = NULL, y = "Within-beetle SD (°)") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        panel.spacing.x = unit(0.9, "lines"))

save_svg(p_q2_sd, "q2-baseline-head-thorax-sd")

# Second plot for the distribution of head and thorax angles
q2_centered_plot_data <- baseline_paired %>%
  group_by(rotation_plane, beetle_id) %>%
  mutate(head_centered = head_angle_clean - mean(head_angle_clean),
         thorax_centered = thorax_angle_clean - mean(thorax_angle_clean)) %>%
  ungroup() %>%
  select(rotation_plane, beetle_id, head_centered, thorax_centered) %>%
  pivot_longer(c(head_centered, thorax_centered),
               names_to = "body_segment", values_to = "angle_centered") %>%
  mutate(body_segment = recode(body_segment,
                               head_centered = "Head",
                               thorax_centered = "Thorax"),
         x = if_else(body_segment == "Head", 0.85, 1.15),
         plot_group = paste(rotation_plane, tolower(body_segment), sep = "_"))

p_q2_centered <- ggplot(q2_centered_plot_data, aes(y = angle_centered)) +
  geom_hline(yintercept = seq(-60, 30, 30), linewidth = 0.3, colour = "grey92") +
  geom_hline(yintercept = 0, linewidth = 0.45, colour = "grey75") +
  gghalves::geom_half_violin(data = filter(q2_centered_plot_data, body_segment == "Head"),
                             aes(x = 0.65, fill = plot_group, colour = plot_group),
                             side = "l", alpha = 0.28, trim = FALSE, width = 0.35, linewidth = 0.55) +
  gghalves::geom_half_violin(data = filter(q2_centered_plot_data, body_segment == "Thorax"),
                             aes(x = 1.35, fill = plot_group, colour = plot_group),
                             side = "r", alpha = 0.28, trim = FALSE, width = 0.35, linewidth = 0.55) +
  geom_boxplot(data = filter(q2_centered_plot_data, body_segment == "Head"),
               aes(x = 0.74, colour = plot_group), width = 0.07,
               outlier.shape = NA, fill = "white", linewidth = 0.5) +
  geom_boxplot(data = filter(q2_centered_plot_data, body_segment == "Thorax"),
               aes(x = 1.26, colour = plot_group), width = 0.07,
               outlier.shape = NA, fill = "white", linewidth = 0.5) +
  geom_point(aes(x = x, colour = plot_group), size = 1.6, alpha = 0.55,
             position = position_jitter(width = 0.025, height = 0)) +
  geom_text(data = q2_sd_labels, aes(x = x, y = Inf, label = label, colour = plot_group),
            inherit.aes = FALSE, vjust = 1.4, fontface = "bold", size = 3.2, show.legend = FALSE) +
  facet_wrap(~ rotation_plane, labeller = as_labeller(c(roll = "Roll", pitch = "Pitch")),
             axes = "all_y", axis.labels = "all_y") +
  scale_colour_manual(values = c(roll_head = "#8CC8E8", roll_thorax = "#0066A6",
                                 pitch_head = "#F3A074", pitch_thorax = "#C94B00"), guide = "none") +
  scale_fill_manual(values = c(roll_head = "#8CC8E8", roll_thorax = "#0066A6",
                               pitch_head = "#F3A074", pitch_thorax = "#C94B00"), guide = "none") +
  scale_x_continuous(breaks = c(0.85, 1.15), labels = c("Head", "Thorax"),
                     limits = c(0.45, 1.55), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-60, 30, 30)) +
  labs(title = "Q2: Baseline centred angle distributions",
       subtitle = "Angles centred within beetle; labels show geometric mean within-beetle SD",
       x = NULL, y = "Deviation from individual mean (°)") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        panel.spacing.x = unit(0.9, "lines"))

save_svg(p_q2_centered, "q2-baseline-head-thorax-distributions")

# 5. Q3 + Q4 - sensory effects on mean orientation and spread --------------

# Prepare data for each axis
roll_data <- analysis_data %>%
  filter(rotation_plane == "roll", !is.na(antennae_grp), !is.na(condition)) %>%
  droplevels()
pitch_data <- analysis_data %>%
  filter(rotation_plane == "pitch", !is.na(antennae_grp), !is.na(horizon_pitch)) %>%
  droplevels()
planes <- list(roll = list(data = roll_data, horizon = "condition"),
               pitch = list(data = pitch_data, horizon = "horizon_pitch"))

# Response variables
segments <- c(head = "head_angle_clean", thorax = "thorax_angle_clean")
runs <- expand_grid(plane = names(planes), segment = names(segments))
run_key <- function(plane, segment) paste(plane, segment)

# Q3: Y ~ Horizon * Antennae + (1 | beetle), with cell-specific residual variance as a nuisance
fit_location <- function(data, response, horizon_var) {
  glmmTMB(
    formula = as.formula(paste0(response, " ~ ", horizon_var, " * antennae_grp + (1 | beetle_id)")),
    dispformula = as.formula(paste0("~ ", horizon_var, " * antennae_grp")),
    family = gaussian, data = data, REML = TRUE
  )
}

# type3 anova and planned post-hoc comparisons (simple and main effects) for the average position
run_q3 <- function(model, horizon_var) {
  emm_h <- emmeans(model, as.formula(paste0("~ ", horizon_var, " | antennae_grp")), component = "cond")
  emm_a <- emmeans(model, as.formula(paste0("~ antennae_grp | ", horizon_var)), component = "cond")
  # Pooled (main effect) contrasts, averaged over the other factor
  emm_hp <- emmeans(model, as.formula(paste0("~ ", horizon_var)), component = "cond")
  emm_ap <- emmeans(model, ~ antennae_grp, component = "cond")
  list(
    anova = car::Anova(model, type = 3, component = "cond"),
    mode = "planned_simple_effects",
    horizon = holm(tidy_contrast(horizon_contrasts(emm_h, horizon_var))),
    antenna = holm(tidy_contrast(pairs(emm_a, adjust = "none"))),
    horizon_pooled = holm(tidy_contrast(horizon_contrasts(emm_hp, horizon_var))),
    antenna_pooled = holm(tidy_contrast(pairs(emm_ap, adjust = "none")))
  )
}

# Extract model estimated mean angles and CI for each conditions
location_emmeans <- function(model, horizon_var) {
  emmeans(model, as.formula(paste0("~ ", horizon_var, " * antennae_grp")), component = "cond") %>%
    summary(infer = c(TRUE, TRUE), type = "response") %>%
    as.data.frame() %>%
    add_ci()
}

# Extract residual variance model(as backup, but will be inspected in q4)
scale_emmeans <- function(model, horizon_var) {
  s <- emmeans(model, as.formula(paste0("~ ", horizon_var, " * antennae_grp")), component = "disp") %>%
    summary(infer = c(TRUE, TRUE), type = "link") %>%
    as.data.frame()
  nm <- ci_names(s)
  mutate(s, residual_sd = exp(emmean / 2),
         sd_low = exp(.data[[nm[1]]] / 2), sd_high = exp(.data[[nm[2]]] / 2))
}

# Fit the 4 models (body segment x axis) and run analysis 
q3_models <- list()
q3 <- list()
for (i in seq_len(nrow(runs))) {
  pl <- runs$plane[i]; sg <- runs$segment[i]
  d <- planes[[pl]]$data
  q3_models[[run_key(pl, sg)]] <- fit_location(d[!is.na(d[[segments[[sg]]]]), ], segments[[sg]], planes[[pl]]$horizon)
}
q3_checks <- tibble(
  model = names(q3_models),
  converged = imap_lgl(q3_models, ~ glmmTMB_ok(.x, paste("Q3", .y)))
)
for (i in seq_len(nrow(runs))) {
  pl <- runs$plane[i]; sg <- runs$segment[i]
  q3[[run_key(pl, sg)]] <- run_q3(q3_models[[run_key(pl, sg)]], planes[[pl]]$horizon)
}
q3_emm <- map(names(planes), function(pl) {
  bind_rows(lapply(names(segments), function(sg) {
    location_emmeans(q3_models[[run_key(pl, sg)]], planes[[pl]]$horizon) %>%
      mutate(segment = sg, rotation_plane = pl)
  }))
}) %>% set_names(names(planes))

# Q4: DHGLM with a beetle random intercept in both location and scale.
# Scale model: log(phi) ~ Horizon * Antennae ("cell") vs Horizon + Antennae ("additive").

# Double hierarchical model, with the location being the same as Q3 but then with scale model added
fit_dhglm <- function(dat, scale_rhs) {
  mean_model <- dhglm::DHGLMMODELING(
    Model = "mean", Link = "identity",
    LinPred = y ~ horizon * antennae_grp + (1 | beetle_id), RandDist = "gaussian"
  )
  scale_model <- dhglm::DHGLMMODELING(
    Model = "dispersion", Link = "log",
    LinPred = as.formula(paste0("phi ~ ", scale_rhs, " + (1 | beetle_id)")), RandDist = "gaussian"
  )
  fit <- NULL
  capture.output(
    fit <- tryCatch(
      dhglm::dhglmfit(RespDist = "gaussian", DataMain = dat, MeanModel = mean_model,
                      DispersionModel = scale_model, REML = TRUE, Maxiter = 200, convergence = 1e-06),
      error = function(e) NULL
    )
  )
  fit
}

# Prepare data for each model to fit (data name and treatments cells)
fit_spread <- function(data, response, horizon_var, scale = "cell") {
  dat <- data %>%
    filter(!is.na(.data[[response]])) %>%
    transmute(y = .data[[response]], beetle_id = factor(beetle_id),
              antennae_grp = droplevels(antennae_grp), horizon = droplevels(.data[[horizon_var]])) %>%
    mutate(scale_cell = interaction(horizon, antennae_grp, sep = "__", drop = TRUE))
  rhs <- c(cell = "scale_cell", additive = "horizon + antennae_grp",
           horizon = "horizon", antennae = "antennae_grp")[[scale]]
  model <- fit_dhglm(dat, rhs)
  list(model = model, ok = !is.null(model), data = dat)
}

# Sanity checks for the models
q4_health <- function(fit) {
  m <- fit$model
  if (is.null(m)) {
    return(tibble(fit_returned = FALSE, finite_ml = FALSE, finite_coefficients = FALSE,
                  finite_scale_random_effect = FALSE, ml = NA_real_, rl = NA_real_, caic = NA_real_))
  }
  coefs <- unlist(lapply(list(m$beta_coeff, m$phi_coeff, m$lambda_coeff, m$tau_coeff),
                         function(x) if (is.null(x)) numeric(0) else as.numeric(x)))
  tibble(
    fit_returned = TRUE,
    finite_ml = is.finite(as.numeric(m$ml)),
    finite_coefficients = length(coefs) > 0 && all(is.finite(coefs)),
    finite_scale_random_effect = !is.null(m$tau_coeff) && all(is.finite(as.numeric(m$tau_coeff))),
    ml = as.numeric(m$ml), rl = as.numeric(m$rl), caic = as.numeric(m$caic)
  )
}

# Population residual SD per Horizon x Antennae cell from the full scale model, with Wald 95% CI.
# SE per cell: refit with that cell as reference level; the phi intercept is then the cell's log-variance.
q4_cell_sd <- function(fit, horizon_var) {
  if (!fit$ok) return(tibble())
  cells <- fit$data %>% distinct(horizon, antennae_grp, scale_cell) %>% arrange(horizon, antennae_grp)
  X <- model.matrix(~ scale_cell, data = cells)
  beta <- fit$model$phi_coeff[, "Estimate"]
  if (ncol(X) != length(beta)) stop("Q4 DHGLM: scale design matrix does not match phi coefficients.")
  log_variance <- as.numeric(X %*% beta)
  
  cell_se <- vapply(as.character(cells$scale_cell), function(cell) {
    d <- fit$data
    d$wald_term <- relevel(factor(as.character(d$scale_cell)), ref = cell)
    contrasts(d$wald_term) <- contr.treatment(nlevels(d$wald_term))
    m <- fit_dhglm(d, "wald_term")
    if (is.null(m)) NA_real_ else as.numeric(m$phi_coeff[1, "Std. Error"])
  }, numeric(1))
  
  cells %>%
    transmute(horizon = as.character(horizon), antennae_grp = as.character(antennae_grp),
              cell_key = as.character(scale_cell), log_variance = log_variance, log_variance_se = cell_se,
              residual_sd = exp(log_variance / 2),
              sd_low = exp((log_variance - 1.96 * cell_se) / 2),
              sd_high = exp((log_variance + 1.96 * cell_se) / 2)) %>%
    rename(!!horizon_var := horizon)
}

# Wald contrast of log-variance: refit with level_b as reference so level_a is a single coefficient
wald_contrast <- function(dat, var, level_a, level_b, scale_rhs) {
  dat$wald_term <- relevel(factor(as.character(dat[[var]])), ref = level_b)
  contrasts(dat$wald_term) <- contr.treatment(nlevels(dat$wald_term))
  fit <- fit_dhglm(dat, scale_rhs)
  if (is.null(fit)) {
    return(tibble(estimate = NA_real_, SE = NA_real_, z_value = NA_real_, p_raw = NA_real_,
                  sd_ratio = NA_real_, sd_ratio_low = NA_real_, sd_ratio_high = NA_real_))
  }
  cf <- as.data.frame(fit$phi_coeff)
  cf$term <- rownames(fit$phi_coeff)
  level <- which(levels(dat$wald_term) == level_a)
  row <- filter(cf, term == paste0("wald_term", level))
  if (nrow(row) != 1) row <- filter(cf, startsWith(term, "wald_term"))[level - 1L, , drop = FALSE]
  
  est <- row$Estimate[1]
  se <- row$`Std. Error`[1]
  z <- est / se
  tibble(estimate = est, SE = se, z_value = z,
         p_raw = 2 * pnorm(abs(z), lower.tail = FALSE),
         sd_ratio = exp(est / 2),
         sd_ratio_low = exp((est - 1.96 * se) / 2),
         sd_ratio_high = exp((est + 1.96 * se) / 2))
}

# Pooled Wald contrast: level_a - level_b of one factor averaged over the other (cell model with sum coding)
wald_pooled <- function(dat, var, level_a, level_b, other) {
  lv <- levels(dat[[var]])
  cvec <- ifelse(lv == level_a, 1, ifelse(lv == level_b, -1, 0))
  rest <- qr.Q(qr(cbind(1, cvec)), complete = TRUE)[, -(1:2), drop = FALSE]
  dat$wald_term <- dat[[var]]
  contrasts(dat$wald_term) <- unname(cbind(cvec, rest))
  dat$other_term <- dat[[other]]
  contrasts(dat$other_term) <- contr.sum(nlevels(dat$other_term))
  fit <- fit_dhglm(dat, "wald_term * other_term")
  if (is.null(fit)) {
    return(tibble(estimate = NA_real_, SE = NA_real_, z_value = NA_real_, p_raw = NA_real_,
                  sd_ratio = NA_real_, sd_ratio_low = NA_real_, sd_ratio_high = NA_real_))
  }
  cf <- as.data.frame(fit$phi_coeff)
  row <- cf[rownames(cf) == "wald_term1", , drop = FALSE]
  
  # The coefficient is half the level difference, so double it
  est <- 2 * row$Estimate[1]
  se <- 2 * row$`Std. Error`[1]
  z <- est / se
  tibble(estimate = est, SE = se, z_value = z,
         p_raw = 2 * pnorm(abs(z), lower.tail = FALSE),
         sd_ratio = exp(est / 2),
         sd_ratio_low = exp((est - 1.96 * se) / 2),
         sd_ratio_high = exp((est + 1.96 * se) / 2))
}

q4_pooled <- function(dat, var, level_a, level_b, other, label, family) {
  wald_pooled(dat, var, level_a, level_b, other) %>%
    mutate(contrast = label, family = family, test_method = "Wald") %>%
    select(contrast, family, estimate, SE, z_value, sd_ratio, sd_ratio_low, sd_ratio_high, p_raw, test_method)
}

# Cleaner output for the Wald contrast
q4_contrast <- function(full, cell_a, cell_b, label, family) {
  wald_contrast(full$data, "scale_cell", cell_a, cell_b, "wald_term") %>%
    mutate(contrast = label, family = family, test_method = "Wald") %>%
    select(contrast, family, estimate, SE, z_value, sd_ratio, sd_ratio_low, sd_ratio_high, p_raw, test_method)
}

# Fit full scale, then additive model, ANOVA, planned comparisons for the test and holm correction for stat family
run_q4 <- function(data, response, horizon_var) {
  full <- fit_spread(data, response, horizon_var, "cell")
  diagnostics <- q4_health(full) %>% mutate(model = "full scale") %>% select(model, everything())
  
  if (!full$ok) {
    return(list(full = full, cells = tibble(), horizon = tibble(), antenna = tibble(),
                horizon_pooled = tibble(), antenna_pooled = tibble(),
                mode = "fit_failed", diagnostics = diagnostics))
  }
  
  n_h <- nlevels(full$data$horizon)
  h_lv <- levels(full$data$horizon)
  a_lv <- levels(full$data$antennae_grp)
  h_pairs <- if (n_h == 2) list(h_lv[c(2, 1)]) else list(h_lv[c(2, 1)], h_lv[c(3, 1)], h_lv[c(2, 3)])
  a_pairs <- list(c("both", "one"), c("both", "none"), c("one", "none"))
  
  horizon <- bind_rows(lapply(a_lv, function(a) bind_rows(lapply(h_pairs, function(p) {
    q4_contrast(full, paste(p[1], a, sep = "__"), paste(p[2], a, sep = "__"),
                paste(p[1], "-", p[2], "| antennae =", a), "horizon")
  })))) %>% holm()
  antenna <- bind_rows(lapply(h_lv, function(h) bind_rows(lapply(a_pairs, function(p) {
    q4_contrast(full, paste(h, p[1], sep = "__"), paste(h, p[2], sep = "__"),
                paste(p[1], "-", p[2], "| horizon =", h), "antenna")
  })))) %>% holm()
  
  # Pooled contrasts, averaged over the other factor
  horizon_pooled <- bind_rows(lapply(h_pairs, function(p) {
    q4_pooled(full$data, "horizon", p[1], p[2], "antennae_grp",
              paste(p[1], "-", p[2], "| pooled over antennae"), "horizon")
  })) %>% holm()
  antenna_pooled <- bind_rows(lapply(a_pairs, function(p) {
    q4_pooled(full$data, "antennae_grp", p[1], p[2], "horizon",
              paste(p[1], "-", p[2], "| pooled over horizon"), "antenna")
  })) %>% holm()
  
  list(full = full, cells = q4_cell_sd(full, horizon_var), horizon = horizon, antenna = antenna,
       horizon_pooled = horizon_pooled, antenna_pooled = antenna_pooled,
       mode = "planned_simple_effects", diagnostics = diagnostics)
}

# Run the whole model fitting and analysis (many models so can take some time to run!!!)
q4 <- list()
for (i in seq_len(nrow(runs))) {
  pl <- runs$plane[i]; sg <- runs$segment[i]
  q4[[run_key(pl, sg)]] <- run_q4(planes[[pl]]$data, segments[[sg]], planes[[pl]]$horizon)
}
q4_cells <- map(names(planes), function(pl) {
  bind_rows(lapply(names(segments), function(sg) {
    mutate(q4[[run_key(pl, sg)]]$cells, segment = sg, rotation_plane = pl)
  }))
}) %>% set_names(names(planes))

# Figures

fig_cols <- c(roll_head = "#8CC8E8", roll_thorax = "#0066A6",
              pitch_head = "#F3A074", pitch_thorax = "#C94B00")
fig_col_labels <- c(roll_head = "Roll head", roll_thorax = "Roll thorax",
                    pitch_head = "Pitch head", pitch_thorax = "Pitch thorax")

theme_fig <- function() {
  theme_classic(base_size = 11) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 11),
          panel.spacing.x = unit(0.9, "lines"),
          legend.position = "top", legend.justification = "left")
}

# Labels (roll first, then pitch, as in Q1/Q2)
horizon_names <- c(flat0 = "Flat", tilted10 = "Tilted", flat0_pitch = "Flat",
                   tilted10_pitchL = "Low tilt", tilted10_pitchR = "High tilt")
horizon_levels <- c("Flat", "Tilted", "Low tilt", "High tilt")
antenna_names <- c(both = "Both", one = "One", none = "None")
plane_f <- function(x) factor(x, levels = c("roll", "pitch"), labels = c("Roll", "Pitch"))
segment_f <- function(x) factor(x, levels = c("head", "thorax"), labels = c("Head", "Thorax"))
horizon_of <- function(d) {
  key <- if ("condition" %in% names(d)) coalesce(as.character(d$condition), as.character(d$horizon_pitch))
  else as.character(d$horizon_pitch)
  unname(horizon_names[key])
}

# 1. Estimated cell values, one panel per segment x plane
# Lines connect horizon levels within an antennal state. Used for Q3 (mean) and Q4 (SD).
plot_cells <- function(d, y, title, subtitle, ylab, ymin = NULL, ymax = NULL, raw = NULL, zero_line = FALSE) {
  dodge <- position_dodge(width = 0.5)
  p <- ggplot(d, aes(horizon, .data[[y]], group = antennae, colour = plot_group,
                     shape = antennae, linetype = antennae))
  if (zero_line) p <- p + geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey75")
  if (!is.null(raw)) {
    p <- p + geom_point(data = raw, aes(horizon, value, group = antennae, colour = plot_group, shape = antennae),
                        inherit.aes = FALSE, position = dodge, size = 1.1, alpha = 0.22, show.legend = FALSE)
  }
  p <- p + geom_line(position = dodge, linewidth = 0.7)
  p <- if (!is.null(ymin)) {
    p + geom_pointrange(aes(ymin = .data[[ymin]], ymax = .data[[ymax]]), position = dodge,
                        linewidth = 0.5, linetype = "solid")
  } else {
    p + geom_point(position = dodge, size = 2.4)
  }
  p +
    facet_grid(segment ~ plane, scales = "free", space = "free_x") +
    scale_colour_manual(values = fig_cols, guide = "none") +
    scale_shape_manual(values = c(Both = 16, One = 17, None = 15), name = "Antennae") +
    scale_linetype_manual(values = c(Both = "solid", One = "dashed", None = "dotted"), name = "Antennae") +
    guides(shape = guide_legend(override.aes = list(colour = "grey30"))) +
    labs(title = title, subtitle = subtitle, x = "Horizon", y = ylab) +
    theme_fig() +
    theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.35))
}

# Q3: model-estimated means (+ beetle-level means as faint points)
q3_cells <- bind_rows(q3_emm) %>%
  mutate(plot_group = paste(rotation_plane, segment, sep = "_"),
         plane = plane_f(rotation_plane), segment = segment_f(segment),
         horizon = factor(horizon_of(.), levels = horizon_levels),
         antennae = factor(antenna_names[as.character(antennae_grp)], levels = antenna_names))

q3_raw <- analysis_data %>%
  filter(!is.na(rotation_plane), !is.na(antennae_grp)) %>%
  mutate(horizon_key = if_else(rotation_plane == "roll", as.character(condition), as.character(horizon_pitch))) %>%
  filter(!is.na(horizon_key)) %>%
  pivot_longer(c(head_angle_clean, thorax_angle_clean), names_to = "segment", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(segment = str_remove(segment, "_angle_clean")) %>%
  group_by(rotation_plane, segment, horizon_key, antennae_grp, beetle_id) %>%
  summarise(value = mean(value), .groups = "drop") %>%
  mutate(plot_group = paste(rotation_plane, segment, sep = "_"),
         plane = plane_f(rotation_plane), segment = segment_f(segment),
         horizon = factor(unname(horizon_names[horizon_key]), levels = horizon_levels),
         antennae = factor(antenna_names[as.character(antennae_grp)], levels = antenna_names))

p_q3_cells <- plot_cells(q3_cells, "emmean", ymin = "lower.CL", ymax = "upper.CL", raw = q3_raw, zero_line = TRUE,
                         title = "Q3: Sensory effects on mean orientation",
                         subtitle = "Large symbols = model estimates with 95% CI; faint points = beetle means",
                         ylab = "Mean angle (°)")
save_svg(p_q3_cells, "q3-mean-orientation-by-condition", width = 8, height = 6.5)

# Q4: model-estimated residual SD
# Per-beetle within-cell SD (raw counterpart of the DHGLM residual SD); needs >= 3 frames per beetle x cell
q4_raw <- analysis_data %>%
  filter(!is.na(rotation_plane), !is.na(antennae_grp)) %>%
  mutate(horizon_key = if_else(rotation_plane == "roll", as.character(condition), as.character(horizon_pitch))) %>%
  filter(!is.na(horizon_key)) %>%
  pivot_longer(c(head_angle_clean, thorax_angle_clean), names_to = "segment", values_to = "angle") %>%
  filter(!is.na(angle)) %>%
  mutate(segment = str_remove(segment, "_angle_clean")) %>%
  group_by(rotation_plane, segment, horizon_key, antennae_grp, beetle_id) %>%
  summarise(n_frames = n(), value = sd(angle), .groups = "drop") %>%
  filter(n_frames >= 3, is.finite(value)) %>%
  mutate(plot_group = paste(rotation_plane, segment, sep = "_"),
         plane = plane_f(rotation_plane), segment = segment_f(segment),
         horizon = factor(unname(horizon_names[horizon_key]), levels = horizon_levels),
         antennae = factor(antenna_names[as.character(antennae_grp)], levels = antenna_names))

q4_cells_plot <- bind_rows(q4_cells) %>%
  mutate(plot_group = paste(rotation_plane, segment, sep = "_"),
         plane = plane_f(rotation_plane), segment = segment_f(segment),
         horizon = factor(horizon_of(.), levels = horizon_levels),
         antennae = factor(antenna_names[as.character(antennae_grp)], levels = antenna_names))

p_q4_cells <- plot_cells(q4_cells_plot, "residual_sd", ymin = "sd_low", ymax = "sd_high", raw = q4_raw,
                         title = "Q4: Sensory effects on angular spread",
                         subtitle = "Large symbols = DHGLM residual SD with Wald 95% CI; faint points = beetle SDs",
                         ylab = "Residual SD (°)") +
  scale_y_continuous(limits = c(0, NA))
save_svg(p_q4_cells, "q4-residual-sd-by-condition", width = 8, height = 6.5)

# 2. Planned contrasts as forest plots
# Antenna contrasts are shown as manipulated minus intact (sign / ratio flipped vs the tables).

context_ant <- c(both = "both antennae", one = "one antenna", none = "no antennae")

forest_data <- function(res_list, kind) {
  imap_dfr(res_list, function(res, key) {
    pl <- word(key, 1); sg <- word(key, 2)
    if (kind == "diff") {   # Q3: differences in mean angle
      h <- res$horizon %>%
        transmute(family = "Horizon effect",
                  pair = recode(contrast, "tilted - flat" = "Tilted vs flat", "low - flat" = "Low vs flat",
                                "high - flat" = "High vs flat", "low - high" = "Low vs high"),
                  context = context_ant[as.character(antennae_grp)],
                  est = estimate, lo = lower.CL, hi = upper.CL, p_holm)
      within <- as.character(res$antenna[[intersect(c("condition", "horizon_pitch"), names(res$antenna))[1]]])
      a <- res$antenna %>%
        transmute(family = "Antenna effect",
                  pair = recode(contrast, "both - one" = "One vs both", "both - none" = "None vs both",
                                "one - none" = "None vs one"),
                  context = tolower(unname(horizon_names[within])),
                  est = -estimate, lo = -upper.CL, hi = -lower.CL, p_holm)
    } else {                # Q4: SD ratios
      h <- res$horizon %>%
        transmute(family = "Horizon effect",
                  pair = recode(str_trim(str_remove(contrast, "\\|.*")),
                                "tilted10 - flat0" = "Tilted vs flat",
                                "tilted10_pitchL - flat0_pitch" = "Low vs flat",
                                "tilted10_pitchR - flat0_pitch" = "High vs flat",
                                "tilted10_pitchL - tilted10_pitchR" = "Low vs high"),
                  context = context_ant[str_match(contrast, "antennae = (\\w+)")[, 2]],
                  est = sd_ratio, lo = sd_ratio_low, hi = sd_ratio_high, p_holm)
      a <- res$antenna %>%
        transmute(family = "Antenna effect",
                  pair = recode(str_trim(str_remove(contrast, "\\|.*")), "both - one" = "One vs both",
                                "both - none" = "None vs both", "one - none" = "None vs one"),
                  context = tolower(unname(horizon_names[str_match(contrast, "horizon = (\\w+)")[, 2]])),
                  est = 1 / sd_ratio, lo = 1 / sd_ratio_high, hi = 1 / sd_ratio_low, p_holm)
    }
    bind_rows(h, a) %>%
      mutate(rotation_plane = pl, segment = sg, label = paste(pair, context, sep = " | "))
  })
}

plot_forest <- function(d, xlab, title, subtitle, ref, logx = FALSE) {
  block_levels <- c("Roll | Horizon effect", "Roll | Antenna effect", "Pitch | Horizon effect", "Pitch | Antenna effect")
  block_labels <- c("Roll:\nhorizon\neffect", "Roll:\nantenna\neffect", "Pitch:\nhorizon\neffect", "Pitch:\nantenna\neffect")
  d <- d %>%
    mutate(plane = plane_f(rotation_plane),
           plot_group = paste(rotation_plane, segment, sep = "_"),
           segment = segment_f(segment),
           block = factor(paste(plane, family, sep = " | "), levels = block_levels, labels = block_labels),
           sig = if_else(p_holm < 0.05, plot_group, "ns"))
  # one y scale per block, first contrast at the top
  d$row <- factor(interaction(d$block, d$label, drop = TRUE, lex.order = TRUE))
  d$row <- factor(d$row, levels = rev(unique(as.character(d$row[order(d$block)]))))
  row_labels <- setNames(as.character(d$label), as.character(d$row))
  
  ggplot(d, aes(x = est, y = row, colour = plot_group, fill = sig, group = segment)) +
    geom_vline(xintercept = ref, linewidth = 0.45, colour = "grey60") +
    geom_linerange(aes(xmin = lo, xmax = hi), orientation = "y", linewidth = 0.55,
                   position = position_dodge(width = 0.65)) +
    geom_point(shape = 21, size = 1.9, stroke = 0.8, position = position_dodge(width = 0.65)) +
    facet_grid(block ~ ., scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = row_labels) +
    {if (logx) scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2)) else scale_x_continuous()} +
    scale_colour_manual(values = fig_cols, breaks = names(fig_cols), labels = fig_col_labels, name = NULL) +
    scale_fill_manual(values = c(fig_cols, ns = "white"), guide = "none") +
    guides(colour = guide_legend(nrow = 1, override.aes = list(fill = fig_cols, shape = 21, linetype = 0, linewidth = 0, size = 2))) +
    labs(title = title, subtitle = subtitle, x = xlab, y = NULL) +
    theme_fig() +
    theme(strip.text.y = element_text(angle = 0, hjust = 0, face = "bold", size = 10),
          panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.35),
          panel.spacing.y = unit(0.8, "lines"))
}

q3_forest <- forest_data(q3, "diff")
p_q3_forest <- plot_forest(q3_forest, ref = 0,
                           xlab = "Difference in mean angle (°)",
                           title = "Q3: Planned contrasts for mean orientation",
                           subtitle = "Estimate with 95% CI; filled = Holm p < 0.05\nAntenna effects = manipulated minus intact (sign flipped vs tables)")
save_svg(p_q3_forest, "q3-planned-contrasts", width = 8, height = 9)

q4_forest <- forest_data(q4, "ratio")
p_q4_forest <- plot_forest(q4_forest, ref = 1, logx = TRUE,
                           xlab = "SD ratio (first / second condition; < 1 = less variable)",
                           title = "Q4: Planned contrasts for angular spread",
                           subtitle = "Ratio of residual SDs with 95% CI; filled = Holm p < 0.05\nAntenna effects = manipulated / intact (ratio inverted vs tables)")
save_svg(p_q4_forest, "q4-sd-ratio-contrasts", width = 8, height = 9)


# 6. Q5 - sensory effects on the head-thorax slope --------------------------
# Full factorial: head ~ thorax_dev * Horizon * Antennae, with the random-effect structure

# Prepare the data (keep frames with both head and thorax)
q5_data <- function(d) {
  d %>%
    filter(!is.na(head_angle_clean), !is.na(thorax_angle_clean)) %>%
    group_by(beetle_id) %>%
    mutate(thorax_dev = thorax_angle_clean - mean(thorax_angle_clean)) %>%
    ungroup()
}

# Model with three way interactions as we want to ask if coupling change with both factors
# Random-effects structure is so that each beetle has its own intercept and own thorax to head coupling
# I checked a diagonal random effect model for correlation between slope and intercept
fit_q5 <- function(data, horizon_var) {
  fixed <- paste0("head_angle_clean ~ thorax_dev * ", horizon_var, " * antennae_grp")
  fit <- glmmTMB(
    as.formula(paste(fixed, "+ (1 + thorax_dev | beetle_id)")),
    dispformula = as.formula(paste0("~ ", horizon_var, " * antennae_grp")),
    family = gaussian, data = data, REML = TRUE
  )
  if (fit$fit$convergence != 0 || !isTRUE(fit$sdr$pdHess))
    stop("Q5 model failed to converge for ", horizon_var)
  fit
}

# Fits the two models for each axes
q5_models <- list(
  roll = fit_q5(q5_data(roll_data), "condition"),
  pitch = fit_q5(q5_data(pitch_data), "horizon_pitch")
)

# Type III test of the effects
q5_anova <- map(q5_models, ~ car::Anova(.x, type = 3, component = "cond"))

# Slope for coupling strength for each conditions
cell_slopes <- function(model, horizon_var) {
  emtrends(
    model,
    specs = as.formula(paste0("~ ", horizon_var, " * antennae_grp")),
    var = "thorax_dev",
    component = "cond"
  ) %>%
    summary(infer = c(TRUE, TRUE)) %>%
    as.data.frame() %>%
    add_ci()
}

# Slope comparisons
slope_contrasts <- function(model, horizon_var) {
  trends <- function(specs)
    emtrends(
      model,
      specs = as.formula(specs),
      var = "thorax_dev",
      component = "cond"
    )
  
  summarise_family <- function(x)
    as.data.frame(summary(x, infer = c(TRUE, TRUE), adjust = "none")) %>%
    mutate(p_holm = p.adjust(p.value, method = "holm"))
  
  list(
    horizon = summarise_family(
      horizon_contrasts(
        trends(paste0("~ ", horizon_var, " | antennae_grp")),
        horizon_var
      )
    ),
    antenna = summarise_family(
      pairs(
        trends(paste0("~ antennae_grp | ", horizon_var)),
        adjust = "none"
      )
    ),
    # Pooled contrasts, averaged over the other factor
    horizon_pooled = summarise_family(horizon_contrasts(trends(paste0("~ ", horizon_var)), horizon_var)),
    antenna_pooled = summarise_family(pairs(trends("~ antennae_grp"), adjust = "none"))
  )
}

# Results
horizon_axis <- c(roll = "condition", pitch = "horizon_pitch")
q5_slopes <- imap(q5_models, ~ cell_slopes(.x, horizon_axis[[.y]]))
q5_conditional_sd <- imap(q5_models, ~ scale_emmeans(.x, horizon_axis[[.y]]))
q5_contrasts <- imap(q5_models, ~ slope_contrasts(.x, horizon_axis[[.y]]))

# Figures
# Q5 main figure: cell-specific head-thorax slopes, roll and pitch together
q5_cells <- imap_dfr(q5_slopes, ~ mutate(.x, rotation_plane = .y)) %>%
  mutate(plot_group = paste0(rotation_plane, "_thorax"),   # dark shade; use "_head" for the light one
         plane = plane_f(rotation_plane),
         segment = factor("Head ~ thorax"),
         horizon = factor(horizon_of(.), levels = horizon_levels),
         antennae = factor(antenna_names[as.character(antennae_grp)], levels = antenna_names))

p_q5_cells <- plot_cells(q5_cells, "thorax_dev.trend", ymin = "lower.CL", ymax = "upper.CL",
                         title = "Q5: Sensory effects on head-thorax coupling",
                         subtitle = "Model-estimated slope with 95% CI; dashed line = 1 (head moves one-to-one with thorax)",
                         ylab = "Head-thorax slope") +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5, colour = "grey45") +
  facet_grid(segment ~ plane, scales = "free_x", space = "free_x")   # shared y axis so roll and pitch are comparable

save_svg(p_q5_cells, "q5-head-thorax-slopes", width = 8, height = 4)

# Regression line per condition
q5_line_data <- function(pl) {
  hv <- horizon_axis[[pl]]
  d <- q5_data(if (pl == "roll") roll_data else pitch_data)
  xs <- seq(min(d$thorax_dev), max(d$thorax_dev), length.out = 60)
  specs <- c(hv, "antennae_grp", "thorax_dev")
  tidy <- function(at) emmeans(q5_models[[pl]], specs, at = at, component = "cond") %>%
    as.data.frame() %>% add_ci() %>% rename(horizon_key = all_of(hv))
  
  rng <- d %>% rename(horizon_key = all_of(hv)) %>%
    group_by(horizon_key, antennae_grp) %>%
    summarise(x_lo = min(thorax_dev), x_hi = max(thorax_dev), .groups = "drop")
  icpt <- tidy(list(thorax_dev = 0)) %>% select(horizon_key, antennae_grp, intercept = emmean)
  
  lines <- tidy(list(thorax_dev = xs)) %>%
    left_join(rng, by = c("horizon_key", "antennae_grp")) %>%
    filter(thorax_dev >= x_lo, thorax_dev <= x_hi) %>%
    left_join(icpt, by = c("horizon_key", "antennae_grp")) %>%
    mutate(ref_11 = intercept + thorax_dev)
  pts <- d %>% rename(horizon_key = all_of(hv))
  
  prep <- function(x) x %>% mutate(
    horizon = factor(unname(horizon_names[as.character(horizon_key)]), levels = horizon_levels),
    antennae = factor(antenna_names[as.character(antennae_grp)], levels = antenna_names))
  list(lines = prep(lines), points = prep(pts))
}

plot_q5_lines <- function(pl) {
  dd <- q5_line_data(pl)
  col <- fig_cols[[paste0(pl, "_thorax")]]
  ggplot(dd$points, aes(thorax_dev, head_angle_clean)) +
    geom_point(colour = col, alpha = 0.25, size = 1.1) +
    geom_line(data = dd$lines, aes(y = ref_11), linetype = "dashed", linewidth = 0.5, colour = "grey45") +
    geom_ribbon(data = dd$lines, aes(y = emmean, ymin = lower.CL, ymax = upper.CL),
                fill = col, alpha = 0.18, colour = NA) +
    geom_line(data = dd$lines, aes(y = emmean), colour = col, linewidth = 1) +
    facet_grid(antennae ~ horizon) +
    coord_equal() +
    labs(title = paste0("Q5: Head-thorax regression per condition - ", plane_f(pl)),
         subtitle = "Solid = model fit with 95% CI; dashed = 1:1 slope through the same intercept; points = frames",
         x = "Thorax deviation from beetle mean (°)", y = "Head angle (°)") +
    theme_fig()
}

save_svg(plot_q5_lines("roll"),  "q5-roll-regression-by-condition",  width = 7, height = 8)
save_svg(plot_q5_lines("pitch"), "q5-pitch-regression-by-condition", width = 9, height = 8)

# 7. Model diagnostics (Normality) ------------------------------------------------------

# Q1: residuals against fitted values and residual QQ plot
par(mfrow = c(2, 2))
for (pl in c("roll", "pitch")) {
  plot(fitted(q1_models[[pl]]), resid(q1_models[[pl]]), main = paste("Q1", pl),
       xlab = "Fitted", ylab = "Residuals")
  abline(h = 0, lty = 2)
  qqnorm(resid(q1_models[[pl]]), main = paste("Q1", pl))
  qqline(resid(q1_models[[pl]]))
}
par(mfrow = c(1, 1))

# Q2: QQ plot of the log SD ratios per beetle
par(mfrow = c(1, 2))
for (pl in c("roll", "pitch")) {
  lr <- filter(q2_sd_by_beetle, rotation_plane == pl)$log_sd_ratio
  qqnorm(lr, main = paste("Q2", pl))
  qqline(lr)
}
par(mfrow = c(1, 1))

# Q3 and Q5: simulated residuals
diag_models <- c(setNames(q3_models, paste("Q3", names(q3_models))),
                   setNames(q5_models, paste("Q5", names(q5_models))))
  
# Residual QQ plot and residuals against predicted values
for (nm in names(diag_models)) {
    message("Diagnostics: ", nm)
    plot(DHARMa::simulateResiduals(diag_models[[nm]], n = 1000, seed = 42))
    mtext(nm, side = 3, outer = TRUE, line = -1.5, font = 2)
}
  
# Beetle random intercepts
par(mfrow = c(2, 3))
for (nm in names(diag_models)) {
    re <- glmmTMB::ranef(diag_models[[nm]])$cond$beetle_id[, 1]
    qqnorm(re, main = nm)
    qqline(re)
}
par(mfrow = c(1, 1))

# 8. Review file -----------------------------------------------------------

# Summaries all data obtained (Created using generative AI: Claude Sonnet 5.5)
q4_fit_valid <- function(result) {
  health <- filter(result$diagnostics, model == "full scale")
  isTRUE(result$full$ok) && nrow(health) == 1 &&
    isTRUE(health$fit_returned[1]) && isTRUE(health$finite_ml[1]) &&
    isTRUE(health$finite_coefficients[1]) && isTRUE(health$finite_scale_random_effect[1])
}

q4_structure <- "random beetle intercept in location and scale; DHGLM fit-validity checks applied"

q5_ok <- map_lgl(
  q5_models,
  ~ .x$fit$convergence == 0 && isTRUE(.x$sdr$pdHess)
)

q5_structure <- "correlated random intercept + random thorax_dev slope"

model_checks <- bind_rows(
  mutate(q3_checks, random_structure = "location random intercept; fixed-effect dispersion model"),
  tibble(
    model = c("Q4 DHGLM roll head", "Q4 DHGLM roll thorax", "Q4 DHGLM pitch head", "Q4 DHGLM pitch thorax",
              "Q5 roll", "Q5 pitch"),
    converged = c(
      map_lgl(q4, q4_fit_valid),
      q5_ok
    ),
    random_structure = c(
      rep(q4_structure, 4),
      rep(q5_structure, 2)
    )
  )
)

banner <- function(title, lead = "\n") {
  cat(lead, strrep("=", 60), "\n", title, "\n", strrep("=", 60), "\n", sep = "")
}

print_q3_block <- function(label, result) {
  cat("\n", label, "\n", sep = "")
  print.data.frame(as.data.frame(result$anova), row.names = TRUE)
  cat("Post-hoc mode:", result$mode, "\n")
  cat("Horizon family:\n")
  print_contrasts(result$horizon)
  cat("Antenna family:\n")
  print_contrasts(result$antenna)
  cat("Pooled horizon contrasts (averaged over antennae):\n")
  print_contrasts(result$horizon_pooled)
  cat("Pooled antenna contrasts (averaged over horizon):\n")
  print_contrasts(result$antenna_pooled)
}

print_q4_block <- function(label, x) {
  cat("\n", label, "\n", sep = "")
  cat("Post-hoc mode:", x$mode, "\n")
  cat("Cell-specific population residual SD estimates:\n")
  if (nrow(x$cells) > 0) print_rounded(x$cells)
  cat("Horizon family:\n")
  print_contrasts(x$horizon)
  cat("Antenna family:\n")
  print_contrasts(x$antenna)
  cat("Pooled horizon contrasts (averaged over antennae):\n")
  print_contrasts(x$horizon_pooled)
  cat("Pooled antenna contrasts (averaged over horizon):\n")
  print_contrasts(x$antenna_pooled)
  
  if (x$full$ok && !is.null(x$full$model$tau_coeff)) {
    cat("Scale random-effect variance model (estimate and SE only; no zero-variance test):\n")
    tau <- as.data.frame(x$full$model$tau_coeff)
    print.data.frame(tau[, intersect(c("Estimate", "Std. Error"), names(tau)), drop = FALSE], row.names = TRUE)
  }
  
  cat("DHGLM fit / likelihood diagnostics:\n")
  print.data.frame(mutate(x$diagnostics, across(c(ml, rl, caic), ~ round(.x, 3))), row.names = FALSE)
}

print_q5_block <- function(pl, lead = "\n") {
  PL <- toupper(pl)
  cat(lead, PL, " - TYPE III TESTS\n", sep = "")
  print.data.frame(as.data.frame(q5_anova[[pl]]), row.names = TRUE)
  cat("\n", PL, " - CELL-SPECIFIC HEAD-THORAX SLOPES\n", sep = "")
  print.data.frame(q5_slopes[[pl]], row.names = FALSE)
  cat("\n", PL, " - CONDITIONAL HEAD RESIDUAL SD (DESCRIPTIVE)\n", sep = "")
  sd_table <- q5_conditional_sd[[pl]] %>%
    select(any_of(c("condition", "horizon_pitch", "antennae_grp", "residual_sd", "sd_low", "sd_high")))
  print_rounded(sd_table)
  cat("\n", PL, " - HORIZON SLOPE FAMILY, HOLM\n", sep = "")
  print.data.frame(q5_contrasts[[pl]]$horizon, row.names = FALSE)
  cat("\n", PL, " - ANTENNA SLOPE FAMILY, HOLM\n", sep = "")
  print.data.frame(q5_contrasts[[pl]]$antenna, row.names = FALSE)
  cat("\n", PL, " - POOLED HORIZON SLOPE CONTRASTS (AVERAGED OVER ANTENNAE), HOLM\n", sep = "")
  print.data.frame(q5_contrasts[[pl]]$horizon_pooled, row.names = FALSE)
  cat("\n", PL, " - POOLED ANTENNA SLOPE CONTRASTS (AVERAGED OVER HORIZON), HOLM\n", sep = "")
  print.data.frame(q5_contrasts[[pl]]$antenna_pooled, row.names = FALSE)
}

report <- capture.output({
  cat("SCRIPT VERSION:", SCRIPT_VERSION, "\n")
  cat("glmmTMB version:", as.character(packageVersion("glmmTMB")), "\n")
  cat("dhglm version:", as.character(packageVersion("dhglm")), "\n\n")
  
  banner("DATA AND PRIMARY ANALYSIS DEFINITIONS", lead = "")
  cat("Rows retained:", nrow(analysis_data), "\n")
  cat("Beetles:", n_distinct(analysis_data$beetle_id), "\n")
  cat("Head outliers removed:", sum(analysis_data$outlier_remove_head, na.rm = TRUE), "\n")
  cat("Thorax outliers removed:", sum(analysis_data$outlier_remove_thorax, na.rm = TRUE), "\n")
  cat("Antennae factor: both / one / none (left and right pooled).\n")
  cat("Primary models contain no external covariates.\n")
  cat("Roll and pitch are analysed separately.\n")
  cat("Pitch horizon levels: flat / low tilt / high tilt.\n")
  cat("Pitch mapping: sideL = low tilt; sideR = high tilt.\n\n")
  
  cat("Q3 model formula\n")
  cat("Y ~ Horizon * antennae_grp + (1 | beetle_id)\n")
  cat("Residual dispersion is allowed to differ across Horizon x Antennae cells as a nuisance variance model.\n\n")
  
  cat("Q4 DHGLM formulas\n")
  cat("Location:   Y ~ Horizon * antennae_grp + (1 | beetle_id)\n")
  cat("Scale: log(phi) ~ Horizon * antennae_grp + (1 | beetle_id)\n")
  cat("Pitch uses horizon_pitch; roll uses condition.\n\n")
  
  cat("Q5 full fixed-effect structure:\n")
  cat("head_angle_clean ~ thorax_dev * Horizon * antennae_grp\n")
  cat("The terms containing thorax_dev determine the head-thorax slope.\n\n")
  
  banner("Q1. BASELINE HEAD-THORAX SLOPE", lead = "")
  cat("Baseline = flat horizon + both antennae.\n")
  cat("Model: head_angle ~ within-beetle thorax deviation + (1 | beetle)\n")
  cat("Slope 1 = one-to-one movement; slope 0 = no linear coupling.\n")
  cat("Scientific reference = 1. Primary interval = beetle-cluster bootstrap.\n\n")
  print_rounded(select(q1_summary, -t_vs_1, -p_vs_1))
  
  banner("Q2. BASELINE HEAD VERSUS THORAX SD")
  cat("Each beetle contributes paired head and thorax SDs from matched frames.\n")
  cat("Test: one-sample t-test of log(head SD / thorax SD) against 0.\n")
  cat("SD ratio < 1 means the head varies less than the thorax.\n\n")
  print_rounded(mutate(q2_summary, p_value = format_p(p_value)))
  
  banner("Q3. DO SENSORY CONDITIONS CHANGE AVERAGE ORIENTATION?")
  cat("Head and thorax are tested separately.\n")
  cat("Planned simple contrasts are reported regardless of the interaction p-value.\n")
  cat("The interaction test describes non-additivity rather than gating which contrasts are shown.\n")
  cat("Horizon and antenna contrast sets are separate Holm families.\n")
  cat("Pooled contrasts (main effects averaged over the other factor) are also given, as separate Holm families.\n")
  print_q3_block("ROLL - HEAD MEAN", q3[["roll head"]])
  print_q3_block("ROLL - THORAX MEAN", q3[["roll thorax"]])
  print_q3_block("PITCH - HEAD MEAN", q3[["pitch head"]])
  print_q3_block("PITCH - THORAX MEAN", q3[["pitch thorax"]])
  cat("\nQ3 MODEL-ESTIMATED MEANS\n")
  cat("Roll:\n")
  print_rounded(q3_emm$roll)
  cat("\nPitch:\n")
  print_rounded(q3_emm$pitch)
  
  banner("Q4. DO SENSORY CONDITIONS CHANGE ANGULAR SPREAD?")
  cat("DHGLM: beetle random intercepts are fitted in both location and scale.\n")
  cat("Location: Y ~ Horizon * Antennae + (1 | beetle)\n")
  cat("Full scale: log(phi) ~ Horizon * Antennae + (1 | beetle)\n")
  cat("Planned simple-effect Wald contrasts are always estimated from the full cell model.\n")
  cat("Horizon and antenna contrast sets are Holm-corrected as separate families.\n")
  cat("Pooled Wald contrasts (averaged over the other factor) come from the same full cell model, as separate Holm families.\n")
  cat("SD ratios are exp(log-variance contrast / 2).\n")
  print_q4_block("ROLL - HEAD SPREAD", q4[["roll head"]])
  print_q4_block("ROLL - THORAX SPREAD", q4[["roll thorax"]])
  print_q4_block("PITCH - HEAD SPREAD", q4[["pitch head"]])
  print_q4_block("PITCH - THORAX SPREAD", q4[["pitch thorax"]])
  
  banner("Q5. DO SENSORY CONDITIONS CHANGE THE HEAD-THORAX SLOPE?")
  cat("Full crossed model: head ~ thorax_dev * Horizon * Antennae.\n")
  cat("The three-way thorax_dev:Horizon:Antennae term tests whether the effect\n")
  cat("of one sensory modality on coupling depends on the other modality.\n")
  cat("Cell-specific slopes are reported for every crossed sensory combination.\n")
  cat("Slope post-hocs are grouped into two Holm families:\n")
  cat("  1) horizon slope contrasts within antenna conditions;\n")
  cat("  2) antenna slope contrasts within horizon conditions.\n")
  cat("Pooled slope contrasts (averaged over the other factor) are also given, as separate Holm families.\n\n")
  cat("Q5 random-effect structure: correlated random intercept + random thorax_dev slope.\n\n")
  cat("Q5 conditional head residual SD is also reported descriptively from the existing dispersion model.\n")
  cat("This is residual head variation after accounting for thorax coupling, not a replacement for Q4.\n\n")
  print_q5_block("roll", lead = "")
  print_q5_block("pitch")
  
  banner("MODEL CHECKS")
  print.data.frame(as.data.frame(model_checks), row.names = FALSE)
  cat("\nTECHNICAL NOTE ON Q4\n")
  cat("Q4 is fitted with dhglm so beetle identity can enter both the mean and dispersion models.\n")
  cat("Q4 planned simple contrasts use DHGLM Wald tests from the full cell model.\n")
  cat("For Q4, the model-check table reports TRUE only when the DHGLM fit returned and its likelihood, coefficients, and scale random-effect estimate are all finite.\n")
  
  cat("\nRandom-effect structures:\n")
  cat("Roll head Q3:\n"); print(VarCorr(q3_models[["roll head"]]))
  cat("Roll thorax Q3:\n"); print(VarCorr(q3_models[["roll thorax"]]))
  cat("Pitch head Q3:\n"); print(VarCorr(q3_models[["pitch head"]]))
  cat("Pitch thorax Q3:\n"); print(VarCorr(q3_models[["pitch thorax"]]))
  cat("Q5 roll:\n"); print(VarCorr(q5_models$roll))
  cat("Q5 pitch:\n"); print(VarCorr(q5_models$pitch))
})

writeLines(report, output_file)

cat("\nOutput files:\n")
cat("-", output_file, "\n")
cat("- raw_figures/\n")
cat("\nAnalysis complete. Review", output_file, "first.\n")