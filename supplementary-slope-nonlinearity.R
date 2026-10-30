# Head stabilisation in Kheper lamarcki - supplementary: non-linear head-thorax slope (Q1)
# Companion to main-analysis.R: tests if a curve explains the baseline head-thorax relationship better than the straight line
# Baseline = flat horizon + both antennae, same data and random effect as Q1 in the main script
# Models (ML fits, likelihood ratio test, AIC and BIC against the linear model):
#   linear (main), odd cubic (flat centre, steeper edges), quadratic, cubic, natural spline
# A GAM smooth is added as a data-driven check (edf close to 1 = straight line)

library(tidyverse)
library(lme4)
library(lmerTest)
library(mgcv)

set.seed(42)

SCRIPT_VERSION <- "v1-q1-nonlinearity"
data_file <- "data_stabilisation_processed.csv"
output_file <- "reviews/q1_nonlinearity_review_v1.txt"

dir.create("reviews", showWarnings = FALSE)

# 1. Data (same as main-analysis.R) ----------------------------------------

# Remove a value if both outlier rules flag it, or one rule flags it and the row has a comment
flag_outlier <- function(z3, iqr, has_comment) {
  z3 <- coalesce(z3, FALSE)
  iqr <- coalesce(iqr, FALSE)
  (z3 & iqr) | ((z3 | iqr) & has_comment)
}

analysis_data <- read.csv(data_file, fileEncoding = "UTF-8") %>%
  mutate(
    rotation_plane = factor(if_else(view == "front", "roll", "pitch"), levels = c("roll", "pitch")),
    has_comment = !is.na(comments) & str_trim(comments) != "",
    head_angle_clean = if_else(flag_outlier(outlier_z3_head, outlier_iqr_head, has_comment),
                               NA_real_, head_angle_corrected),
    thorax_angle_clean = if_else(flag_outlier(outlier_z3_thorax, outlier_iqr_thorax, has_comment),
                                 NA_real_, thorax_angle_corrected)
  ) %>%
  filter(is.na(comments) | !str_detect(comments, regex("discard", ignore_case = TRUE)))

# Baseline = flat horizon + both antennae, thorax is centred within beetle (x in units of 10 degrees)
q1_data <- analysis_data %>%
  filter(condition == "flat0", antennae == "both",
         !is.na(head_angle_clean), !is.na(thorax_angle_clean)) %>%
  group_by(rotation_plane, beetle_id) %>%
  mutate(x = (thorax_angle_clean - mean(thorax_angle_clean)) / 10) %>%
  ungroup()

# 2. Models ----------------------------------------------------------------

# Fixed-effect shapes, the random beetle intercept is the same as in the main script
shapes <- list(
  linear = "x",
  odd_cubic = "x + I(x^3)",
  quadratic = "x + I(x^2)",
  cubic = "x + I(x^2) + I(x^3)",
  spline = "splines::ns(x, df = 4)"
)

# Fitted with ML so the fixed effects can be compared
fit_shapes <- function(d) {
  map(shapes, ~ lmer(as.formula(paste0("head_angle_clean ~ ", .x, " + (1 | beetle_id)")),
                     data = d, REML = FALSE))
}

fits <- map(c(roll = "roll", pitch = "pitch"), ~ fit_shapes(filter(q1_data, rotation_plane == .x)))

# Likelihood ratio test, AIC and BIC against the linear model
comparison <- imap_dfr(fits, function(f, pl) {
  imap_dfr(f, function(m, nm) {
    chisq <- 2 * (as.numeric(logLik(m)) - as.numeric(logLik(f$linear)))
    df <- attr(logLik(m), "df") - attr(logLik(f$linear), "df")
    tibble(plane = pl, model = nm, n_par = attr(logLik(m), "df"), AIC = AIC(m), BIC = BIC(m),
           Chisq = chisq, Df = df, p_value = if (df > 0) pchisq(chisq, df, lower.tail = FALSE) else NA_real_)
  })
}) %>%
  group_by(plane) %>%
  mutate(delta_AIC = AIC - min(AIC), delta_BIC = BIC - min(BIC)) %>%
  ungroup()

# Local slope (head degrees per thorax degree) from the cubic model: centre and the 10th / 90th percentile
local_slopes <- imap_dfr(fits, function(f, pl) {
  b <- fixef(f$cubic)[c("x", "I(x^2)", "I(x^3)")]
  x <- filter(q1_data, rotation_plane == pl)$x
  at <- c("10th percentile" = unname(quantile(x, 0.10)), "centre" = 0, "90th percentile" = unname(quantile(x, 0.90)))
  tibble(plane = pl, position = names(at), thorax_dev_deg = unname(at) * 10,
         slope = (b[[1]] + 2 * b[[2]] * at + 3 * b[[3]] * at^2) / 10,
         linear_slope = fixef(f$linear)[["x"]] / 10)
})

# GAM smooth with a beetle random intercept (edf close to 1 = straight line)
gam_edf <- map_dfr(c("roll", "pitch"), function(pl) {
  d <- filter(q1_data, rotation_plane == pl) %>% mutate(beetle_id = factor(beetle_id))
  g <- gam(head_angle_clean ~ s(x, k = 5) + s(beetle_id, bs = "re"), data = d, method = "REML")
  tibble(plane = pl, edf = summary(g)$s.table["s(x)", "edf"])
})

# 3. Plot (shown on screen, not saved) ------------------------------------

# Curves with 95% pointwise confidence bands (fixed effects only, beetle intercept left out)
curves <- imap_dfr(fits, function(f, pl) {
  x <- filter(q1_data, rotation_plane == pl)$x
  grid <- tibble(x = seq(min(x), max(x), length.out = 100))
  map_dfr(c("linear", "odd_cubic", "spline"), function(nm) {
    m <- f[[nm]]
    X <- model.matrix(delete.response(terms(m, fixed.only = TRUE)), data = mutate(grid, head_angle_clean = 0))
    fit <- drop(X %*% fixef(m))
    se <- sqrt(rowSums((X %*% as.matrix(vcov(m))) * X))
    mutate(grid, rotation_plane = pl, model = nm, fit = fit, lower = fit - 1.96 * se, upper = fit + 1.96 * se)
  })
})

print(
  ggplot(q1_data, aes(x * 10, head_angle_clean)) +
    geom_point(alpha = 0.3) +
    geom_ribbon(data = curves, aes(x * 10, ymin = lower, ymax = upper, fill = model),
                inherit.aes = FALSE, alpha = 0.15) +
    geom_line(data = curves, aes(x * 10, fit, colour = model), linewidth = 1) +
    facet_wrap(~ rotation_plane, scales = "free") +
    labs(x = "Thorax deviation (deg)", y = "Head angle (deg)", colour = NULL, fill = NULL) +
    theme_classic()
)

# 4. Review file -----------------------------------------------------------

banner <- function(title, lead = "\n") {
  cat(lead, strrep("=", 60), "\n", title, "\n", strrep("=", 60), "\n", sep = "")
}
format_p <- function(p) ifelse(is.na(p), NA_character_, formatC(p, format = "e", digits = 2))
print_rounded <- function(x, digits = 3) {
  print.data.frame(as.data.frame(mutate(x, across(where(is.numeric), ~ round(.x, digits)))), row.names = FALSE)
}

report <- capture.output({
  cat("SCRIPT VERSION:", SCRIPT_VERSION, "\n\n")
  
  banner("DEFINITIONS", lead = "")
  cat("Baseline = flat horizon + both antennae, head ~ thorax deviation + (1 | beetle), as in Q1.\n")
  cat("Frames: roll", sum(q1_data$rotation_plane == "roll"), ", pitch", sum(q1_data$rotation_plane == "pitch"), "\n")
  cat("Beetles: roll", n_distinct(filter(q1_data, rotation_plane == "roll")$beetle_id),
      ", pitch", n_distinct(filter(q1_data, rotation_plane == "pitch")$beetle_id), "\n")
  cat("Odd cubic (x + x^3) = flat centre and steeper edges, the shape expected if stabilisation is best near the middle.\n")
  cat("Likelihood ratio test, AIC and BIC are against the linear model; delta = difference to the best model (0 = best).\n")
  
  banner("MODEL COMPARISON")
  print_rounded(mutate(comparison, p_value = format_p(p_value)))
  
  banner("LOCAL SLOPE FROM THE CUBIC MODEL")
  cat("Head degrees per thorax degree at the centre and at the edges; linear_slope is the straight-line slope.\n\n")
  print_rounded(local_slopes)
  
  banner("GAM SMOOTH")
  cat("Effective degrees of freedom of the thorax smooth (1 = straight line).\n\n")
  print_rounded(gam_edf)
})

writeLines(report, output_file)

cat("\nOutput files:\n")
cat("-", output_file, "\n")
cat("\nAnalysis complete. Review", output_file, "first.\n")