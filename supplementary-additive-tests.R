# Head stabilisation in Kheper lamarcki - supplementary: additive models versus interactions
# Companion to main-analysis.R: tests if the interactions kept in the main models are worth it
# Q3: location model, Horizon * Antennae vs Horizon + Antennae (ML fits, likelihood ratio test)
# Q4: DHGLM scale model, full vs additive plus the additive omnibus tests (adjusted-profile LRT)
# Q5: head ~ thorax_dev * Horizon * Antennae, step-wise reductions (ML fits, likelihood ratio test)
# Sensitivity: the planned contrasts of the main script, from the full model and from the reduced one

library(tidyverse)
library(glmmTMB)
library(emmeans)

# Check if dhglm is installed (functions directly called later)
if (!requireNamespace("dhglm", quietly = TRUE)) {
  stop("Package 'dhglm' is required. Install it with install.packages('dhglm') and rerun the script.")
}

set.seed(42)

# Setting the factor contrast
options(contrasts = c("contr.sum", "contr.poly"))

SCRIPT_VERSION <- "v2-supplementary"
data_file <- "data_stabilisation_processed.csv"
output_file <- "reviews/supplementary_review_v2.txt"

dir.create("reviews", showWarnings = FALSE)

# 1. Data (same as main-analysis.R) ----------------------------------------

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

# 2. Helpers ---------------------------------------------------------------

# Format p values
format_p <- function(p) ifelse(is.na(p), NA_character_, formatC(p, format = "e", digits = 2))
# Round values for output
print_rounded <- function(x, digits = 3) {
  print.data.frame(as.data.frame(mutate(x, across(where(is.numeric), ~ round(.x, digits)))), row.names = FALSE)
}

# Likelihood ratio test between two nested glmmTMB models (both fitted with ML), with AIC and BIC
lrt <- function(small, big, label, tested_in) {
  chisq <- 2 * (as.numeric(logLik(big)) - as.numeric(logLik(small)))
  df <- attr(logLik(big), "df") - attr(logLik(small), "df")
  tibble(
    effect = label, tested_in = tested_in,
    Chisq = chisq, Df = df, p_value = pchisq(chisq, df = df, lower.tail = FALSE),
    AIC_reduced = AIC(small), AIC_full = AIC(big), BIC_reduced = BIC(small), BIC_full = BIC(big)
  )
}

# 3. Q3 - interaction in the mean orientation model ------------------------
# Same model as the main script (cell-specific residual variance), but with ML so fixed effects can be compared

fit_location <- function(data, response, horizon_var, sign, reml = FALSE) {
  glmmTMB(
    formula = as.formula(paste0(response, " ~ ", horizon_var, sign, "antennae_grp + (1 | beetle_id)")),
    dispformula = as.formula(paste0("~ ", horizon_var, " * antennae_grp")),
    family = gaussian, data = data, REML = reml
  )
}

q3_tests <- map_dfr(seq_len(nrow(runs)), function(i) {
  pl <- runs$plane[i]; sg <- runs$segment[i]
  hv <- planes[[pl]]$horizon
  d <- planes[[pl]]$data
  d <- d[!is.na(d[[segments[[sg]]]]), ]
  full <- fit_location(d, segments[[sg]], hv, " * ")
  additive <- fit_location(d, segments[[sg]], hv, " + ")
  lrt(additive, full, paste0(hv, ":antennae_grp"), "full model") %>%
    mutate(model = paste(pl, sg), .before = 1)
})

# 4. Q4 - interaction in the angular spread model --------------------------

# Double hierarchical model (same as the main script)
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
spread_data <- function(data, response, horizon_var) {
  data %>%
    filter(!is.na(.data[[response]])) %>%
    transmute(y = .data[[response]], beetle_id = factor(beetle_id),
              antennae_grp = droplevels(antennae_grp), horizon = droplevels(.data[[horizon_var]])) %>%
    mutate(scale_cell = interaction(horizon, antennae_grp, sep = "__", drop = TRUE))
}

fit_spread <- function(data, response, horizon_var, scale = "cell") {
  dat <- spread_data(data, response, horizon_var)
  rhs <- c(cell = "scale_cell", additive = "horizon + antennae_grp",
           horizon = "horizon", antennae = "antennae_grp")[[scale]]
  model <- fit_dhglm(dat, rhs)
  list(model = model, ok = !is.null(model), data = dat)
}

# Adjusted-profile likelihood ratio test between two nested DHGLM fits
q4_lrt <- function(full, reduced, df, effect) {
  raw <- as.numeric(reduced$model$ml) - as.numeric(full$model$ml)
  tibble(
    effect = effect, full_ml = as.numeric(full$model$ml), reduced_ml = as.numeric(reduced$model$ml),
    raw_Chisq = raw, Chisq = max(0, raw), Df = df,
    p_value = pchisq(max(0, raw), df = df, lower.tail = FALSE),
    caic_full = as.numeric(full$model$caic), caic_reduced = as.numeric(reduced$model$caic),
    test_method = "adjusted-profile LRT"
  )
}

# Wald test of the horizon effect in the additive scale model (only for 2 horizon levels)
wald_horizon <- function(dat) {
  dat$wald_term <- dat$horizon
  contrasts(dat$wald_term) <- contr.treatment(nlevels(dat$wald_term))
  cf <- fit_dhglm(dat, "wald_term + antennae_grp")$phi_coeff
  cf <- cf[startsWith(rownames(cf), "wald_term"), ]
  z <- cf[["Estimate"]] / cf[["Std. Error"]]
  list(chisq = z^2, p_value = 2 * pnorm(abs(z), lower.tail = FALSE))
}

# Full scale model vs additive, then the additive model with horizon or antennae dropped
run_q4 <- function(data, response, horizon_var) {
  full <- fit_spread(data, response, horizon_var, "cell")
  additive <- fit_spread(data, response, horizon_var, "additive")
  no_horizon <- fit_spread(data, response, horizon_var, "antennae")
  no_antenna <- fit_spread(data, response, horizon_var, "horizon")
  n_h <- nlevels(full$data$horizon)
  n_a <- nlevels(full$data$antennae_grp)
  
  tests <- bind_rows(
    q4_lrt(full, additive, (n_h - 1) * (n_a - 1), paste0(horizon_var, ":antennae_grp")),
    q4_lrt(additive, no_horizon, n_h - 1, horizon_var),
    q4_lrt(additive, no_antenna, n_a - 1, "antennae_grp")
  )
  
  # A one-df LRT can be slightly negative in DHGLM; use the equivalent Wald test then
  if (n_h == 2 && tests$raw_Chisq[2] < 0) {
    w <- wald_horizon(additive$data)
    tests$Chisq[2] <- w$chisq
    tests$p_value[2] <- w$p_value
    tests$test_method[2] <- "Wald fallback"
  }
  tests
}

# Run all four (many DHGLM fits so can take some time to run!!!)
q4_tests <- map_dfr(seq_len(nrow(runs)), function(i) {
  pl <- runs$plane[i]; sg <- runs$segment[i]
  run_q4(planes[[pl]]$data, segments[[sg]], planes[[pl]]$horizon) %>%
    mutate(model = paste(pl, sg), .before = 1)
})

# 5. Q5 - interactions in the head-thorax slope model ----------------------

# Prepare the data (keep frames with both head and thorax)
q5_data <- function(d) {
  d %>%
    filter(!is.na(head_angle_clean), !is.na(thorax_angle_clean)) %>%
    group_by(beetle_id) %>%
    mutate(thorax_dev = thorax_angle_clean - mean(thorax_angle_clean)) %>%
    ungroup()
}

# Same random effects and dispersion model as the main script, only the fixed effects change
fit_q5 <- function(data, horizon_var, fixed, reml = FALSE) {
  glmmTMB(
    as.formula(paste0("head_angle_clean ~ ", fixed, " + (1 + thorax_dev | beetle_id)")),
    dispformula = as.formula(paste0("~ ", horizon_var, " * antennae_grp")),
    family = gaussian, data = data, REML = reml
  )
}

# Full model, without the three-way term, without one slope interaction, and with one common slope
q5_tests <- map_dfr(names(planes), function(pl) {
  hv <- planes[[pl]]$horizon
  d <- q5_data(planes[[pl]]$data)
  full <- fit_q5(d, hv, paste0("thorax_dev * ", hv, " * antennae_grp"))
  no_three_way <- fit_q5(d, hv, paste0("(thorax_dev + ", hv, " + antennae_grp)^2"))
  no_slope_horizon <- fit_q5(d, hv, paste0("thorax_dev * antennae_grp + ", hv, " * antennae_grp"))
  no_slope_antennae <- fit_q5(d, hv, paste0("thorax_dev * ", hv, " + ", hv, " * antennae_grp"))
  common_slope <- fit_q5(d, hv, paste0("thorax_dev + ", hv, " * antennae_grp"))
  bind_rows(
    lrt(no_three_way, full, paste0("thorax_dev:", hv, ":antennae_grp"), "full model"),
    lrt(no_slope_horizon, no_three_way, paste0("thorax_dev:", hv), "no-three-way model"),
    lrt(no_slope_antennae, no_three_way, "thorax_dev:antennae_grp", "no-three-way model"),
    lrt(common_slope, no_three_way, "both thorax_dev interactions (common slope)", "no-three-way model")
  ) %>%
    mutate(model = paste(pl, "head"), .before = 1)
})

# 6. Sensitivity - do the planned contrasts change without the interactions? ---
# Same contrasts and Holm families as the main script, from the full model (as in main) and the reduced one
# Q3 and Q4 reduced = additive model, Q5 reduced = model without the three-way term (REML like the main models)

horizon_contrasts <- function(emm, horizon_var) {
  methods <- if (horizon_var == "condition") {
    list("tilted - flat" = c(-1, 1))
  } else {
    list("low - flat" = c(-1, 1, 0), "high - flat" = c(-1, 0, 1), "low - high" = c(0, 1, -1))
  }
  contrast(emm, method = methods, adjust = "none")
}

# Horizon and antenna contrasts in one table, Holm within each family
effects_table <- function(horizon, antenna) {
  one <- function(x, family) {
    s <- as.data.frame(summary(x, infer = c(TRUE, TRUE), adjust = "none"))
    g <- intersect(c("condition", "horizon_pitch", "antennae_grp"), names(s))[1]
    tibble(family = family, contrast = as.character(s$contrast), group = as.character(s[[g]]),
           estimate = s$estimate, p_raw = s$p.value) %>%
      mutate(p_holm = p.adjust(p_raw, method = "holm"))
  }
  bind_rows(one(horizon, "horizon"), one(antenna, "antenna"))
}

# Same contrasts for a Q3 model (mean orientation) or a Q5 model (thorax slope, trends = TRUE)
model_effects <- function(model, horizon_var, trends = FALSE) {
  emm <- function(specs) {
    if (trends) emtrends(model, as.formula(specs), var = "thorax_dev", component = "cond")
    else emmeans(model, as.formula(specs), component = "cond")
  }
  effects_table(
    horizon_contrasts(emm(paste0("~ ", horizon_var, " | antennae_grp")), horizon_var),
    pairs(emm(paste0("~ antennae_grp | ", horizon_var)), adjust = "none")
  )
}

# Wald contrast on the DHGLM scale model (same as the main script)
wald_contrast <- function(dat, var, level_a, level_b, scale_rhs) {
  dat$wald_term <- relevel(factor(as.character(dat[[var]])), ref = level_b)
  contrasts(dat$wald_term) <- contr.treatment(nlevels(dat$wald_term))
  fit <- fit_dhglm(dat, scale_rhs)
  if (is.null(fit)) return(tibble(sd_ratio = NA_real_, p_raw = NA_real_))
  cf <- as.data.frame(fit$phi_coeff)
  cf$term <- rownames(fit$phi_coeff)
  level <- which(levels(dat$wald_term) == level_a)
  row <- filter(cf, term == paste0("wald_term", level))
  if (nrow(row) != 1) row <- filter(cf, startsWith(term, "wald_term"))[level - 1L, , drop = FALSE]
  est <- row$Estimate[1]
  z <- est / row$`Std. Error`[1]
  tibble(sd_ratio = exp(est / 2), p_raw = 2 * pnorm(abs(z), lower.tail = FALSE))
}

# Q4 contrasts (SD ratios) from the full scale model (one term per cell) or the additive one
q4_effects <- function(data, response, horizon_var, scale) {
  dat <- spread_data(data, response, horizon_var)
  h_lv <- levels(dat$horizon)
  a_lv <- levels(dat$antennae_grp)
  h_pairs <- if (length(h_lv) == 2) list(h_lv[c(2, 1)]) else list(h_lv[c(2, 1)], h_lv[c(3, 1)], h_lv[c(2, 3)])
  a_pairs <- list(c("both", "one"), c("both", "none"), c("one", "none"))
  
  effect <- function(var, a, b, cell_a, cell_b) {
    if (scale == "cell") return(wald_contrast(dat, "scale_cell", cell_a, cell_b, "wald_term"))
    other <- c(horizon = "antennae_grp", antennae_grp = "horizon")[[var]]
    wald_contrast(dat, var, a, b, paste("wald_term +", other))
  }
  horizon <- bind_rows(lapply(a_lv, function(a) bind_rows(lapply(h_pairs, function(p) {
    mutate(effect("horizon", p[1], p[2], paste(p[1], a, sep = "__"), paste(p[2], a, sep = "__")),
           family = "horizon", contrast = paste(p[1], "-", p[2]), group = a)
  }))))
  antenna <- bind_rows(lapply(h_lv, function(h) bind_rows(lapply(a_pairs, function(p) {
    mutate(effect("antennae_grp", p[1], p[2], paste(h, p[1], sep = "__"), paste(h, p[2], sep = "__")),
           family = "antenna", contrast = paste(p[1], "-", p[2]), group = h)
  }))))
  bind_rows(horizon, antenna) %>%
    group_by(family) %>%
    mutate(p_holm = p.adjust(p_raw, method = "holm")) %>%
    ungroup() %>%
    transmute(family, contrast, group, estimate = sd_ratio, p_raw, p_holm)
}

# Put the full and reduced results side by side
compare_models <- function(full, reduced, label) {
  inner_join(full, reduced, by = c("family", "contrast", "group"), suffix = c("_full", "_reduced")) %>%
    mutate(model = label, same_call = (p_holm_full < 0.05) == (p_holm_reduced < 0.05), .before = 1) %>%
    select(model, family, contrast, group, estimate_full, estimate_reduced,
           p_holm_full, p_holm_reduced, same_call)
}

q3_sens <- map_dfr(seq_len(nrow(runs)), function(i) {
  pl <- runs$plane[i]; sg <- runs$segment[i]
  hv <- planes[[pl]]$horizon
  d <- planes[[pl]]$data
  d <- d[!is.na(d[[segments[[sg]]]]), ]
  compare_models(
    model_effects(fit_location(d, segments[[sg]], hv, " * ", reml = TRUE), hv),
    model_effects(fit_location(d, segments[[sg]], hv, " + ", reml = TRUE), hv),
    paste(pl, sg)
  )
})

# Many DHGLM refits again so can take some time to run!!!
q4_sens <- map_dfr(seq_len(nrow(runs)), function(i) {
  pl <- runs$plane[i]; sg <- runs$segment[i]
  compare_models(
    q4_effects(planes[[pl]]$data, segments[[sg]], planes[[pl]]$horizon, "cell"),
    q4_effects(planes[[pl]]$data, segments[[sg]], planes[[pl]]$horizon, "additive"),
    paste(pl, sg)
  )
})

q5_sens <- map_dfr(names(planes), function(pl) {
  hv <- planes[[pl]]$horizon
  d <- q5_data(planes[[pl]]$data)
  compare_models(
    model_effects(fit_q5(d, hv, paste0("thorax_dev * ", hv, " * antennae_grp"), reml = TRUE), hv, trends = TRUE),
    model_effects(fit_q5(d, hv, paste0("(thorax_dev + ", hv, " + antennae_grp)^2"), reml = TRUE), hv, trends = TRUE),
    paste(pl, "head")
  )
})

# Short sentence for the review
sens_line <- function(x) {
  cat(sum(!x$same_call, na.rm = TRUE), "of", nrow(x),
      "contrasts change between significant and non-significant (Holm p < 0.05) without the interaction.\n\n")
}

# 7. Review file -----------------------------------------------------------

banner <- function(title, lead = "\n") {
  cat(lead, strrep("=", 60), "\n", title, "\n", strrep("=", 60), "\n", sep = "")
}

report <- capture.output({
  cat("SCRIPT VERSION:", SCRIPT_VERSION, "\n")
  cat("glmmTMB version:", as.character(packageVersion("glmmTMB")), "\n")
  cat("dhglm version:", as.character(packageVersion("dhglm")), "\n\n")
  
  banner("DEFINITIONS", lead = "")
  cat("Same data and model structures as main-analysis.R.\n")
  cat("Each test removes interaction term(s) from a model and compares the two fits.\n")
  cat("Tests are not corrected for multiple testing, read them as sensitivity checks.\n\n")
  
  banner("Q3. IS THE HORIZON x ANTENNAE INTERACTION NEEDED FOR MEAN ORIENTATION?")
  cat("Full: Y ~ Horizon * antennae_grp + (1 | beetle_id); additive: Y ~ Horizon + antennae_grp + (1 | beetle_id).\n")
  cat("Both fitted with ML and the same cell-specific residual variance. Likelihood ratio test, AIC and BIC.\n\n")
  print_rounded(mutate(q3_tests, p_value = format_p(p_value)))
  
  banner("Q4. IS THE HORIZON x ANTENNAE INTERACTION NEEDED FOR ANGULAR SPREAD?")
  cat("DHGLM scale model, full: log(phi) ~ Horizon * antennae_grp; additive: log(phi) ~ Horizon + antennae_grp.\n")
  cat("The interaction compares full and additive. Horizon and antennae main effects are tested in the additive model.\n")
  cat("Adjusted-profile LRT; a one-df test that comes out negative uses a Wald test instead.\n\n")
  print_rounded(mutate(q4_tests, p_value = format_p(p_value)))
  
  banner("Q5. ARE THE THORAX SLOPE INTERACTIONS NEEDED FOR HEAD-THORAX COUPLING?")
  cat("Full: head ~ thorax_dev * Horizon * antennae_grp + (1 + thorax_dev | beetle_id).\n")
  cat("The three-way term is tested in the full model, the other terms in the model without the three-way term.\n")
  cat("Both fitted with ML and the same cell-specific residual variance. Likelihood ratio test, AIC and BIC.\n\n")
  print_rounded(mutate(q5_tests, p_value = format_p(p_value)))
  
  banner("SENSITIVITY. DO THE PLANNED CONTRASTS CHANGE WITHOUT THE INTERACTIONS?")
  cat("Same contrasts and Holm families as main-analysis.R, from the full model and the reduced model (REML).
")
  cat("Q3: mean orientation difference (deg). Q4: SD ratio. Q5: head-thorax slope difference.
")
  cat("Reduced = additive model (Q3, Q4) or model without the three-way term (Q5).
")
  cat("In an additive model the horizon and antenna differences are the same in every group.
")
  
  banner("Q3 contrasts, full vs additive", lead = "
")
  sens_line(q3_sens)
  print_rounded(mutate(q3_sens, across(starts_with("p_holm"), format_p)))
  
  banner("Q4 contrasts, full vs additive")
  sens_line(q4_sens)
  print_rounded(mutate(q4_sens, across(starts_with("p_holm"), format_p)))
  
  banner("Q5 contrasts, full vs no three-way term")
  sens_line(q5_sens)
  print_rounded(mutate(q5_sens, across(starts_with("p_holm"), format_p)))
})

writeLines(report, output_file)

cat("\nOutput files:\n")
cat("-", output_file, "\n")
cat("\nAnalysis complete. Review", output_file, "first.\n")