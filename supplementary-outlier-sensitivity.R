# Head stabilisation in Kheper lamarcki - supplementary: outlier sensitivity
# Companion to main-analysis.R: repeats the main models with all values kept and compares them to the
# outlier-removed version used in the main analysis (same models, contrasts and Holm families)
# Q1: baseline slope, Q2: baseline SD ratio, Q3: mean orientation, Q4: angular spread, Q5: head-thorax slope

library(tidyverse)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(emmeans)

# Check if dhglm is installed (functions directly called later)
if (!requireNamespace("dhglm", quietly = TRUE)) {
  stop("Package 'dhglm' is required. Install it with install.packages('dhglm') and rerun the script.")
}

set.seed(42)

# Setting the factor contrast
options(contrasts = c("contr.sum", "contr.poly"))

SCRIPT_VERSION <- "v1-outlier-sensitivity"
data_file <- "data_stabilisation_processed.csv"
output_file <- "reviews/outlier_sensitivity_review_v1.txt"

dir.create("reviews", showWarnings = FALSE)

# 1. Data (same as main-analysis.R) ----------------------------------------

# Remove a value if both outlier rules flag it, or one rule flags it and the row has a comment
flag_outlier <- function(z3, iqr, has_comment) {
  z3 <- coalesce(z3, FALSE)
  iqr <- coalesce(iqr, FALSE)
  (z3 & iqr) | ((z3 | iqr) & has_comment)
}

raw_data <- read.csv(data_file, fileEncoding = "UTF-8")

# Prepare the data, keep_all = TRUE keeps the flagged values
prep_data <- function(keep_all) {
  raw_data %>%
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
      remove_head = flag_outlier(outlier_z3_head, outlier_iqr_head, has_comment) & !keep_all,
      remove_thorax = flag_outlier(outlier_z3_thorax, outlier_iqr_thorax, has_comment) & !keep_all,
      head_angle_clean = if_else(remove_head, NA_real_, head_angle_corrected),
      thorax_angle_clean = if_else(remove_thorax, NA_real_, thorax_angle_corrected),
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
}

datasets <- list(removed = prep_data(FALSE), all_kept = prep_data(TRUE))

# Prepare data for each axis
planes_of <- function(d) {
  list(
    roll = list(data = d %>% filter(rotation_plane == "roll", !is.na(antennae_grp), !is.na(condition)) %>% droplevels(),
                horizon = "condition"),
    pitch = list(data = d %>% filter(rotation_plane == "pitch", !is.na(antennae_grp), !is.na(horizon_pitch)) %>% droplevels(),
                 horizon = "horizon_pitch")
  )
}

# Response variables
segments <- c(head = "head_angle_clean", thorax = "thorax_angle_clean")
runs <- expand_grid(plane = c("roll", "pitch"), segment = names(segments))

# 2. Helpers ---------------------------------------------------------------

# Format p values
format_p <- function(p) ifelse(is.na(p), NA_character_, formatC(p, format = "e", digits = 2))
# Round values for output
print_rounded <- function(x, digits = 3) {
  print.data.frame(as.data.frame(mutate(x, across(where(is.numeric), ~ round(.x, digits)))), row.names = FALSE)
}

# Adjust which horizon contrasts are tested
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

# Contrasts for a Q3 model (mean orientation) or a Q5 model (thorax slope, trends = TRUE)
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

# Same models as the main script
fit_location <- function(data, response, horizon_var) {
  glmmTMB(
    formula = as.formula(paste0(response, " ~ ", horizon_var, " * antennae_grp + (1 | beetle_id)")),
    dispformula = as.formula(paste0("~ ", horizon_var, " * antennae_grp")),
    family = gaussian, data = data, REML = TRUE
  )
}

q5_data <- function(d) {
  d %>%
    filter(!is.na(head_angle_clean), !is.na(thorax_angle_clean)) %>%
    group_by(beetle_id) %>%
    mutate(thorax_dev = thorax_angle_clean - mean(thorax_angle_clean)) %>%
    ungroup()
}

fit_q5 <- function(data, horizon_var) {
  glmmTMB(
    as.formula(paste0("head_angle_clean ~ thorax_dev * ", horizon_var, " * antennae_grp + (1 + thorax_dev | beetle_id)")),
    dispformula = as.formula(paste0("~ ", horizon_var, " * antennae_grp")),
    family = gaussian, data = data, REML = TRUE
  )
}

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

# Q4 contrasts (SD ratios) from the full scale model, one term per cell
q4_effects <- function(data, response, horizon_var) {
  dat <- data %>%
    filter(!is.na(.data[[response]])) %>%
    transmute(y = .data[[response]], beetle_id = factor(beetle_id),
              antennae_grp = droplevels(antennae_grp), horizon = droplevels(.data[[horizon_var]])) %>%
    mutate(scale_cell = interaction(horizon, antennae_grp, sep = "__", drop = TRUE))
  h_lv <- levels(dat$horizon)
  a_lv <- levels(dat$antennae_grp)
  h_pairs <- if (length(h_lv) == 2) list(h_lv[c(2, 1)]) else list(h_lv[c(2, 1)], h_lv[c(3, 1)], h_lv[c(2, 3)])
  a_pairs <- list(c("both", "one"), c("both", "none"), c("one", "none"))
  
  cell <- function(a, b) wald_contrast(dat, "scale_cell", a, b, "wald_term")
  horizon <- bind_rows(lapply(a_lv, function(a) bind_rows(lapply(h_pairs, function(p) {
    mutate(cell(paste(p[1], a, sep = "__"), paste(p[2], a, sep = "__")),
           family = "horizon", contrast = paste(p[1], "-", p[2]), group = a)
  }))))
  antenna <- bind_rows(lapply(h_lv, function(h) bind_rows(lapply(a_pairs, function(p) {
    mutate(cell(paste(h, p[1], sep = "__"), paste(h, p[2], sep = "__")),
           family = "antenna", contrast = paste(p[1], "-", p[2]), group = h)
  }))))
  bind_rows(horizon, antenna) %>%
    group_by(family) %>%
    mutate(p_holm = p.adjust(p_raw, method = "holm")) %>%
    ungroup() %>%
    transmute(family, contrast, group, estimate = sd_ratio, p_raw, p_holm)
}

# 3. Run all models on one version of the data -----------------------------

run_all <- function(d) {
  planes <- planes_of(d)
  
  # Q1 baseline slope and Q2 baseline SD ratio (flat horizon, both antennae)
  base <- d %>%
    filter(condition == "flat0", antennae == "both",
           !is.na(head_angle_clean), !is.na(thorax_angle_clean), !is.na(rotation_plane))
  baseline <- map_dfr(c("roll", "pitch"), function(pl) {
    b <- filter(base, rotation_plane == pl)
    x <- b %>% group_by(beetle_id) %>% mutate(thorax_dev = thorax_angle_clean - mean(thorax_angle_clean)) %>% ungroup()
    fit <- lmer(head_angle_clean ~ thorax_dev + (1 | beetle_id), data = x, REML = TRUE)
    ci <- confint(fit, parm = "thorax_dev", method = "Wald")
    sds <- b %>%
      group_by(beetle_id) %>%
      summarise(n = n(), h = sd(head_angle_clean), t = sd(thorax_angle_clean), .groups = "drop") %>%
      filter(n >= 3, is.finite(h), is.finite(t), h > 0, t > 0)
    test <- t.test(log(sds$h / sds$t), mu = 0)
    tibble(question = c("Q1 slope", "Q2 SD ratio"), model = pl,
           estimate = c(fixef(fit)[["thorax_dev"]], exp(unname(test$estimate))),
           ci_low = c(ci[1], exp(test$conf.int[1])), ci_high = c(ci[2], exp(test$conf.int[2])))
  })
  
  q3 <- map_dfr(seq_len(nrow(runs)), function(i) {
    pl <- runs$plane[i]; sg <- runs$segment[i]
    x <- planes[[pl]]$data
    x <- x[!is.na(x[[segments[[sg]]]]), ]
    model_effects(fit_location(x, segments[[sg]], planes[[pl]]$horizon), planes[[pl]]$horizon) %>%
      mutate(model = paste(pl, sg), .before = 1)
  })
  
  # Many DHGLM fits so can take some time to run!!!
  q4 <- map_dfr(seq_len(nrow(runs)), function(i) {
    pl <- runs$plane[i]; sg <- runs$segment[i]
    q4_effects(planes[[pl]]$data, segments[[sg]], planes[[pl]]$horizon) %>%
      mutate(model = paste(pl, sg), .before = 1)
  })
  
  q5 <- map_dfr(names(planes), function(pl) {
    hv <- planes[[pl]]$horizon
    model_effects(fit_q5(q5_data(planes[[pl]]$data), hv), hv, trends = TRUE) %>%
      mutate(model = paste(pl, "head"), .before = 1)
  })
  
  list(baseline = baseline, q3 = q3, q4 = q4, q5 = q5)
}

results <- map(datasets, run_all)

# 4. Compare the two versions ----------------------------------------------

# Put outlier-removed (main) and all-values-kept results side by side
compare_models <- function(x) {
  inner_join(results$removed[[x]], results$all_kept[[x]],
             by = c("model", "family", "contrast", "group"), suffix = c("_removed", "_all")) %>%
    mutate(same_call = (p_holm_removed < 0.05) == (p_holm_all < 0.05),
           change = estimate_all - estimate_removed)
}

comparisons <- list(Q3 = compare_models("q3"), Q4 = compare_models("q4"), Q5 = compare_models("q5"))

# One line per question
summary_table <- imap_dfr(comparisons, function(x, q) {
  tibble(question = q, n_contrasts = nrow(x), n_change_significance = sum(!x$same_call, na.rm = TRUE),
         max_abs_estimate_change = max(abs(x$change), na.rm = TRUE))
})

baseline_table <- inner_join(results$removed$baseline, results$all_kept$baseline,
                             by = c("question", "model"), suffix = c("_removed", "_all"))

# Contrasts that change between significant and not significant
changed <- bind_rows(comparisons, .id = "question") %>%
  filter(!same_call) %>%
  select(question, model, family, contrast, group, estimate_removed, estimate_all, p_holm_removed, p_holm_all)

# 5. Review file -----------------------------------------------------------

banner <- function(title, lead = "\n") {
  cat(lead, strrep("=", 60), "\n", title, "\n", strrep("=", 60), "\n", sep = "")
}

report <- capture.output({
  cat("SCRIPT VERSION:", SCRIPT_VERSION, "\n\n")
  
  banner("DEFINITIONS", lead = "")
  cat("Removed = main analysis (outliers removed), all = same models with every value kept.\n")
  cat("Head values removed in the main analysis:", sum(is.na(datasets$removed$head_angle_clean) & !is.na(datasets$all_kept$head_angle_clean)), "\n")
  cat("Thorax values removed in the main analysis:", sum(is.na(datasets$removed$thorax_angle_clean) & !is.na(datasets$all_kept$thorax_angle_clean)), "\n")
  cat("Contrasts, families and Holm correction as in main-analysis.R.\n")
  cat("A contrast changes if its Holm p crosses 0.05 between the two versions.\n")
  
  banner("SUMMARY PER QUESTION")
  print_rounded(summary_table)
  cat("Estimate units: Q3 degrees, Q4 SD ratio, Q5 slope difference.\n")
  
  banner("Q1 AND Q2 BASELINE")
  cat("Q1 = head-thorax slope, Q2 = head / thorax SD ratio, with 95% intervals.\n\n")
  print_rounded(baseline_table)
  
  banner("CONTRASTS THAT CHANGE SIGNIFICANCE")
  if (nrow(changed) == 0) {
    cat("None.\n")
  } else {
    print_rounded(mutate(changed, across(starts_with("p_holm"), format_p)))
  }
})

writeLines(report, output_file)

cat("\nOutput files:\n")
cat("-", output_file, "\n")
cat("\nAnalysis complete. Review", output_file, "first.\n")