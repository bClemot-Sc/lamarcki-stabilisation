# Head stabilisation in Kheper lamarcki - supplementary: covariate sensitivity
# Companion to main-analysis.R: checks if time of day, date, dance entry axis, dance progress and antenna side
# change the sensory results (same models, contrasts and Holm families as the main script)
# Covariates enter the mean models as additive terms:
#   time = time of day (centred hours), date = recording day, dance = entry axis * (progress + progress^2)
# Main and adjusted models are fitted on the same rows (complete covariates), so only the covariates differ
# Q1: baseline slope, Q2: baseline SD ratio, Q3: mean orientation, Q4: angular spread, Q5: head-thorax slope
# Q4 has the covariates in the mean model, Q1 and Q2 use values with the covariate part removed

library(tidyverse)
library(glmmTMB)
library(lme4)
library(lmerTest)
library(emmeans)
library(svglite)

# Check if dhglm is installed (functions directly called later)
if (!requireNamespace("dhglm", quietly = TRUE)) {
  stop("Package 'dhglm' is required. Install it with install.packages('dhglm') and rerun the script.")
}

set.seed(42)

# Setting the factor contrast
options(contrasts = c("contr.sum", "contr.poly"))

SCRIPT_VERSION <- "v2-covariate-sensitivity"
data_file <- "data_stabilisation_processed.csv"
output_file <- "reviews/covariate_sensitivity_review_v2.txt"

dir.create("raw_figures", showWarnings = FALSE)
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
    condition = factor(condition, levels = c("flat0", "tilted10")),
    antennae = factor(antennae, levels = c("both", "right", "left", "none")),
    view = factor(recode(view, "left" = "sideL", "right" = "sideR"),
                  levels = c("front", "sideL", "sideR")),
    rotation_plane = factor(
      case_when(view == "front" ~ "roll", view %in% c("sideL", "sideR") ~ "pitch"),
      levels = c("roll", "pitch")
    ),
    has_comment = !is.na(comments) & str_trim(comments) != "",
    head_angle_clean = if_else(flag_outlier(outlier_z3_head, outlier_iqr_head, has_comment),
                               NA_real_, head_angle_corrected),
    thorax_angle_clean = if_else(flag_outlier(outlier_z3_thorax, outlier_iqr_thorax, has_comment),
                                 NA_real_, thorax_angle_corrected),
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
    ),
    # Covariates
    axis_entry = factor(axis_entry, levels = c("post", "head")),
    date = factor(date),
    time_c = time_decimal - mean(time_decimal, na.rm = TRUE)
  ) %>%
  filter(is.na(comments) | !str_detect(comments, regex("discard", ignore_case = TRUE)))

# Rows with every covariate, main and adjusted models both use these
covariate_data <- analysis_data %>%
  filter(!is.na(time_c), !is.na(date), !is.na(axis_entry), !is.na(progress_c))

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

# Covariate sets added to the mean models
dance <- "axis_entry * (progress_c + progress_c2)"
versions <- c(main = "", time = "time_c", date = "date", dance = dance,
              all = paste("date + time_c +", dance))
q4_versions <- c("main", "all")

# 2. Helpers ---------------------------------------------------------------

# Format p values
format_p <- function(p) ifelse(is.na(p), NA_character_, formatC(p, format = "e", digits = 2))
# Round values for output
print_rounded <- function(x, digits = 3) {
  print.data.frame(as.data.frame(mutate(x, across(where(is.numeric), ~ round(.x, digits)))), row.names = FALSE)
}
save_svg <- function(plot, name, width = 8, height = 5) {
  ggsave(paste0("raw_figures/", name, ".svg"), plot,
         width = width, height = height, device = svglite, bg = "white")
}

# Add covariate terms to a formula string
with_cov <- function(rhs, cov) if (nzchar(cov)) paste(rhs, "+", cov) else rhs

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

# Same models as the main script, with optional covariates and antenna factor
fit_location <- function(data, response, horizon_var, cov = "", antenna = "antennae_grp", reml = TRUE) {
  glmmTMB(
    formula = as.formula(paste0(response, " ~ ", with_cov(paste0(horizon_var, " * ", antenna), cov), " + (1 | beetle_id)")),
    dispformula = as.formula(paste0("~ ", horizon_var, " * antennae_grp")),
    family = gaussian, data = data, REML = reml
  )
}

q5_data <- function(d) {
  d %>%
    filter(!is.na(head_angle_clean), !is.na(thorax_angle_clean)) %>%
    group_by(beetle_id) %>%
    mutate(thorax_dev = thorax_angle_clean - mean(thorax_angle_clean)) %>%
    ungroup()
}

fit_q5 <- function(data, horizon_var, cov = "", antenna = "antennae_grp", reml = TRUE) {
  glmmTMB(
    as.formula(paste0("head_angle_clean ~ ", with_cov(paste0("thorax_dev * ", horizon_var, " * ", antenna), cov),
                      " + (1 + thorax_dev | beetle_id)")),
    dispformula = as.formula(paste0("~ ", horizon_var, " * antennae_grp")),
    family = gaussian, data = data, REML = reml
  )
}

# Double hierarchical model (same as the main script, covariates in the mean model)
fit_dhglm <- function(dat, scale_rhs, cov = "") {
  mean_model <- dhglm::DHGLMMODELING(
    Model = "mean", Link = "identity",
    LinPred = as.formula(paste0("y ~ ", with_cov("horizon * antennae_grp", cov), " + (1 | beetle_id)")),
    RandDist = "gaussian"
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
wald_contrast <- function(dat, var, level_a, level_b, scale_rhs, cov = "") {
  dat$wald_term <- relevel(factor(as.character(dat[[var]])), ref = level_b)
  contrasts(dat$wald_term) <- contr.treatment(nlevels(dat$wald_term))
  fit <- fit_dhglm(dat, scale_rhs, cov)
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

# Data for the DHGLM, one scale term per horizon x antennae cell
q4_data <- function(data, response, horizon_var, antenna = "antennae_grp") {
  data %>%
    filter(!is.na(.data[[response]])) %>%
    transmute(y = .data[[response]], beetle_id = factor(beetle_id), antennae_grp = droplevels(antennae_grp),
              horizon = droplevels(.data[[horizon_var]]), antenna = .data[[antenna]],
              date, time_c, axis_entry, progress_c, progress_c2) %>%
    mutate(scale_cell = interaction(horizon, antenna, sep = "__", drop = TRUE))
}

# Q4 contrasts (SD ratios) from the full scale model, one term per cell
q4_effects <- function(data, response, horizon_var, cov = "") {
  dat <- q4_data(data, response, horizon_var)
  h_lv <- levels(dat$horizon)
  a_lv <- levels(dat$antennae_grp)
  h_pairs <- if (length(h_lv) == 2) list(h_lv[c(2, 1)]) else list(h_lv[c(2, 1)], h_lv[c(3, 1)], h_lv[c(2, 3)])
  a_pairs <- list(c("both", "one"), c("both", "none"), c("one", "none"))
  
  cell <- function(a, b) wald_contrast(dat, "scale_cell", a, b, "wald_term", cov)
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

# 3. Do the covariates matter? ---------------------------------------------

# Likelihood ratio tests on ML fits, each block is added to the model without covariates
lrt <- function(m0, m1) {
  a <- anova(m0, m1)
  tibble(df = a$`Chi Df`[2], chisq = a$Chisq[2], p_raw = a$`Pr(>Chisq)`[2])
}

block_tests <- function(fit_fun) {
  bind_rows(
    mutate(lrt(fit_fun(""), fit_fun(versions[["time"]])), block = "time"),
    mutate(lrt(fit_fun(""), fit_fun(versions[["date"]])), block = "date"),
    mutate(lrt(fit_fun(""), fit_fun(dance)), block = "dance (entry * progress)"),
    mutate(lrt(fit_fun("axis_entry + progress_c + progress_c2"), fit_fun(dance)), block = "entry changes progress curve"),
    mutate(lrt(fit_fun("axis_entry * progress_c + progress_c2"), fit_fun(dance)), block = "entry changes curvature (entry x progress^2)")
  )
}

planes_cc <- planes_of(covariate_data)

covariate_tests <- bind_rows(
  map_dfr(seq_len(nrow(runs)), function(i) {
    pl <- runs$plane[i]; sg <- runs$segment[i]
    x <- planes_cc[[pl]]$data
    x <- x[!is.na(x[[segments[[sg]]]]), ]
    block_tests(function(cov) fit_location(x, segments[[sg]], planes_cc[[pl]]$horizon, cov, reml = FALSE)) %>%
      mutate(model = paste("Q3", pl, sg), .before = 1)
  }),
  map_dfr(names(planes_cc), function(pl) {
    x <- q5_data(planes_cc[[pl]]$data)
    block_tests(function(cov) fit_q5(x, planes_cc[[pl]]$horizon, cov, reml = FALSE)) %>%
      mutate(model = paste("Q5", pl, "head"), .before = 1)
  })
) %>%
  mutate(p_holm = p.adjust(p_raw, method = "holm"))

# Are the covariates mixed up with the experimental design? (frames, not beetles)
design_tables <- list(
  "Date by antennae" = table(covariate_data$date, covariate_data$antennae_grp),
  "Date by horizon" = table(covariate_data$date, covariate_data$condition),
  "Entry axis by horizon" = table(covariate_data$axis_entry, covariate_data$condition),
  "Entry axis by antennae" = table(covariate_data$axis_entry, covariate_data$antennae_grp)
)

# 4. Sensitivity: main versus adjusted contrasts ---------------------------

run_version <- function(cov, with_q4) {
  q3 <- map_dfr(seq_len(nrow(runs)), function(i) {
    pl <- runs$plane[i]; sg <- runs$segment[i]
    x <- planes_cc[[pl]]$data
    x <- x[!is.na(x[[segments[[sg]]]]), ]
    model_effects(fit_location(x, segments[[sg]], planes_cc[[pl]]$horizon, cov), planes_cc[[pl]]$horizon) %>%
      mutate(model = paste(pl, sg), .before = 1)
  })
  
  # Many DHGLM fits so can take some time to run!!!
  q4 <- if (with_q4) {
    map_dfr(seq_len(nrow(runs)), function(i) {
      pl <- runs$plane[i]; sg <- runs$segment[i]
      q4_effects(planes_cc[[pl]]$data, segments[[sg]], planes_cc[[pl]]$horizon, cov) %>%
        mutate(model = paste(pl, sg), .before = 1)
    })
  }
  
  q5 <- map_dfr(names(planes_cc), function(pl) {
    hv <- planes_cc[[pl]]$horizon
    model_effects(fit_q5(q5_data(planes_cc[[pl]]$data), hv, cov), hv, trends = TRUE) %>%
      mutate(model = paste(pl, "head"), .before = 1)
  })
  
  list(Q3 = q3, Q4 = q4, Q5 = q5)
}

results <- imap(versions, function(cov, nm) run_version(cov, with_q4 = nm %in% q4_versions))

# Each adjusted version against the main version, same rows
comparison <- imap_dfr(results[names(versions) != "main"], function(res, nm) {
  imap_dfr(res, function(x, q) {
    if (is.null(x)) return(NULL)
    inner_join(results$main[[q]], x, by = c("model", "family", "contrast", "group"),
               suffix = c("_main", "_adj")) %>%
      mutate(version = nm, question = q, change = estimate_adj - estimate_main,
             same_call = (p_holm_main < 0.05) == (p_holm_adj < 0.05))
  })
})

summary_table <- comparison %>%
  group_by(version, question) %>%
  summarise(n_contrasts = n(), n_significant_main = sum(p_holm_main < 0.05, na.rm = TRUE),
            n_change_significance = sum(!same_call, na.rm = TRUE),
            max_abs_estimate_change = max(abs(change), na.rm = TRUE), .groups = "drop") %>%
  mutate(version = factor(version, levels = names(versions))) %>%
  arrange(version, question)

changed <- comparison %>%
  filter(!same_call) %>%
  mutate(version = factor(version, levels = names(versions))) %>%
  arrange(version, question) %>%
  select(version, question, model, family, contrast, group, estimate_main, estimate_adj, p_holm_main, p_holm_adj)

# Q1 and Q2 baseline (flat horizon, both antennae), 74 frames are too few for the covariates in the model
# so the covariate part is taken from the full-data mean model (all covariates) and removed from the values
cov_part <- function(model, d) {
  X <- model.matrix(as.formula(paste("~", versions[["all"]])), data = d)[, -1, drop = FALSE]
  b <- glmmTMB::fixef(model)$cond[colnames(X)]
  stopifnot(!anyNA(b))
  drop(X %*% b)
}

baseline_effects <- function(adjust) {
  map_dfr(c("roll", "pitch"), function(pl) {
    full <- planes_cc[[pl]]$data
    b <- filter(full, condition == "flat0", antennae == "both",
                !is.na(head_angle_clean), !is.na(thorax_angle_clean))
    if (adjust) {
      for (sg in segments) {
        x <- full[!is.na(full[[sg]]), ]
        b[[sg]] <- b[[sg]] - cov_part(fit_location(x, sg, planes_cc[[pl]]$horizon, versions[["all"]]), b)
      }
    }
    x <- b %>% group_by(beetle_id) %>% mutate(thorax_dev = thorax_angle_clean - mean(thorax_angle_clean)) %>% ungroup()
    fit <- lmerTest::lmer(head_angle_clean ~ thorax_dev + (1 | beetle_id), data = x, REML = TRUE)
    ci <- confint(fit, parm = "thorax_dev", method = "Wald")
    sds <- b %>%
      group_by(beetle_id) %>%
      summarise(n = n(), h = sd(head_angle_clean), t = sd(thorax_angle_clean), .groups = "drop") %>%
      filter(n >= 3, is.finite(h), is.finite(t), h > 0, t > 0)
    test <- t.test(log(sds$h / sds$t), mu = 0)
    tibble(question = c("Q1 slope", "Q2 SD ratio"), model = pl,
           estimate = c(lme4::fixef(fit)[["thorax_dev"]], exp(unname(test$estimate))),
           ci_low = c(ci[1], exp(test$conf.int[1])), ci_high = c(ci[2], exp(test$conf.int[2])))
  })
}

baseline_table <- inner_join(baseline_effects(FALSE), baseline_effects(TRUE),
                             by = c("question", "model"), suffix = c("_main", "_adj"))

# 5. Antenna side (left versus right in the one-antenna group) -------------

# Right and left antenna are pooled in the main analysis, here they are split (all rows, no covariates)
planes_all <- planes_of(analysis_data)

side_lrt <- function(fit_fun) lrt(fit_fun("antennae_grp"), fit_fun("antennae"))

side_contrast <- function(model, horizon_var, trends = FALSE) {
  specs <- as.formula(paste0("~ antennae | ", horizon_var))
  emm <- if (trends) emtrends(model, specs, var = "thorax_dev", component = "cond")
  else emmeans(model, specs, component = "cond")
  s <- as.data.frame(summary(contrast(emm, method = list("right - left" = c(0, 1, -1, 0)), adjust = "none"),
                             infer = c(TRUE, TRUE)))
  tibble(group = as.character(s[[horizon_var]]), estimate = s$estimate, p_raw = s$p.value) %>%
    mutate(p_holm = p.adjust(p_raw, method = "holm"))
}

side_q3 <- map_dfr(seq_len(nrow(runs)), function(i) {
  pl <- runs$plane[i]; sg <- runs$segment[i]
  x <- planes_all[[pl]]$data
  x <- x[!is.na(x[[segments[[sg]]]]), ]
  hv <- planes_all[[pl]]$horizon
  bind_cols(
    tibble(model = paste("Q3", pl, sg)),
    side_lrt(function(a) fit_location(x, segments[[sg]], hv, antenna = a, reml = FALSE)) %>% rename(lrt_p = p_raw),
    side_contrast(fit_location(x, segments[[sg]], hv, antenna = "antennae"), hv) %>% nest(by_horizon = everything())
  )
})

side_q5 <- map_dfr(names(planes_all), function(pl) {
  x <- q5_data(planes_all[[pl]]$data)
  hv <- planes_all[[pl]]$horizon
  bind_cols(
    tibble(model = paste("Q5", pl, "head")),
    side_lrt(function(a) fit_q5(x, hv, antenna = a, reml = FALSE)) %>% rename(lrt_p = p_raw),
    side_contrast(fit_q5(x, hv, antenna = "antennae"), hv, trends = TRUE) %>% nest(by_horizon = everything())
  )
})

# Q4: right / left SD ratio from the scale model, one cell per horizon x antennae (four antenna levels)
side_q4 <- map_dfr(seq_len(nrow(runs)), function(i) {
  pl <- runs$plane[i]; sg <- runs$segment[i]
  dat <- q4_data(mutate(planes_all[[pl]]$data, side = antennae), segments[[sg]], planes_all[[pl]]$horizon, "side")
  map_dfr(levels(dat$horizon), function(h) {
    mutate(wald_contrast(dat, "scale_cell", paste0(h, "__right"), paste0(h, "__left"), "wald_term"), group = h)
  }) %>%
    mutate(model = paste("Q4", pl, sg), estimate = sd_ratio, p_holm = p.adjust(p_raw, method = "holm"), .before = 1) %>%
    select(model, group, estimate, p_raw, p_holm)
})

side_table <- bind_rows(side_q3, side_q5) %>%
  select(-by_horizon) %>%
  mutate(test = "all side terms (LRT)") %>%
  select(model, test, df, chisq, p = lrt_p)

side_contrasts <- bind_rows(
  bind_rows(side_q3, side_q5) %>%
    select(model, by_horizon) %>%
    unnest(by_horizon),
  side_q4
) %>%
  mutate(question = substr(model, 1, 2), .before = 1)

# 6. Plots (saved) ---------------------------------------------------------

theme_fig <- function() {
  theme_classic(base_size = 11) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 11),
          panel.spacing.x = unit(0.9, "lines"),
          legend.position = "top", legend.justification = "left",
          panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.35))
}

fig_cols <- c(roll_head = "#8CC8E8", roll_thorax = "#0066A6",
              pitch_head = "#F3A074", pitch_thorax = "#C94B00")

# One row per angle (head, thorax), cleaned values as used in the models
plot_data <- analysis_data %>%
  select(plane = rotation_plane, date, time_decimal, axis_entry, progress, head = head_angle_clean,
         thorax = thorax_angle_clean) %>%
  pivot_longer(c(head, thorax), names_to = "segment", values_to = "angle") %>%
  filter(!is.na(angle)) %>%
  mutate(key = paste(plane, segment, sep = "_"),
         date_label = format(as.Date(as.character(date)), "%d %b"),
         plane = factor(plane, levels = c("roll", "pitch"), labels = c("Roll", "Pitch")),
         segment = factor(segment, levels = c("head", "thorax"), labels = c("Head", "Thorax")))

plot_covariate <- function(d, x, xlab) {
  ggplot(d, aes(.data[[x]], angle, colour = key, fill = key)) +
    facet_grid(segment ~ plane, scales = "free_y") +
    scale_colour_manual(values = fig_cols, guide = "none") +
    scale_fill_manual(values = fig_cols, guide = "none") +
    labs(x = xlab, y = "Angle (deg)") +
    theme_fig()
}

# Angle by date
p_date <- plot_covariate(plot_data, "date_label", "Date") +
  geom_boxplot(alpha = 0.25, outlier.shape = NA, width = 0.6, linewidth = 0.4) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.25)
save_svg(p_date, "sup-date-effect", width = 9)

# Angle by time of day
p_time <- plot_covariate(plot_data, "time_decimal", "Time of day (h)") +
  geom_point(size = 1, alpha = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, linewidth = 0.7, alpha = 0.2)
save_svg(p_time, "sup-time-effect")

# Angle by entry axis
p_entry <- plot_covariate(plot_data, "axis_entry", "Entry axis") +
  geom_boxplot(alpha = 0.25, outlier.shape = NA, width = 0.5, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 1, alpha = 0.25)
save_svg(p_entry, "sup-axis-effect")

# Angle by dance progress, one quadratic curve per entry axis
p_progress <- plot_covariate(filter(plot_data, !is.na(progress), !is.na(axis_entry)), "progress", "Dance progress") +
  aes(shape = axis_entry, linetype = axis_entry, group = axis_entry) +
  geom_point(size = 1, alpha = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 0.7, alpha = 0.15) +
  scale_shape_manual(values = c(post = 16, head = 17), name = "Entry axis") +
  scale_linetype_manual(values = c(post = "solid", head = "dashed"), name = "Entry axis") +
  guides(shape = guide_legend(override.aes = list(colour = "grey30")),
         linetype = guide_legend(override.aes = list(colour = "grey30")))
save_svg(p_progress, "sup-progress-effect")

print(p_date)
print(p_time)
print(p_entry)
print(p_progress)

# 7. Review file -----------------------------------------------------------

banner <- function(title, lead = "\n") {
  cat(lead, strrep("=", 60), "\n", title, "\n", strrep("=", 60), "\n", sep = "")
}

report <- capture.output({
  cat("SCRIPT VERSION:", SCRIPT_VERSION, "\n\n")
  
  banner("DEFINITIONS", lead = "")
  cat("Main = main-analysis.R models, adjusted = same models with covariates added to the mean part.\n")
  cat("Both are fitted on the same rows:", nrow(covariate_data), "of", nrow(analysis_data),
      "(rows with time, date, entry axis and progress). Beetles:", n_distinct(covariate_data$beetle_id), "\n")
  cat("time = time of day (centred hours), date = recording day (fixed factor),\n")
  cat("dance = entry axis * (progress + progress^2), progress = 0 at dance start, 1 at the end.\n")
  cat("all = date + time + dance. Q4 is run for", paste(q4_versions, collapse = " and "), "only (slow).\n")
  cat("Contrasts, families and Holm correction as in main-analysis.R.\n")
  cat("A contrast changes if its Holm p crosses 0.05 between main and adjusted.\n")
  
  banner("SUMMARY PER VERSION AND QUESTION")
  print_rounded(summary_table)
  cat("Estimate units: Q3 degrees, Q4 SD ratio, Q5 slope difference.\n")
  
  banner("CONTRASTS THAT CHANGE SIGNIFICANCE")
  if (nrow(changed) == 0) {
    cat("None.\n")
  } else {
    print_rounded(mutate(changed, across(starts_with("p_holm"), format_p)))
  }
  
  banner("Q1 AND Q2 BASELINE")
  cat("Q1 = head-thorax slope, Q2 = head / thorax SD ratio, with 95% intervals.\n")
  cat("Adjusted = values with the date, time and dance terms (estimated in the full-data models) removed.\n\n")
  print_rounded(baseline_table)
  
  banner("DO THE COVARIATES MATTER?")
  cat("Likelihood ratio test of each block added to the model without covariates (ML fits).\n")
  cat("Q3 = mean orientation models, Q5 = head model with thorax deviation. Holm over all rows.\n\n")
  print_rounded(mutate(covariate_tests, across(c(p_raw, p_holm), format_p)))
  
  banner("ARE THE COVARIATES MIXED UP WITH THE DESIGN?")
  cat("Number of frames (not beetles), rows used in the adjusted models.\n")
  iwalk(design_tables, function(tab, nm) {
    cat("\n", nm, "\n", sep = "")
    print(tab)
  })
  
  banner("ANTENNA SIDE: RIGHT VERSUS LEFT")
  cat("Right and left antenna split in the one-antenna group, all rows, no covariates.\n")
  cat("LRT = does separating right and left improve the model, Q3 and Q5 on the mean, Q4 on the scale.\n\n")
  print_rounded(mutate(side_table, p = format_p(p)))
  cat("\nRight - left within each horizon (Holm over horizons). Q3 degrees, Q4 SD ratio, Q5 slope.\n\n")
  print_rounded(mutate(side_contrasts, across(c(p_raw, p_holm), format_p)))
})

writeLines(report, output_file)

cat("\nOutput files:\n")
cat("-", output_file, "\n")
cat("- raw_figures/sup-date-effect.svg, sup-time-effect.svg, sup-axis-effect.svg, sup-progress-effect.svg\n")
cat("\nAnalysis complete. Review", output_file, "first.\n")