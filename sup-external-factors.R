### External factor screening for head stabilisation
### Head stabilisation in the diurnal dung beetle Kheper lamarcki
### Bastien CLEMOT, Lund Vision Group, Biology Department, Lund University


# Packages and Data -------------------------------------------------------

# install.packages(c(
#   "tidyverse", "lme4", "lmerTest", "performance", "svglite"
# ))

library(tidyverse)
library(lme4)
library(lmerTest)
library(performance)
library(svglite)

# load data
df <- read.csv("data_stabilisation_raw.csv", fileEncoding = "UTF-8")


# 0 - Helper functions -----------------------------------------------------

# angular distance from start to target following rotation direction
path_distance <- function(start, target, rotation) {
  case_when(
    rotation == "CW"  ~ (target - start) %% 360,
    rotation == "CCW" ~ (start - target) %% 360,
    TRUE ~ NA_real_
  )
}

# significance labels
p_stars <- function(p) {
  case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    p < 0.1   ~ ".",
    TRUE      ~ ""
  )
}

# categorical association test
cat_assoc <- function(data, x, y) {
  
  tab <- table(data[[x]], data[[y]])
  test <- suppressWarnings(chisq.test(tab))
  
  tibble(
    external_var = x,
    experimental_var = y,
    chi_sq = round(unname(test$statistic), 3),
    df = unname(test$parameter),
    p_value = signif(unname(test$p.value), 3),
    signif = p_stars(unname(test$p.value))
  )
}

# numeric association test
num_assoc <- function(data, numeric_var, group_var) {
  
  model <- lm(
    as.formula(paste(numeric_var, "~", group_var)),
    data = data
  )
  
  anova(model) %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    filter(term == group_var) %>%
    transmute(
      external_var = numeric_var,
      experimental_var = group_var,
      F_value = round(`F value`, 3),
      p_value = signif(`Pr(>F)`, 3),
      signif = p_stars(`Pr(>F)`)
    )
}


# 1 - Prepare external predictors -----------------------------------------

# external data
external_data <- df %>%
  mutate(
    date = factor(date),
    
    datetime = as.POSIXct(
      paste(as.character(date), time),
      format = "%Y-%m-%d %H:%M:%OS",
      tz = "UTC"
    ),
    
    time_decimal = as.numeric(format(datetime, "%H")) +
      as.numeric(format(datetime, "%M")) / 60 +
      as.numeric(format(datetime, "%OS")) / 3600,
    
    time_decimal_c = as.numeric(scale(time_decimal, center = TRUE, scale = FALSE)),
    
    condition = factor(condition, levels = c("flat0", "tilted10")),
    antennae  = factor(antennae, levels = c("both", "right", "left", "none")),
    
    view = recode(view, "left" = "sideL", "right" = "sideR"),
    view = factor(view, levels = c("front", "sideL", "sideR")),
    
    view_group = case_when(
      view == "front" ~ "front",
      view %in% c("sideL", "sideR") ~ "side",
      TRUE ~ NA_character_
    ),
    view_group = factor(view_group, levels = c("front", "side")),
    
    head_angle_corrected = case_when(
      view == "sideL" ~ -head_angle,
      TRUE ~ head_angle
    ),
    thorax_angle_corrected = case_when(
      view == "sideL" ~ -thorax_angle,
      TRUE ~ thorax_angle
    ),
    abdomen_angle_corrected = case_when(
      view == "sideL" ~ -abdomen_angle,
      TRUE ~ abdomen_angle
    ),
    
    axis_entry = factor(axis_entry, levels = c("butt", "head")),
    
    view_angle = case_when(
      view == "front" ~ 180,
      view == "sideL" ~ 270,
      view == "sideR" ~ 90,
      TRUE ~ NA_real_
    ),
    
    entry_angle = case_when(
      axis_entry == "butt" ~ climb_angle,
      axis_entry == "head" ~ (climb_angle + 180) %% 360,
      TRUE ~ NA_real_
    ),
    
    view_distance = path_distance(entry_angle, view_angle, rotation),
    progress_raw = view_distance / angular_coverage,
    
    progress_status = case_when(
      is.na(view_distance) ~ "missing_view_distance",
      is.na(angular_coverage) ~ "missing_coverage",
      angular_coverage <= 0 ~ "zero_or_negative_coverage",
      progress_raw < 0 ~ "outside_coverage",
      progress_raw > 1 ~ "outside_coverage",
      TRUE ~ "ok"
    ),
    progress_status = factor(
      progress_status,
      levels = c(
        "ok",
        "missing_view_distance",
        "missing_coverage",
        "zero_or_negative_coverage",
        "outside_coverage"
      )
    ),
    
    progress = if_else(progress_status == "ok", progress_raw, NA_real_),
    progress_c = progress - 0.5,
    progress_c2 = progress_c^2
  )

# check whether progress could be computed from data
progress_check <- external_data %>%
  count(progress_status, name = "n") %>%
  mutate(prop = n / sum(n))

print(progress_check)

# save data 
write.csv(
  external_data,
  "data_stabilisation_processed.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# filter
external_data <- external_data %>%
  filter(
    !is.na(beetle_id),
    !is.na(date),
    !is.na(time_decimal_c),
    !is.na(condition),
    !is.na(antennae),
    !is.na(view),
    !is.na(view_group),
    !is.na(axis_entry)
  )


# 2 - External factor models -----------------------------------------------

# valid response variables by recording view
response_plan <- tribble(
  ~view_group, ~response,
  "front",     "head_angle",
  "front",     "thorax_angle",
  "side",      "head_angle_corrected",
  "side",      "thorax_angle_corrected",
  "side",      "abdomen_angle_corrected"
)

# model function
fit_external_model <- function(view_group, response, data) {
  
  model_data <- data %>%
    filter(
      view_group == .env$view_group,
      !is.na(.data[[response]]),
      !is.na(progress_c),
      !is.na(progress_c2)
    )
  
  model_formula <- as.formula(
    paste0(
      response,
      " ~ date + time_decimal_c + ",
      "axis_entry * (progress_c + progress_c2) + ",
      "(1 | beetle_id)"
    )
  )
  
  fit <- lmer(
    model_formula,
    data = model_data,
    REML = FALSE
  )
  
  anova(fit, type = 3) %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    transmute(
      response = response,
      view_group = view_group,
      term,
      F_value = round(`F value`, 3),
      p_value = signif(`Pr(>F)`, 3),
      signif = p_stars(`Pr(>F)`),
      singular = performance::check_singularity(fit),
      n_rows = nrow(model_data),
      n_beetles = n_distinct(model_data$beetle_id)
    )
}

# compact model output
external_results <- pmap_dfr(
  response_plan,
  fit_external_model,
  data = external_data
)

print(external_results)


# 3 - Confounding checks ---------------------------------------------------

# categorical external variables against experimental variables
categorical_checks <- bind_rows(
  cat_assoc(external_data, "date", "condition"),
  cat_assoc(external_data, "date", "antennae"),
  cat_assoc(external_data, "date", "view"),
  cat_assoc(external_data, "axis_entry", "condition"),
  cat_assoc(external_data, "axis_entry", "antennae"),
  cat_assoc(external_data, "axis_entry", "view"),
  cat_assoc(external_data, "view_group", "condition"),
  cat_assoc(external_data, "view_group", "antennae")
)

print(categorical_checks)

# numeric external variables against experimental variables
numeric_checks <- bind_rows(
  num_assoc(external_data, "time_decimal", "condition"),
  num_assoc(external_data, "time_decimal", "antennae"),
  num_assoc(external_data, "time_decimal", "view"),
  num_assoc(external_data, "progress", "condition"),
  num_assoc(external_data, "progress", "antennae"),
  num_assoc(external_data, "progress", "view")
)

print(numeric_checks)


# 4 - Exploratory plots ----------------------------------------------------

# reshape response data
external_long <- external_data %>%
  select(
    beetle_id,
    date,
    time_decimal,
    view_group,
    axis_entry,
    progress,
    progress_c,
    progress_c2,
    progress_status,
    head_angle,
    thorax_angle,
    abdomen_angle,
    head_angle_corrected,
    thorax_angle_corrected,
    abdomen_angle_corrected
  ) %>%
  pivot_longer(
    cols = c(
      head_angle,
      thorax_angle,
      abdomen_angle,
      head_angle_corrected,
      thorax_angle_corrected,
      abdomen_angle_corrected
    ),
    names_to = "response",
    values_to = "angle"
  ) %>%
  filter(!is.na(angle)) %>%
  semi_join(response_plan, by = c("view_group", "response"))

# angle by date
p_date <- ggplot(
  external_long,
  aes(
    x = date,
    y = angle
  )
) +
  geom_boxplot(outlier.alpha = 0.4) +
  facet_grid(response ~ view_group, scales = "free_y") +
  labs(
    title = "Angles by date",
    x = "Date",
    y = "Angle"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = "raw_figures/sup-date-effect.svg",
  plot = p_date,
  width = 10,
  height = 7,
  device = svglite,
  bg = "white"
)

# angle by time
p_time <- ggplot(
  external_long,
  aes(
    x = time_decimal,
    y = angle
  )
) +
  geom_point(alpha = 0.35) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_grid(response ~ view_group, scales = "free_y") +
  labs(
    title = "Angles by time",
    x = "Time of day",
    y = "Angle"
  ) +
  theme_classic()

ggsave(
  filename = "raw_figures/sup-time-effect.svg",
  plot = p_time,
  width = 10,
  height = 7,
  device = svglite,
  bg = "white"
)

# angle by axis entry
p_axis_entry <- ggplot(
  external_long,
  aes(
    x = axis_entry,
    y = angle,
    fill = axis_entry
  )
) +
  geom_boxplot(outlier.alpha = 0.4) +
  facet_grid(response ~ view_group, scales = "free_y") +
  labs(
    title = "Angles by axis entry",
    x = "Axis entry",
    y = "Angle"
  ) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

ggsave(
  filename = "raw_figures/sup-axis-effect.svg",
  plot = p_axis_entry,
  width = 10,
  height = 7,
  device = svglite,
  bg = "white"
)

# angle by dance progress and axis entry
p_progress <- ggplot(
  external_long %>% filter(!is.na(progress)),
  aes(
    x = progress,
    y = angle,
    colour = axis_entry
  )
) +
  geom_point(alpha = 0.35) +
  geom_smooth(
    method = "lm",
    formula = y ~ x + I(x^2),
    se = TRUE
  ) +
  facet_grid(response ~ view_group, scales = "free_y") +
  labs(
    title = "Angles by dance progress",
    subtitle = "Progress is reconstructed from view distance divided by dance coverage",
    x = "Dance progress",
    y = "Angle",
    colour = "Axis entry"
  ) +
  theme_classic() +
  theme(
    legend.position = "top"
  )

ggsave(
  filename = "raw_figures/sup-progress-effect.svg",
  plot = p_progress,
  width = 10,
  height = 7,
  device = svglite,
  bg = "white"
)

print(p_date)
print(p_time)
print(p_axis_entry)
print(p_progress)
