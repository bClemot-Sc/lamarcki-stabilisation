# Head stabilisation in Kheper lamarcki - data preparation
# Raw tracking data -> processed data used by all other scripts
# Adds recording time, sign-corrected angles and dance progress, run before outliers-check.R
# Dance progress = angular distance from dance entry to the recorded view, divided by the total dance coverage

library(tidyverse)

raw_file <- "data_stabilisation_raw.csv"
data_file <- "data_stabilisation_processed.csv"

# 1. Data ------------------------------------------------------------------

raw_data <- read.csv(raw_file, fileEncoding = "UTF-8")

# 2. Helpers ---------------------------------------------------------------

# Angular distance from start to target, following the rotation direction
path_distance <- function(start, target, rotation) {
  case_when(
    rotation == "CW" ~ (target - start) %% 360,
    rotation == "CCW" ~ (start - target) %% 360,
    TRUE ~ NA_real_
  )
}

# 3. Prepare data ----------------------------------------------------------

processed <- raw_data %>%
  mutate(
    view = recode(view, "left" = "sideL", "right" = "sideR"),
    view_group = if_else(view == "front", "front", "side"),
    
    # Time of day in decimal hours
    datetime = as.POSIXct(paste(date, time), format = "%Y-%m-%d %H:%M:%OS", tz = "UTC"),
    time_decimal = as.numeric(format(datetime, "%H")) + as.numeric(format(datetime, "%M")) / 60 +
      as.numeric(format(datetime, "%OS")) / 3600,
    time_decimal_c = time_decimal - mean(time_decimal, na.rm = TRUE),
    
    # Side views are mirrored, so the sign of sideL angles is flipped
    across(c(head_angle, thorax_angle, abdomen_angle),
           ~ if_else(view == "sideL", -.x, .x), .names = "{.col}_corrected"),
    
    # Dance progress: where in the dance the frame was recorded (0 = entry, 1 = end of the dance)
    view_angle = case_when(view == "front" ~ 180, view == "sideL" ~ 270, view == "sideR" ~ 90),
    entry_angle = case_when(axis_entry == "post" ~ climb_angle, axis_entry == "head" ~ (climb_angle + 180) %% 360),
    view_distance = path_distance(entry_angle, view_angle, rotation),
    progress_raw = view_distance / angular_coverage,
    progress_status = case_when(
      is.na(view_distance) ~ "missing_view_distance",
      is.na(angular_coverage) ~ "missing_coverage",
      angular_coverage <= 0 ~ "zero_or_negative_coverage",
      progress_raw < 0 | progress_raw > 1 ~ "outside_coverage",
      TRUE ~ "ok"
    ),
    progress = if_else(progress_status == "ok", progress_raw, NA_real_),
    progress_c = progress - 0.5,
    progress_c2 = progress_c^2
  )

# Frames with a usable progress value
processed %>%
  count(progress_status) %>%
  mutate(prop = round(n / sum(n), 3)) %>%
  print()

# 4. Save processed data ---------------------------------------------------

# This overwrites the processed file, so run outliers-check.R again afterwards to add the outlier flags
write.csv(processed, data_file, row.names = FALSE, fileEncoding = "UTF-8")