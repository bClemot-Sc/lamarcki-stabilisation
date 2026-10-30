# Head stabilisation in Kheper lamarcki - sample size figures
# Number of frames per condition, antennae group and view (treemap) and number of unique beetles (bar chart)
# Uses the raw data without the rows marked "discard", so it matches the data analysed in main-analysis.R

library(tidyverse)
library(treemapify)
library(colorspace)
library(svglite)

data_file <- "data_stabilisation_raw.csv"

dir.create("raw_figures", showWarnings = FALSE)

# 1. Data ------------------------------------------------------------------

df <- read.csv(data_file, fileEncoding = "UTF-8") %>%
  filter(is.na(comments) | !str_detect(comments, regex("discard", ignore_case = TRUE))) %>%
  mutate(
    condition = factor(condition, levels = c("flat0", "tilted10"), labels = c("Flat", "Tilted")),
    antennae = factor(antennae, levels = c("both", "right", "left", "none")),
    view = factor(recode(view, "left" = "sideL", "right" = "sideR"), levels = c("front", "sideL", "sideR"))
  )

# 2. Plot style (same as main-analysis.R) ----------------------------------

theme_fig <- function() {
  theme_classic(base_size = 11) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 11),
          legend.position = "top", legend.justification = "left")
}

save_svg <- function(plot, name, width = 8, height = 5) {
  ggsave(paste0("raw_figures/", name, ".svg"), plot,
         width = width, height = height, device = svglite, bg = "white")
}

# Roll (front) in blue, pitch (sides) in orange, as in the main figures
view_cols <- c(front = "#0072B2", sideL = "#D55E00", sideR = "#F3A074")
cond_cols <- c(Flat = "grey70", Tilted = "grey35")

# 3. Frames per condition, antennae and view -------------------------------

tree_data <- df %>%
  count(condition, antennae, view, name = "n_leaf") %>%
  mutate(n_mid = sum(n_leaf), .by = c(condition, antennae)) %>%
  mutate(n_top = sum(n_leaf), .by = condition) %>%
  mutate(
    view_label = paste0(view, "\nn=", n_leaf),
    ant_label = paste0(str_to_title(antennae), "\nn=", n_mid),
    cond_label = paste0(condition, "\nn=", n_top)
  )

p_tree <- ggplot(tree_data, aes(area = n_leaf, fill = view, label = view_label,
                                subgroup = cond_label, subgroup2 = ant_label)) +
  geom_treemap(colour = "white", size = 1) +
  geom_treemap_subgroup2_border(colour = "white", size = 2) +
  geom_treemap_subgroup_border(colour = "white", size = 4) +
  geom_treemap_text(colour = "grey15", size = 11, place = "centre", reflow = TRUE) +
  geom_treemap_subgroup2_text(colour = "grey15", size = 10, place = "topleft", fontface = "italic", alpha = 0.8) +
  geom_treemap_subgroup_text(colour = "grey15", size = 12, place = "bottomleft", fontface = "bold", grow = FALSE) +
  scale_fill_manual(values = view_cols, guide = "none") +
  labs(title = "Sample distribution",
       subtitle = paste0("Frames by horizon, antennae and view (N = ", nrow(df), ")")) +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(colour = "grey30"),
        plot.margin = margin(8, 8, 8, 8))
save_svg(p_tree, "sup1-raw-sample", width = 9, height = 6)
print(p_tree)

# 4. Unique beetles per condition and antennae group -----------------------

# Right and left antenna are stacked in the "One" group, flat and tilted bars sit side by side
beetle_data <- df %>%
  mutate(group = factor(if_else(antennae %in% c("right", "left"), "One", str_to_title(antennae)),
                        levels = c("Both", "One", "None"))) %>%
  summarise(n = n_distinct(beetle_id), .by = c(condition, group, antennae)) %>%
  mutate(x = as.numeric(group) + if_else(condition == "Flat", -0.22, 0.22),
         label = case_when(antennae == "right" ~ paste0("R = ", n),
                           antennae == "left" ~ paste0("L = ", n),
                           TRUE ~ as.character(n)))

beetle_totals <- beetle_data %>%
  filter(group == "One") %>%
  summarise(n = sum(n), .by = c(condition, x))

p_beetles <- ggplot(beetle_data, aes(x, n, fill = condition, alpha = antennae)) +
  geom_col(width = 0.4, colour = "white") +
  geom_text(data = filter(beetle_data, group == "One"), aes(label = label), alpha = 1,
            position = position_stack(vjust = 0.5), size = 3.5) +
  geom_text(data = filter(beetle_data, group != "One"), aes(label = label), alpha = 1,
            vjust = -0.4, size = 3.5) +
  geom_text(data = beetle_totals, aes(x = x, y = n, label = n), inherit.aes = FALSE, vjust = -0.4, size = 3.5) +
  scale_x_continuous(breaks = 1:3, labels = c("Both", "One", "None")) +
  scale_fill_manual(values = cond_cols, name = "Horizon") +
  scale_alpha_manual(values = c(both = 1, right = 0.55, left = 1, none = 1), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Unique beetles", subtitle = paste0("Total = ", n_distinct(df$beetle_id), " beetles"),
       x = "Antennae", y = "Beetles") +
  theme_fig() +
  theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.35))
save_svg(p_beetles, "sup1-unique-beetles", width = 6, height = 5)
print(p_beetles)