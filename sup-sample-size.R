### Head stabilisation in the diurnal dung beetle Kheper lamarcki
### Data sample size visualisation
### Bastien CLEMOT, Lund Vision Group, Biology Department, Lund University


# Packages and Data -------------------------------------------------------

# install.packages(c(
#   "tidyverse", "treemapify", "ggtext", "colorspace", "svglite"
# ))

library(tidyverse)
library(treemapify)
library(ggtext)
library(colorspace)
library(svglite)

# load data
df <- read.csv("data_stabilisation_raw.csv", fileEncoding = "UTF-8")


# Raw sample size treemap -------------------------------------------------

# colour palette
condition_cols <- c(
  "flat0"    = "#3F7351",
  "tilted10" = "#B88A3F"
)

# treemap data
tree_data <- df %>%
  mutate(
    condition = factor(condition, levels = c("flat0", "tilted10")),
    view      = recode(view, "left" = "sideL", "right" = "sideR")
  ) %>%
  count(condition, antennae, view, name = "n_leaf") %>%
  group_by(condition, antennae) %>%
  mutate(n_mid = sum(n_leaf)) %>%
  group_by(condition) %>%
  mutate(n_top = sum(n_leaf)) %>%
  ungroup() %>%
  group_by(condition) %>%
  mutate(antennae_rank = dense_rank(antennae)) %>%
  ungroup() %>%
  mutate(
    view_label = paste0(view, "\nn=", n_leaf),
    ant_label  = paste0(antennae, "\nn=", n_mid),
    cond_label = paste0(condition, "\nn=", n_top),
    fill_color = map2_chr(
      condition_cols[as.character(condition)],
      antennae_rank,
      ~ lighten(.x, amount = (.y - 1) * 0.10)
    )
  )

# plot
p0a <- ggplot(
  tree_data,
  aes(
    area      = n_leaf,
    fill      = fill_color,
    label     = view_label,
    subgroup  = cond_label,
    subgroup2 = ant_label
  )
) +
  geom_treemap(
    colour = "white",
    size   = 1.2
  ) +
  geom_treemap_subgroup2_border(
    colour = "white",
    size   = 1.8
  ) +
  geom_treemap_subgroup_border(
    colour = "white",
    size   = 3.0
  ) +
  
  geom_treemap_text(
    colour   = "#1a1a1a",
    size     = 18,
    place    = "centre",
    reflow   = TRUE,
    fontface = "plain"
  ) +
  geom_treemap_subgroup2_text(
    colour   = "#1a1a1a",
    alpha    = 0.75,
    size     = 15,
    place    = "topleft",
    fontface = "italic"
  ) +
  geom_treemap_subgroup_text(
    colour   = "#1a1a1a",
    size     = 18,
    place    = "bottomleft",
    fontface = "bold",
    grow     = FALSE
  ) +
  
  scale_fill_identity() +
  labs(
    title    = "Sample distribution",
    subtitle = "Nested by **condition** → **antennae** → **view**",
    caption  = paste0("Total N = ", nrow(df))
  ) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title = element_text(
      colour = "#1a1a1a",
      size   = 26,
      face   = "bold",
      margin = margin(12, 0, 4, 12)
    ),
    plot.subtitle = element_markdown(
      colour = "#444444",
      size   = 15,
      margin = margin(0, 0, 8, 12)
    ),
    plot.caption = element_text(
      colour = "#666666",
      size   = 10,
      margin = margin(6, 12, 8, 0),
      hjust  = 1
    ),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(
  filename = "raw_figures/sup1-raw-sample.svg",
  plot     = p0a,
  width    = 12,
  height   = 8,
  device   = svglite,
  bg       = "white"
)

print(p0a)


# Unique beetles per condition ------------------------------------------------------

# beetle data
beetle_data <- df %>%
  mutate(
    condition = factor(condition, levels = c("flat0", "tilted10")),
    antennae  = factor(antennae, levels = c("both", "right", "left", "none")),
    antennae_group = case_when(
      antennae %in% c("right", "left") ~ "one antenna",
      TRUE ~ as.character(antennae)
    ),
    antennae_group = factor(antennae_group, levels = c("both", "one antenna", "none")),
    antennae_part  = as.character(antennae)
  ) %>%
  summarise(
    n_beetles = n_distinct(beetle_id),
    .by = c(condition, antennae_group, antennae_part)
  ) %>%
  mutate(
    x_base = as.numeric(antennae_group),
    x_pos = case_when(
      condition == "flat0"    ~ x_base - 0.22,
      condition == "tilted10" ~ x_base + 0.22
    ),
    label = case_when(
      antennae_part == "right" ~ paste0("r=", n_beetles),
      antennae_part == "left"  ~ paste0("l=", n_beetles),
      TRUE ~ as.character(n_beetles)
    )
  )

# totals for stacked one-antenna bars
beetle_totals <- beetle_data %>%
  summarise(
    n_total = sum(n_beetles),
    .by = c(condition, antennae_group, x_pos)
  ) %>%
  filter(antennae_group == "one antenna")

# plot
p0b <- ggplot(
  beetle_data,
  aes(
    x     = x_pos,
    y     = n_beetles,
    fill  = condition,
    alpha = antennae_part
  )
) +
  geom_col(
    width  = 0.38,
    colour = "white"
  ) +
  geom_text(
    data = beetle_data %>% filter(antennae_group == "one antenna"),
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  geom_text(
    data = beetle_data %>% filter(antennae_group != "one antenna"),
    aes(label = label),
    vjust = -0.4,
    size = 4
  ) +
  geom_text(
    data = beetle_totals,
    aes(
      x     = x_pos,
      y     = n_total,
      label = n_total
    ),
    inherit.aes = FALSE,
    vjust = -0.4,
    size = 4
  ) +
  scale_x_continuous(
    breaks = c(1, 2, 3),
    labels = c("both", "one antenna", "none")
  ) +
  scale_fill_manual(values = condition_cols) +
  scale_alpha_manual(
    values = c(
      "right" = 0.55,
      "left"  = 1.00,
      "both"  = 1.00,
      "none"  = 0.85
    ),
    guide = "none"
  ) +
  labs(
    title    = "Unique beetles",
    subtitle = "Right and left antenna treatments shown together",
    x        = "Antennae treatment",
    y        = "Unique beetles",
    fill     = "Condition",
    caption  = paste0("Total unique beetles = ", n_distinct(df$beetle_id))
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      colour = "#1a1a1a",
      size   = 22,
      face   = "bold"
    ),
    plot.subtitle = element_text(
      colour = "#444444",
      size   = 13
    ),
    plot.caption = element_text(
      colour = "#666666",
      size   = 10,
      hjust  = 1
    ),
    legend.position = "top"
  )

ggsave(
  filename = "raw_figures/sup1-unique-beetles.svg",
  plot     = p0b,
  width    = 7,
  height   = 5,
  device   = svglite,
  bg       = "white"
)

print(p0b)

