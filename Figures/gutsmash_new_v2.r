library(tidyverse)
library(patchwork)
library(viridis)
library(scales)

# ====== Load file ======
file_path <- "/Users/zisan/SGI_Paper/Gutsmash_analysis/gutsmash_pathway_matrix.genome_normalized_R_expanded_v2.csv"
df <- read_csv(file_path, show_col_types = FALSE)
names(df) <- tolower(names(df))  # ensure lowercase

# ====== Category color setup ======
category_colors <- c(
  "SCFA" = "darkred",
  "SCFA-other" = "darkblue",
  "Aliphatic amines" = "cornflowerblue",
  "Aromatic" = "coral1",
  "npAA" = "#F5D300",
  "E-MGC" = "darkgreen",
  "Wood-Ljungdahl" = "#007C91",
  "Other" = "grey80"
)

# ====== NEW: Arrange pathways by category ======
pathway_order_by_category <- df %>%
  distinct(pathway_reformatted, category_manual) %>%
  drop_na(category_manual) %>%
  mutate(category_manual = factor(category_manual, levels = names(category_colors))) %>%
  arrange(category_manual, pathway_reformatted) %>%
  pull(pathway_reformatted)

category_df <- df %>%
  distinct(pathway_reformatted, category_manual) %>%
  drop_na(category_manual) %>%
  mutate(
    category_manual = factor(category_manual, levels = names(category_colors)),
    pathway_reformatted = factor(pathway_reformatted, levels = pathway_order_by_category)
  )

# ---- Function to add compact inline color bar ----
category_bar_fun <- function(x_levels) {
  bar_data <- category_df %>%
    mutate(pathway_reformatted = factor(pathway_reformatted, levels = x_levels))
  
  ggplot(bar_data, aes(x = pathway_reformatted, y = "bar", fill = category_manual)) +
    geom_tile(height = 1) +
    scale_fill_manual(values = category_colors, name = "Category") +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin = margin(0, 5, -5, 5)
    )
}

# ====== Prepare main data ======
agg_base <- df %>%
  filter(presence == 1) %>%
  distinct(order, genus, pathway_reformatted) %>%
  group_by(order, pathway_reformatted) %>%
  summarise(n_genus = n(), .groups = "drop")

order_levels <- agg_base %>%
  count(order, wt = n_genus, sort = TRUE) %>%
  pull(order)

# Use the category-ordered pathway levels
agg_base <- agg_base %>%
  mutate(
    order = factor(order, levels = order_levels),
    pathway_reformatted = factor(pathway_reformatted, levels = pathway_order_by_category)
  )

# ==========
# 1️⃣ Raw counts
# ==========
p1_body <- ggplot(agg_base, aes(x = pathway_reformatted, y = order, fill = n_genus)) +
  geom_tile(color = "white") +
  scale_fill_viridis(
    option = "E",
    trans = "log10",
    name = expression(atop("Genus Count", "(log"[10]*")"))
  )+
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12, face = 'italic'),
    axis.title.y = element_text(size = 25),
    panel.grid = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
    legend.position = "none",
    axis.ticks = element_line(color="black")
  ) +
  labs(x = NULL, y = "Order")

p1 <- category_bar_fun(pathway_order_by_category) / p1_body + plot_layout(heights = c(0.05, 1))
print(p1)
ggsave("/Users/zisan/SGI_Paper/Gutsmash_analysis/order_pathway_heatmap_new.png", p1, width = 18.5, height = 10.5,dpi=600)
