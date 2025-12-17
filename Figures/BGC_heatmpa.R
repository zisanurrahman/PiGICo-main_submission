#!/usr/bin/env Rscript
# ==========================================================
# Heatmap of BGC counts per Order × Category
# + Top category color bar (NRPS/PKS/TERPENE/RIPP/OTHER)
# ==========================================================

suppressPackageStartupMessages({
  library(tidyverse)   # dplyr, ggplot2, tidyr, readr, forcats
  library(patchwork)   # for stacking the top bar over the heatmap
})

# ---------------- user params ----------------
in_csv   <- "/Users/zisan/SGI_Paper/BGC_analysis/BiG-SCAPE/v2/output_files/2025-09-05_22-07-12_c0.3/Category_Count_per_Order.csv"
out_tiff  <- "/Users/zisan/SGI_Paper/Full_paper_clean/Figures/BGC_category_heatmap.tiff"
out_pdf  <- "/Users/zisan/SGI_Paper/Full_paper_clean/Figures/BGC_category_heatmap.pdf"
out_mat  <- "/Users/zisan/SGI_Paper/Full_paper_clean/Figures/BGC_category_heatmap_matrix.csv"

# Colors for category type bar (ALL CAPS; keep in sync with your palette)
pal_cat <- c(
  "NRPS" = "cornflowerblue",
  "PKS" = "#E6AB02",
  "TERPENE" = "#1B9E77",
  "RIPP" = "darkslateblue",
  "HYBRID" = "#D95F02",
  "OTHER" = "#BEBEBE"
)

# Heatmap fill colors (counts)
heat_low  <- "white"
heat_high <- "#1B9E77"  # steelblue-like

# Bucket order (columns)
BUCKETS <- c("NRPS", "PKS", "TERPENE", "RIPP", "HYBRID", "OTHER")

`%||%` <- function(x, y) if (is.null(x)) y else x

# Map raw Category -> 5 buckets (ALL CAPS). Anything else (including decimals/hybrids) -> OTHER
to_bucket <- function(x) {
  x <- tolower(trimws(x %||% ""))
  if (grepl("\\.", x)) return("OTHER")  # decimals => hybrid => OTHER
  if (x %in% c("nrps","pks","terpene","ripp")) toupper(x) else "OTHER"
}

# ---------------- load & prep ----------------
# Expect columns like: Order, Category, Avg_Cat, Avg_total (we use Avg_Cat as counts)
df_raw <- readr::read_csv(in_csv, show_col_types = FALSE)

df <- df_raw %>%
  mutate(Category_bucket = vapply(Category, to_bucket, character(1))) %>%
  mutate(Category_bucket = factor(Category_bucket, levels = BUCKETS)) %>%
  group_by(Order, Category_bucket) %>%
  summarise(count = sum(Avg_Cat, na.rm = TRUE), .groups = "drop") %>%
  tidyr::complete(Order, Category_bucket = factor(BUCKETS, levels = BUCKETS), fill = list(count = 0))

# Order rows (Orders) by total descending
order_totals <- df %>% group_by(Order) %>% summarise(Total = sum(count), .groups = "drop") %>% arrange(desc(Total))
df <- df %>% mutate(Order = factor(Order, levels = order_totals$Order))

# Save the matrix (wide) used in the heatmap (optional but handy)
heat_mat <- df %>%
  mutate(Category_bucket = as.character(Category_bucket)) %>%
  pivot_wider(names_from = Category_bucket, values_from = count, values_fill = 0) %>%
  arrange(factor(Order, levels = order_totals$Order))
readr::write_csv(heat_mat, out_mat)
message("Saved heatmap matrix CSV: ", out_mat)

# ---------------- plots ----------------
# ---------------- plot heatmap with x-axis categories ----------------
#Z-Score normalization
df <- df %>%
  group_by(Order) %>%
  mutate(count_z = scale(count)[,1]) %>%
  ungroup()
#df <- df %>%
#  mutate(count = ifelse(count == 0, 0.1, count))  # Avoid log10(0)

p_heat <- ggplot(df, aes(x = Category_bucket, y = Order, fill = count_z)) +
  geom_tile(color = "grey85", linewidth = 0.25) +
  scale_fill_gradient2(low = "darkslateblue", mid = "white", high = "#1B9E77",
                      midpoint = 0,
                      name = "Enrichment (Z-score)") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size = 07, color = "black"),
    legend.text      = element_text(size = 06, color = "black"),
    axis.text.x      = element_text(angle = 45, face="italic", hjust = 1, vjust = 1, 
                                    color = "black", size = 6),
    axis.text.y      = element_text(color = "black", size = 6),
    axis.title = element_text(color = "black", face = "bold", size=9),
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.7),
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    plot.margin      = margin(8, 12, 8, 12)
  ) +
  xlab("BGC Categories")+
  ylab("Orders")+
  coord_flip()
print(p_heat)

# ---------------- save ----------------
ggsave(out_tiff, p_heat, width = 6.5, height = 2, dpi = 500, bg = "white")
ggsave(out_pdf, p_heat, width = 8.5, height = 10, device = cairo_pdf)
