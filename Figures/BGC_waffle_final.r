#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# ========== USER PARAMS ==========
input_file <- "/Users/zisan/SGI_Paper/BGC_analysis/BiG-SCAPE/v2/output_files/2025-09-05_22-07-12_c0.3/record_annotations.tsv"
out_dir    <- "/Users/zisan/SGI_Paper/Full_paper_clean/Figures"
cutoff     <- "record_annotations"

out_width <- 7.0
out_height <- 7.0
dpi_tiff <- 600

# ========== WAFFLE GRID ==========
n_cols <- 10
n_rows <- 10
n_cells <- n_cols * n_rows

bucket_priority <- c("NRPS", "PKS", "TERPENE", "RIPP")
all_buckets <- c(bucket_priority, "HYBRID", "OTHER")  # All caps

class_cols <- c(
  "NRPS" = "cornflowerblue",
  "PKS" = "#E6AB02",
  "TERPENE" = "#1B9E77",
  "RIPP" = "darkslateblue",
  "HYBRID" = "#D95F02",
  "OTHER" = "#BEBEBE"
)

# ========== HELPERS ==========
`%||%` <- function(x, y) if (is.null(x)) y else x
tokenize <- function(s) {
  s <- tolower(trimws(s %||% ""))
  unlist(strsplit(s, "[/;,+:&|\\-_. ]+"), use.names = FALSE)
}

map_to_bucket <- function(class_val, category_val) {
  if (grepl("\\.", category_val)) return("HYBRID")
  
  toks <- unique(c(tokenize(class_val), tokenize(category_val)))
  if (length(toks) == 0L) return("OTHER")
  
  for (p in bucket_priority) {
    p_low <- tolower(p)
    if (any(grepl(paste0("^", p_low, "$"), toks, ignore.case = TRUE)) ||
        any(grepl(p_low, toks, fixed = TRUE))) {
      return(p)
    }
  }
  return("OTHER")
}

# ========== LOAD & CLASSIFY ==========
df <- read_tsv(input_file, show_col_types = FALSE)
df <- df %>%
  mutate(Class = as.character(Class),
         Category = as.character(Category),
         bucket = mapply(map_to_bucket, Class, Category),
         bucket = factor(bucket, levels = all_buckets))

# ========== CALCULATE WAFFLE COUNTS ==========
cnt <- df %>%
  count(bucket, name = "n", sort = TRUE) %>%
  complete(bucket = factor(all_buckets, levels = all_buckets), fill = list(n = 0)) %>%
  mutate(prop = if (sum(n) > 0) n / sum(n) else 0)

# Largest Remainder Method
raw_cells   <- cnt$prop * n_cells
cells_floor <- floor(raw_cells)
remainder   <- raw_cells - cells_floor
cells_alloc <- cells_floor
remaining   <- n_cells - sum(cells_alloc)

if (remaining > 0) {
  give_idx <- order(remainder, decreasing = TRUE)[seq_len(remaining)]
  cells_alloc[give_idx] <- cells_alloc[give_idx] + 1L
} else if (remaining < 0) {
  take_idx <- order(remainder, decreasing = FALSE)[seq_len(abs(remaining))]
  cells_alloc[take_idx] <- pmax(0L, cells_alloc[take_idx] - 1L)
}

stopifnot(sum(cells_alloc) == n_cells)

waffle_df <- cnt %>%
  mutate(cells = cells_alloc) %>%
  select(bucket, cells)

# ========== EXPAND GRID ==========
grid <- tibble(idx = 1:n_cells) %>%
  mutate(
    row = (idx - 1) %% n_rows + 1,
    col = (idx - 1) %/% n_rows + 1
  )

grid$buckets <- rep(waffle_df$bucket, times = waffle_df$cells)
grid$buckets <- factor(grid$buckets, levels = all_buckets)

label_df <- cnt %>%
  mutate(pct_label = sprintf("%.1f%%", prop * 100),
         legend_lbl = paste0(as.character(bucket), " (", n, ", ", pct_label, ")"))

label_map <- setNames(label_df$legend_lbl, label_df$bucket)
grid$bucket_lbl <- factor(label_map[as.character(grid$buckets)],
                          levels = label_map[all_buckets])

# ========== PLOT ==========
p_waffle <- ggplot(grid, aes(x = col, y = row, fill = buckets)) +
  geom_tile(width = 0.9, height = 0.9, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = class_cols, name = "Class",
                    labels = label_map[levels(grid$buckets)], drop = FALSE) +
  coord_equal(expand = FALSE) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme_classic(base_size = 11) +
  theme(
    legend.position   = "right",
    legend.key.height = unit(0.5, "cm"),
    legend.text       = element_text(size = 9),
    axis.title        = element_blank(),
    axis.text         = element_blank(),
    axis.ticks        = element_blank(),
    panel.border      = element_blank(),
    plot.margin       = margin(5, 10, 5, 10)
  )

print(p_waffle)

# ========== SAVE ==========
pdf_fp <- file.path(out_dir, sprintf("fig_bgc_class_waffle_%dx%d_%s.pdf", n_cols, n_rows, cutoff))
tif_fp <- file.path(out_dir, sprintf("fig_bgc_class_waffle_%dx%d_%s.tiff", n_cols, n_rows, cutoff))

ggsave(pdf_fp, p_waffle, width = out_width, height = out_height, device = cairo_pdf)
ggsave(tif_fp, p_waffle, width = out_width, height = out_height, dpi = dpi_tiff, compression = "lzw")

message("Saved: ", pdf_fp, " and ", tif_fp)
