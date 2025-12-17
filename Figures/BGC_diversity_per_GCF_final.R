#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tidyr)
  library(ggplot2); library(ggpubr); library(scales); library(purrr)
  library(vegan)   # diversity(), vegdist(), adonis2()
})

# ----------------- INPUTS -----------------
members_csv <- "/Users/zisan/SGI_Paper/BGC_analysis/BiG-SCAPE/v2/plots_slim_from_tsv/novel_gfc_members/novel_gcf_members.csv"
ann_path  <- "/Users/zisan/SGI_Paper/BGC_analysis/BiG-SCAPE/v2/output_files/2025-09-05_22-07-12_c0.3/record_annotations.tsv"
abund_long <- "/Users/zisan/SGI_Paper/BGC_analysis/BiG-SCAPE/v2/plots_slim_from_tsv/novel_gfc_members/abundance_join/novel_gcf_node_abundance_long.csv"
outdir    <- "/Users/zisan/SGI_Paper/BGC_analysis/BiG-SCAPE/v2/plots_slim_from_tsv/novel_gfc_members/abundance_join/diversity_outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---------- SETTINGS ----------
MIN_MEMBERS <- 5          # keep GCFs with > MIN_MEMBERS
SWAP_LABELS <- TRUE      # set TRUE to display Healthy as "PWD" and vice-versa

# ---------- LOAD ----------
members <- readr::read_csv(members_csv, show_col_types = FALSE) %>%
  select(novel_number, node_id)
long_df <- readr::read_csv(abund_long, show_col_types = FALSE) %>%
  filter(group %in% c("Healthy","PWD")) %>%
  select(group, sample, novel_number, node_id, abundance)

# keep GCFs with > MIN_MEMBERS
gcf_sizes <- members %>% count(novel_number, name = "n_members")
keep_gcfs <- gcf_sizes %>% filter(n_members > MIN_MEMBERS)

df <- long_df %>%
  semi_join(keep_gcfs, by = "novel_number") %>%
  mutate(group_disp = if (SWAP_LABELS) recode(group, Healthy="PWD", PWD="Healthy") else group)

# ---------- METRICS PER SAMPLE × GCF ----------
# total abundance
tot <- df %>%
  group_by(sample, group_disp, novel_number) %>%
  summarise(total_abundance = sum(abundance, na.rm = TRUE), .groups = "drop")

# richness (q0) inside each GCF
rich <- df %>%
  group_by(sample, group_disp, novel_number) %>%
  summarise(q0 = sum(abundance > 0, na.rm = TRUE), .groups = "drop")

mat <- tot %>% left_join(rich, by = c("sample","group_disp","novel_number"))

# scale richness (q0) to abundance axis *within each GCF* so the slope overlays
scales_tbl <- mat %>%
  group_by(novel_number) %>%
  summarise(maxA = max(total_abundance, na.rm = TRUE),
            max_q0 = max(q0, na.rm = TRUE), .groups = "drop") %>%
  mutate(maxA = ifelse(is.finite(maxA) & maxA > 0, maxA, 1),
         max_q0 = ifelse(is.finite(max_q0) & max_q0 > 0, max_q0, 1))

mat <- mat %>%
  left_join(scales_tbl, by = "novel_number") %>%
  mutate(q0_scaled = q0 / max_q0 * maxA)

# ---------- Add per-GCF p-values ----------
pval_df <- mat %>%
  group_by(novel_number) %>%
  summarise(
    pval_abund = tryCatch({
      wilcox.test(total_abundance ~ group_disp)$p.value
    }, error = function(e) NA_real_),
    
    pval_q0 = tryCatch({
      wilcox.test(q0 ~ group_disp)$p.value
    }, error = function(e) NA_real_)
  )


# ---------- MEDIANS (to draw slope lines) ----------
med_tbl <- mat %>%
  group_by(novel_number, group_disp) %>%
  summarise(
    med_abund = median(total_abundance, na.rm = TRUE),
    med_q0_scaled = median(q0_scaled, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(x = ifelse(group_disp == "Healthy", 1, 2))  # x-position for slope

# build two datasets: one for abundance median slope, one for richness (scaled) median slope
slope_abund <- med_tbl %>%
  select(novel_number, group_disp, x, y = med_abund) %>%
  arrange(novel_number, x)

slope_q0 <- med_tbl %>%
  select(novel_number, group_disp, x, y = med_q0_scaled) %>%
  arrange(novel_number, x)

# convenience: labels to drop inside each facet (who's higher by median)
dir_labels <- med_tbl %>%
  tidyr::pivot_wider(names_from = group_disp, values_from = c(med_abund, med_q0_scaled)) %>%
  mutate(
    ab_dir  = ifelse(med_abund_PWD > med_abund_Healthy, "Abundance: PWD higher", "Abundance: Healthy higher"),
    q0_dir  = ifelse(med_q0_scaled_PWD > med_q0_scaled_Healthy, "Richness: PWD higher", "Richness: Healthy higher")
  ) %>%
  select(novel_number, ab_dir, q0_dir)

ypos <- mat %>%
  group_by(novel_number) %>%
  summarise(y.max = max(total_abundance, na.rm = TRUE), .groups = "drop")

annot <- dir_labels %>%
  left_join(ypos, by = "novel_number") %>%
  mutate(x = 1.5, y1 = y.max * 1.05, y2 = y.max * 1.13)

annot <- annot %>%
  left_join(pval_df, by = "novel_number") %>%
  mutate(
    pval_abund_label = case_when(
      is.na(pval_abund) ~ "n.s.",
      pval_abund <= 0.001 ~ "***",
      pval_abund <= 0.01  ~ "**",
      pval_abund <= 0.05  ~ "*",
      TRUE ~ "n.s."
    ),
    pval_q0_label = case_when(
      is.na(pval_q0) ~ "n.s.",
      pval_q0 <= 0.001 ~ "***",
      pval_q0 <= 0.01  ~ "**",
      pval_q0 <= 0.05  ~ "*",
      TRUE ~ "n.s."
    )
  )

# ---------- PLOT ----------
gcf_lab <- function(v) paste0("GCF ", v)

p <- ggplot(mat, aes(x = group_disp, y = total_abundance)) +
  # abundance distribution
  geom_violin(aes(fill = group_disp), width = 0.9, trim = FALSE, color = NA, alpha = 0.70) +
  geom_boxplot(width = 0.15, outlier.size = 0.7, alpha = 0.95) +
  # median slope for abundance (solid line + endpoints)
  geom_line(data = slope_abund, aes(x = x, y = y, group = novel_number, color = "Median abundance"),
            linewidth = 0.5, inherit.aes = FALSE) +
  geom_point(data = slope_abund, aes(x = x, y = y, color = "Median abundance"),
             size = 1.5, inherit.aes = FALSE, shape = 16) +
  # median slope for richness q0 (scaled) (dashed line + endpoints)
  geom_line(data = slope_q0, aes(x = x, y = y, group = novel_number, color = "Median richness (q0, scaled)"),
            linewidth = 0.5, linetype = "22", inherit.aes = FALSE) +
  geom_point(data = slope_q0, aes(x = x, y = y, color = "Median richness (q0, scaled)"),
             size = 1.5, inherit.aes = FALSE, shape = 16) +
  # direction annotations
  geom_text(data = annot, aes(x = x, y = y1, label = ab_dir),
            size = 3.1, fontface = "bold", inherit.aes = FALSE) +
  geom_text(data = annot, aes(x = x, y = y2, label = q0_dir),
            size = 3.1, fontface = "italic", inherit.aes = FALSE) +
  facet_wrap(~ novel_number, scales = "free_y", labeller = labeller(novel_number = gcf_lab)) +
  scale_fill_manual(values = c(Healthy = "darkslateblue", PWD = "#1F77B4")) +
  scale_color_manual(
    name = NULL,
    values = c("Median abundance" = "grey",
               "Median richness (q0, scaled)" = "black")
  ) +
  labs(
    x = NULL,
    y = expression("Abundance (Sample"^-1* ")"),
    fill = NULL,
    #subtitle = "Dashed red slope = median richness (q0) scaled to abundance axis; Solid green = median abundance"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    axis.title.y = element_text(size = 08),
    axis.text = element_blank(),
    legend.position = "top",
    legend.box = "horizontal",
    strip.text = element_text(size = 9, face = "bold"),
    panel.grid.minor = element_blank()
  )+
  geom_text(data = annot, aes(x = 1.5, y = 100, label = pval_abund_label), inherit.aes = FALSE)
print(p)


q<-p+geom_text(data = annot, aes(x = 1.5, y = 10, label = pval_abund_label), inherit.aes = FALSE)
print(q)
ggsave(file.path(outdir, "combo_abundance_plus_richnessSlope_perGCF.png"), p,
       width = 12, height = 9, dpi = 300, bg = "white")
ggsave(file.path(outdir, "combo_abundance_plus_richnessSlope_perGCF.pdf"), p,
       width = 12, height = 9, bg = "white")

message("Saved: combo_abundance_plus_richnessSlope_perGCF.(png|pdf)")

# ============== SAVE EACH FACET AS ITS OWN PLOT ==============
unique_gcfs <- sort(unique(mat$novel_number))

make_one_plot <- function(gcf_id) {
  dat_g   <- dplyr::filter(mat,          novel_number == gcf_id)
  sab_g   <- dplyr::filter(slope_abund,  novel_number == gcf_id)
  sq0_g   <- dplyr::filter(slope_q0,     novel_number == gcf_id)
  annot_g <- dplyr::filter(annot,        novel_number == gcf_id)
  
  ggplot(dat_g, aes(x = group_disp, y = total_abundance)) +
    geom_violin(aes(fill = group_disp), width = 0.9, trim = FALSE, color = NA, alpha = 0.70) +
    geom_boxplot(width = 0.15, outlier.size = 0.7, alpha = 0.95) +
    geom_line(data = sab_g, aes(x = x, y = y, group = novel_number, color = "Median abundance"),
              linewidth = 0.5, inherit.aes = FALSE, na.rm = TRUE) +
    geom_point(data = sab_g, aes(x = x, y = y, color = "Median abundance"),
               size = 1.5, inherit.aes = FALSE, shape = 16, na.rm = TRUE) +
    geom_line(data = sq0_g, aes(x = x, y = y, group = novel_number, color = "Median richness (q0, scaled)"),
              linewidth = 0.5, linetype = "22", inherit.aes = FALSE, na.rm = TRUE) +
    geom_point(data = sq0_g, aes(x = x, y = y, color = "Median richness (q0, scaled)"),
               size = 1.5, inherit.aes = FALSE, shape = 16, na.rm = TRUE) +
    geom_text(data = annot_g, aes(x = x, y = y1, label = ab_dir),
              size = 3.1, fontface = "bold", inherit.aes = FALSE) +
    geom_text(data = annot_g, aes(x = x, y = y2, label = q0_dir),
              size = 3.1, fontface = "italic", inherit.aes = FALSE) +
    scale_fill_manual(values = c(Healthy = "darkslateblue", PWD = "#1F77B4")) +
    scale_color_manual(
      name = NULL,
      values = c("Median abundance" = "grey",
                 "Median richness (q0, scaled)" = "black")
    ) +
    labs(
      title = paste0("GCF ", gcf_id),
      x = NULL,
      y = expression("Abundance (Sample"^-1* ")"),
      fill = NULL
    ) +
    theme_minimal(base_size = 08) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      panel.grid.major = element_blank(),
      axis.text = element_blank(),
      axis.title.y = element_text(size = 6),
      legend.position = "top",
      legend.box = "horizontal",
      legend.strip.text = element_text(size = 6, face = "bold"),
      panel.grid.minor = element_blank()
    )
}

for (gcf in unique_gcfs) {
  p_one <- make_one_plot(gcf)
  png_file <- file.path(outdir, sprintf("combo_abundance_plus_richnessSlope_GCF_%03d.png", as.integer(gcf)))
  pdf_file <- file.path(outdir, sprintf("combo_abundance_plus_richnessSlope_GCF_%03d.pdf", as.integer(gcf)))
  ggsave(png_file, p_one, width = 1.5, height = 2, dpi = 600, bg = "white")
  ggsave(pdf_file, p_one, width = 4.0, height = 3.2, bg = "white")
}
message("Also saved individual GCF panels to: ", outdir)

