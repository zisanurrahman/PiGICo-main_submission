#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(purrr)
  library(scales)
})

# ==================== INPUTS ====================
members_csv <- "/Users/zisan/SGI_Paper/BGC_analysis/BiG-SCAPE/v2/plots_slim_from_tsv/novel_gfc_members/novel_gcf_members.csv"
long_csv    <- "/Users/zisan/SGI_Paper/BGC_analysis/BiG-SCAPE/v2/plots_slim_from_tsv/novel_gfc_members/abundance_join/novel_gcf_node_abundance_long.csv"
outdir <- "/Users/zisan/SGI_Paper/BGC_analysis/BiG-SCAPE/v2/plots_slim_from_tsv/novel_gfc_members/abundance_join"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -------- load --------
members <- readr::read_csv(members_csv, show_col_types = FALSE)
long_df <- readr::read_csv(long_csv,    show_col_types = FALSE)
stopifnot(all(c("novel_number","node_id") %in% names(members)))
stopifnot(all(c("group","sample","novel_number","node_id","abundance") %in% names(long_df)))

# -------- keep GCFs with >5 members --------
gcf_sizes <- members %>% count(novel_number, name = "n_members")
keep_gcfs <- gcf_sizes %>% filter(n_members > 5)

# Per-sample summed abundance (Healthy vs PWD only)
gcf_sample <- long_df %>%
  semi_join(keep_gcfs, by = "novel_number") %>%
  filter(group %in% c("Healthy","PWD")) %>%
  group_by(group, sample, novel_number) %>%
  summarise(rel_abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  left_join(keep_gcfs, by = "novel_number")

# -------- Wilcoxon p-values per GCF --------
pval_tbl <- gcf_sample %>%
  group_by(novel_number) %>%
  summarise(
    p = tryCatch({
      x <- rel_abundance[group == "Healthy"]
      y <- rel_abundance[group == "PWD"]
      if (length(x) > 0 && length(y) > 0 && (sd(x) > 0 || sd(y) > 0)) {
        stats::wilcox.test(x, y, exact = FALSE)$p.value
      } else NA_real_
    }, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(p_label = ifelse(is.na(p), "n/a", pvalue(p, accuracy = 0.001)))

# y-positions per facet (computed within each GCF)
ypos_tbl <- gcf_sample %>%
  group_by(novel_number) %>%
  summarise(y.pos = max(rel_abundance, na.rm = TRUE) * 1.05, .groups = "drop")

pval_annot <- pval_tbl %>%
  inner_join(ypos_tbl, by = "novel_number") %>%
  transmute(
    novel_number,
    group1 = "Healthy", group2 = "PWD",
    y.position = ifelse(is.finite(y.pos), y.pos, 0),
    p, p_label
  )

# -------- plot (independent y per facet) --------
gcf_labs <- function(v) paste0("GCF ", v)

p <- ggplot(gcf_sample, aes(x = group, y = rel_abundance)) +
  geom_violin(aes(fill = group), width = 0.9, trim = FALSE, color = NA, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.size = 0.8, fatten = 1, alpha = 0.95) +
  ggpubr::stat_pvalue_manual(
    data = pval_annot,
    label = "p_label",
    xmin = "group1", xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01, step.increase = 0, size = 3.1
  ) +
  facet_wrap(~ novel_number, labeller = labeller(novel_number = gcf_labs),
             scales = "free_y") +                     # << independent y-axis
  scale_fill_manual(values = c(Healthy = "#4C78A8", PWD = "#F58518")) +
  labs(x = NULL, y = "Relative abundance (sum per GCF per sample)", fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    legend.box = "horizontal",
    strip.text = element_text(size = 9, face = "bold"),
    panel.grid.minor = element_blank()
  )
print(p)
png_path <- file.path(outdir, "violin_GCF_abundance_faceted_freeY_nMembers_gt5.png")
pdf_path <- file.path(outdir, "violin_GCF_abundance_faceted_freeY_nMembers_gt5.pdf")
ggsave(png_path, p, width = 12, height = 9, dpi = 300, bg = "white")
ggsave(pdf_path, p, width = 12, height = 9, bg = "white")

message("Saved:\n- ", png_path, "\n- ", pdf_path)
