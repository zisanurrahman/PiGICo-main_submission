#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(forcats)
  library(readr)
  library(stringr)
  library(RColorBrewer)
})

# =================== Hardcoded paths ===================
dbcan_file  <- "/Users/zisan/SGI_Paper/Full_paper_clean/Data/Cazyme_analysis/dbcan_cazy_long.tsv"        # Genome | Family | Count
genome_file <- "/Users/zisan/SGI_Paper/Full_paper_clean/Data/genome_stats_out/assembly_stats_with_group.tsv"  # file | sum_len | Group
gtdb_path   <- "/Users/zisan/SGI_Paper/New_ML_Analysis_Trans/gtdb_taxonomy.csv"                           # columns: genome, domain..species
outdir      <- "/Users/zisan/SGI_Paper/Full_paper_clean/Data/Cazyme_analysis/plots_taxon"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ===== Choose GTDB rank here (domain|phylum|class|order|family|genus|species) =====
TAX_LEVEL <- "family"
TOP_N_TAXA <- 25   # how many taxa to display (by # genomes)

# ===== Group colors =====
group_colors <- c("SGI" = "#1B9E77", "PiBac" = "darkslateblue")

# =================== Helpers ===================
infer_cazy_class <- function(fam) str_extract(fam, "^(GH|GT|CE|PL|AA|CBM)")

read_gtdb_ranked <- function(path) {
  gt <- suppressMessages(readr::read_csv(path, show_col_types = FALSE))
  names(gt) <- tolower(names(gt))
  if (!"genome" %in% names(gt)) {
    stop("GTDB file must contain a 'genome' column. Found: ", paste(names(gt), collapse=", "))
  }
  ranks <- c("domain","phylum","class","order","family","genus","species")
  gt %>%
    transmute(
      Genome = as.character(genome),
      across(all_of(intersect(ranks, names(gt))), as.character)
    )
}

ensure_group_colors <- function(groups, base_map) {
  extras <- setdiff(groups, names(base_map))
  if (length(extras) == 0) return(base_map)
  extra_cols <- setNames(RColorBrewer::brewer.pal(max(3, length(extras)), "Set2")[seq_along(extras)], extras)
  c(base_map, extra_cols)
}

# =================== Load data ===================
db <- read_tsv(dbcan_file, show_col_types = FALSE) %>%
  mutate(Family = as.character(Family),
         Class  = infer_cazy_class(Family)) %>%
  filter(!is.na(Class))

gs <- read_tsv(genome_file, show_col_types = FALSE) %>%
  mutate(
    Genome    = basename(file),
    Genome    = gsub("\\.fasta$", "", Genome),
    Length_bp = sum_len
  ) %>%
  select(Genome, Length_bp, Group)

gtdb <- read_gtdb_ranked(gtdb_path)

tax_col <- TAX_LEVEL
if (!tax_col %in% names(gtdb)) stop("Chosen TAX_LEVEL '", TAX_LEVEL, "' not found in GTDB file.")

meta <- gs %>%
  left_join(gtdb %>% select(Genome, all_of(tax_col)), by = "Genome") %>%
  mutate(
    Taxon = .data[[tax_col]],
    Taxon = if_else(is.na(Taxon) | Taxon == "" | Taxon == "NA", paste0(substr(TAX_LEVEL,1,1), "__Unclassified"), Taxon),
    Taxon_label = str_remove(Taxon, "^[a-z]__")
  )

# =================== Aggregate ===================
per_genome_class <- db %>%
  group_by(Genome, Class) %>%
  summarise(n_class = sum(Count), .groups = "drop") %>%
  left_join(meta, by = "Genome") %>%
  mutate(cazy_per_10kb = if_else(Length_bp > 0, n_class / (Length_bp/10000), NA_real_))

genome_totals <- db %>%
  group_by(Genome) %>%
  summarise(total_cazy = sum(Count), .groups = "drop") %>%
  left_join(meta, by = "Genome") %>%
  mutate(total_per_10kb = if_else(Length_bp > 0, total_cazy / (Length_bp/10000), NA_real_))

tax_sizes <- meta %>%
  distinct(Genome, Taxon, Taxon_label) %>%
  count(Taxon, Taxon_label, name = "n_genomes") %>%
  arrange(desc(n_genomes), Taxon_label)

top_taxa <- tax_sizes %>%
  slice_head(n = TOP_N_TAXA) %>%
  pull(Taxon)

tax_order <- tax_sizes %>%
  filter(Taxon %in% top_taxa) %>%
  arrange(desc(n_genomes), Taxon_label) %>%
  pull(Taxon)

# ----- composition (by group) -----
tax_class_comp <- db %>%
  left_join(meta, by = "Genome") %>%
  filter(Taxon %in% top_taxa) %>%
  group_by(Taxon, Taxon_label, Group, Class) %>%
  summarise(n = sum(Count), .groups = "drop_last") %>%
  mutate(total = sum(n), prop = if_else(total > 0, n/total, 0)) %>%
  ungroup() %>%
  mutate(
    Taxon = factor(Taxon, levels = rev(tax_order)),
    Taxon_label = fct_reorder(Taxon_label, as.numeric(Taxon))
  )

# ----- composition (combined) -----
tax_class_comp_comb <- db %>%
  left_join(meta, by = "Genome") %>%
  filter(Taxon %in% top_taxa) %>%
  group_by(Taxon, Taxon_label, Class) %>%
  summarise(n = sum(Count), .groups = "drop_last") %>%
  mutate(total = sum(n), prop = if_else(total > 0, n/total, 0)) %>%
  ungroup() %>%
  mutate(
    Taxon = factor(Taxon, levels = rev(tax_order)),
    Taxon_label = fct_reorder(Taxon_label, as.numeric(Taxon))
  )

# colors
all_groups <- sort(unique(meta$Group))
cols_map   <- ensure_group_colors(all_groups, group_colors)

# =================== PLOTS ===================

# 1) CAZy class proportions by taxon — BY GROUPcazy_colors <- c(
# Define your custom CAZy class colors
cazy_colors <- c(
"GH"  = "#1B9E77",
"GT"  = "#d95f02",
"PL"  = "#7570b3",
"CE"  = "#e7298a",
"CBM" = "#66a61e",
"AA"  = "#e6ab02",
"Other" = "grey60"  # adjust based on your data
)
p_comp_tax <- tax_class_comp %>%
  ggplot(aes(y = Taxon_label, x = prop, fill = Class)) +
  geom_col(width = 0.8) +
  facet_wrap(~ Group, ncol = 1, scales = "free_y") +
  scale_x_continuous(labels = percent_format(), expand = expansion(mult = c(0,0.02))) +
  scale_fill_manual(values = cazy_colors) +
  labs(
    title = paste0("CAZyme class composition by GTDB ", TAX_LEVEL, " (top taxa) — by group"),
    x = "Proportion of CAZymes", y = NULL, fill = "CAZy Class"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.y = element_blank())
print(p_comp_tax)
ggsave(file.path(outdir, paste0("taxon_", TAX_LEVEL, "_cazy_class_proportions_by_group.png")),
       p_comp_tax, width = 10, height = max(6, length(top_taxa) * 0.25 + 2), dpi = 300)

# 1b) CAZy class proportions by taxon — COMBINED
p_comp_tax_all <- tax_class_comp_comb %>%
  ggplot(aes(y = Taxon_label, x = prop, fill = Class)) +
  geom_col(width = 0.8) +
  scale_x_continuous(labels = percent_format(), expand = expansion(mult = c(0,0.02))) +
  labs(
    #title = paste0("CAZyme class composition by GTDB ", TAX_LEVEL, " (top taxa) — combined"),
    x = "Proportion of CAZymes", y = NULL, fill = "CAZy Class"
  ) +
  theme_bw(base_size = 13) +
  scale_fill_manual(values = cazy_colors) +
  theme(panel.grid.major.y = element_blank())+
  theme(axis.text.y = element_text(face = "italic"))
print(p_comp_tax_all)
ggsave(file.path(outdir, paste0("taxon_", TAX_LEVEL, "_cazy_class_proportions_combined.png")),
       p_comp_tax_all, width = 6, height = max(6, length(top_taxa) * 0.05 + 2), dpi = 600)

# 2) Total CAZy per 10 kb — BY GROUP
tax_density <- genome_totals %>%
  filter(Taxon %in% top_taxa) %>%
  mutate(
    Taxon       = factor(Taxon, levels = rev(tax_order)),
    Taxon_label = fct_reorder(Taxon_label, as.numeric(Taxon))
  )

p_density_tax <- tax_density %>%
  ggplot(aes(y = Taxon_label, x = total_per_10kb, fill = Group, color = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.15) +
  geom_jitter(height = 0.15, width = 0, size = 1.2, alpha = 0.75) +
  facet_wrap(~ Group, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = cols_map) +
  scale_color_manual(values = cols_map) +
  labs(
    title = paste0("Total CAZyme density (per 10 kb) by GTDB ", TAX_LEVEL, " — by group"),
    x = "CAZy per 10 kb (per genome)", y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

ggsave(file.path(outdir, paste0("taxon_", TAX_LEVEL, "_cazy_total_per10kb_by_group.png")),
       p_density_tax, width = 10, height = max(6, length(top_taxa) * 0.25 + 2), dpi = 300)

# 2b) Total CAZy per 10 kb — COMBINED
p_density_tax_all <- tax_density %>%
  ggplot(aes(y = Taxon_label, x = total_per_10kb)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.2, fill = "grey80", color = "grey40") +
  geom_jitter(height = 0.15, width = 0, size = 1.2, alpha = 0.6) +
  labs(
    title = paste0("Total CAZyme density (per 10 kb) by GTDB ", TAX_LEVEL, " — combined"),
    x = "CAZy per 10 kb (per genome)", y = NULL
  ) +
  theme_bw(base_size = 11)
print(p_density_tax_all)
ggsave(file.path(outdir, paste0("taxon_", TAX_LEVEL, "_cazy_total_per10kb_combined.png")),
       p_density_tax_all, width = 9, height = max(6, length(top_taxa) * 0.25 + 2), dpi = 300)

# 3) Bubble matrix — BY GROUP
bubble_df <- per_genome_class %>%
  filter(Taxon %in% top_taxa) %>%
  group_by(Taxon, Taxon_label, Group, Class) %>%
  summarise(
    mean_per10kb = mean(cazy_per_10kb, na.rm = TRUE),
    prevalence   = mean(n_class > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Taxon       = factor(Taxon, levels = tax_order),
    Taxon_label = fct_reorder(Taxon_label, as.numeric(Taxon))
  )

p_bubble <- bubble_df %>%
  ggplot(aes(x = Class, y = Taxon_label)) +
  geom_point(aes(size = mean_per10kb, alpha = prevalence, color = Group)) +
  scale_size_continuous(name = "Mean CAZy per 10 kb") +
  scale_alpha(range = c(0.2, 0.95), name = "Prevalence") +
  scale_color_manual(values = cols_map) +
  facet_wrap(~ Group, ncol = 2, scales = "free_y") +
  labs(
    title = paste0("GTDB ", TAX_LEVEL, " × CAZy class: mean density & prevalence — by group"),
    x = "CAZy class", y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")

ggsave(file.path(outdir, paste0("taxon_", TAX_LEVEL, "_class_bubble_matrix_by_group.png")),
       p_bubble, width = 12, height = max(7, length(top_taxa) * 0.3 + 2), dpi = 300)

# 3b) Bubble matrix — COMBINED
bubble_df_all <- per_genome_class %>%
  filter(Taxon %in% top_taxa) %>%
  group_by(Taxon, Taxon_label, Class) %>%
  summarise(
    mean_per10kb = mean(cazy_per_10kb, na.rm = TRUE),
    prevalence   = mean(n_class > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Taxon       = factor(Taxon, levels = tax_order),
    Taxon_label = fct_reorder(Taxon_label, as.numeric(Taxon))
  )
min_val <- min(bubble_df_all$mean_per10kb, na.rm = TRUE)
max_val <- max(bubble_df_all$mean_per10kb, na.rm = TRUE)
p_bubble_all <- bubble_df_all %>%
  ggplot(aes(x = Class, y = Taxon_label)) +
  #geom_point(aes(size = mean_per10kb, alpha = prevalence), color = "#1B9E77") +
  geom_point(aes(size = mean_per10kb), alpha = 1, color = "#1B9E77") +
  #scale_size_continuous(name = "Mean CAZy/10kb") +
  #scale_alpha(range = c(0.2, 0.95), name = "Prevalence") +
  scale_size_area(
    name     = "Mean CAZy\nper 10kb",
    max_size = 5,   # controls largest bubble on plot + legend
    limits   = c(min_val, max_val),
    breaks   = pretty(c(min_val, max_val), n = 6)
  ) +
  labs(
    #title = paste0("GTDB ", TAX_LEVEL, " × CAZy class: mean density & prevalence — combined"),
    x = "CAZy Class", y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        legend.title = element_text(size = 11))+
  theme(axis.text.y = element_text(face = "italic"))
print(p_bubble_all)
ggsave(file.path(outdir, paste0("taxon_", TAX_LEVEL, "_class_bubble_matrix_combined.png")),
       p_bubble_all, width = 4.5, height = max(5, length(top_taxa) * 0.1 + 2), dpi = 600)

# =================== PREVALENCE PLOTS ===================
# Prevalence = fraction of genomes in a Taxon (and Group) with ≥1 gene from a CAZy class

# Stable axis orders
cls_levels <- sort(unique(per_genome_class$Class))
tax_label_order <- tax_sizes %>%
  filter(Taxon %in% top_taxa) %>%
  arrange(desc(n_genomes), Taxon_label) %>%
  pull(Taxon_label)

# ---- By group (SGI, PiBac) ----
prev_by_group <- per_genome_class %>%
  filter(Taxon %in% top_taxa) %>%
  group_by(Taxon, Taxon_label, Group, Class) %>%
  summarise(
    prevalence = mean(n_class > 0, na.rm = TRUE),
    n_genomes  = n_distinct(Genome),
    .groups = "drop"
  ) %>%
  mutate(
    Class        = factor(Class, levels = cls_levels),
    Taxon_label  = factor(Taxon_label, levels = tax_label_order)
  )

p_prev_by_group <- ggplot(prev_by_group, aes(x = Class, y = Taxon_label, fill = prevalence)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_wrap(~ Group, ncol = 2, scales = "free_y") +
  scale_fill_viridis_c(
    name   = "Prevalence\n(share of genomes)",
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    title = paste0("GTDB ", TAX_LEVEL, " × CAZy class — prevalence by group"),
    x = "CAZy Class", y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.y     = element_text(face = "italic"),
    panel.grid      = element_blank()
  )
print(p_prev_by_group)
ggsave(file.path(outdir, paste0("taxon_", TAX_LEVEL, "_prevalence_heatmap_by_group.png")),
       p_prev_by_group, width = 11, height = max(7, length(top_taxa) * 0.28 + 2), dpi = 300)

# ---- Combined (SGI + PiBac pooled) ----
prev_combined <- per_genome_class %>%
  filter(Taxon %in% top_taxa) %>%
  group_by(Taxon, Taxon_label, Class) %>%
  summarise(
    prevalence = mean(n_class > 0, na.rm = TRUE),
    n_genomes  = n_distinct(Genome),
    .groups = "drop"
  ) %>%
  mutate(
    Class        = factor(Class, levels = cls_levels),
    Taxon_label  = factor(Taxon_label, levels = tax_label_order)
  )

p_prev_combined <- ggplot(prev_combined, aes(x = Class, y = Taxon_label, fill = prevalence)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(
    name   = "Prevalence\n(share of genomes)",
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    #title = paste0("GTDB ", TAX_LEVEL, " × CAZy class — prevalence (combined)"),
    x = "CAZy Class", y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.y     = element_text(face = "italic"),
    panel.grid      = element_blank()
  )
print(p_prev_combined)
ggsave(file.path(outdir, paste0("taxon_", TAX_LEVEL, "_prevalence_heatmap_combined.png")),
       p_prev_combined, width = 5, height = max(2.1, length(top_taxa) * 0.18 + 2), dpi = 600)

# ===== Build/confirm the combined prevalence table (if not already built) =====
cls_levels <- sort(unique(per_genome_class$Class))
tax_label_order <- tax_sizes %>%
  filter(Taxon %in% top_taxa) %>%
  arrange(desc(n_genomes), Taxon_label) %>%
  pull(Taxon_label)

prev_combined <- per_genome_class %>%
  filter(Taxon %in% top_taxa) %>%
  group_by(Taxon, Taxon_label, Class) %>%
  summarise(
    prevalence = mean(n_class > 0, na.rm = TRUE),  # fraction of genomes with ≥1 gene from that class
    n_genomes  = n_distinct(Genome),
    .groups    = "drop"
  ) %>%
  mutate(
    Class       = factor(Class, levels = cls_levels),
    Taxon_label = factor(Taxon_label, levels = tax_label_order)
  )

# ===== Bubble plot (combined): size ~ prevalence =====
# By default, show legend from 0→1. If you prefer min→max of observed values, set USE_MINMAX = TRUE.
USE_MINMAX <- FALSE
min_prev <- min(prev_combined$prevalence, na.rm = TRUE)
max_prev <- max(prev_combined$prevalence, na.rm = TRUE)
lims_prev <- if (USE_MINMAX) c(min_prev, max_prev) else c(0, 1)
breaks_prev <- if (USE_MINMAX) pretty(lims_prev, n = 4) else c(0, 0.25, 0.5, 0.75, 1)

p_prev_bubble_combined <- ggplot(prev_combined, aes(x = Class, y = Taxon_label)) +
  # size encodes prevalence; single fill for clean look
  geom_point(aes(size = prevalence), shape = 21, stroke = 0.2,
             fill = "#1B9E77", color = "grey20", alpha = 0.9) +
  scale_size_area(
    name     = "Prevalence\n(share of genomes)",
    max_size = 4,                        # adjust bubble max size if needed
    limits   = lims_prev,
    breaks   = breaks_prev,
    labels   = scales::percent_format(accuracy = 1)
  ) +
  labs(
    #title = paste0("GTDB ", TAX_LEVEL, " × CAZy class — prevalence (combined)"),
    x = "CAZy Class", y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11),
    axis.text.y     = element_text(face = "italic")
  )
print(p_prev_bubble_combined)
ggsave(file.path(outdir, paste0("taxon_", TAX_LEVEL, "_prevalence_bubbles_combined.png")),
       p_prev_bubble_combined,
       width = 5.05,
       height = max(4, length(top_taxa) * 0.18 + 2),
       dpi = 600)

# =================== Exports ===================
write_tsv(tax_class_comp,       file.path(outdir, paste0("taxon_", TAX_LEVEL, "_class_proportions_by_group.tsv")))
write_tsv(tax_class_comp_comb,  file.path(outdir, paste0("taxon_", TAX_LEVEL, "_class_proportions_combined.tsv")))
write_tsv(tax_density,          file.path(outdir, paste0("taxon_", TAX_LEVEL, "_total_per10kb_by_group.tsv")))
write_tsv(bubble_df,            file.path(outdir, paste0("taxon_", TAX_LEVEL, "_class_bubble_matrix_by_group.tsv")))
write_tsv(bubble_df_all,        file.path(outdir, paste0("taxon_", TAX_LEVEL, "_class_bubble_matrix_combined.tsv")))

message("Done. Plots in: ", outdir)
