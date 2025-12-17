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

# ---- FIX: slice_head must use a constant 'n' ----
top_taxa <- tax_sizes %>%
  slice_head(n = TOP_N_TAXA) %>%    # returns fewer if not enough rows; no need for min(., n())
  pull(Taxon)

tax_order <- tax_sizes %>%
  filter(Taxon %in% top_taxa) %>%
  arrange(desc(n_genomes), Taxon_label) %>%
  pull(Taxon)

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

all_groups <- sort(unique(meta$Group))
cols_map   <- ensure_group_colors(all_groups, group_colors)

# =================== Plots ===================

# 1) CAZy class proportions by taxon
p_comp_tax <- tax_class_comp %>%
  ggplot(aes(y = Taxon_label, x = prop, fill = Class)) +
  geom_col(width = 0.8) +
  facet_wrap(~ Group, ncol = 1, scales = "free_y") +
  scale_x_continuous(labels = percent_format(), expand = expansion(mult = c(0,0.02))) +
  labs(
    title = paste0("CAZyme class composition by GTDB ", TAX_LEVEL, " (top taxa)"),
    x = "Proportion of CAZymes", y = NULL, fill = "CAZy class"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.y = element_blank())
print(p_comp_tax)
ggsave(file.path(outdir, paste0("taxon_", TAX_LEVEL, "_cazy_class_proportions.png")),
       p_comp_tax, width = 10, height = max(6, length(top_taxa) * 0.25 + 2), dpi = 300)

# 2) Total CAZy per 10 kb by taxon
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
    title = paste0("Total CAZyme density (per 10 kb) by GTDB ", TAX_LEVEL),
    x = "CAZy per 10 kb (per genome)", y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")
print(p_density_tax)
ggsave(file.path(outdir, paste0("taxon_", TAX_LEVEL, "_cazy_total_per10kb.png")),
       p_density_tax, width = 10, height = max(6, length(top_taxa) * 0.25 + 2), dpi = 300)

# 3) Bubble matrix: taxon × class (size = mean per-10kb, alpha = prevalence)
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
    title = paste0("GTDB ", TAX_LEVEL, " × CAZy class: mean density & prevalence"),
    x = "CAZy class", y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")
print(p_bubble)
ggsave(file.path(outdir, paste0("taxon_", TAX_LEVEL, "_class_bubble_matrix.png")),
       p_bubble, width = 12, height = max(7, length(top_taxa) * 0.3 + 2), dpi = 300)

# =================== Exports ===================
write_tsv(tax_class_comp, file.path(outdir, paste0("taxon_", TAX_LEVEL, "_class_proportions.tsv")))
write_tsv(tax_density,    file.path(outdir, paste0("taxon_", TAX_LEVEL, "_total_per10kb.tsv")))
write_tsv(bubble_df,      file.path(outdir, paste0("taxon_", TAX_LEVEL, "_class_bubble_matrix.tsv")))

message("Done. Plots in: ", outdir)
