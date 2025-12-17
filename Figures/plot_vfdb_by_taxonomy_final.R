#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(forcats)
})

# ========================
# 🔹 Input paths
# ========================
vfdb_path   <- "/Users/zisan/SGI_Paper/Full_paper_clean/Data/vfdb/combined_virulencefinder_results.csv"
gtdb_path   <- "/Users/zisan/SGI_Paper/New_ML_Analysis_Trans/gtdb_taxonomy.csv"
genome_path <- "/Users/zisan/SGI_Paper/Full_paper_clean/Data/genome_stats_out/assembly_stats_with_group.tsv"
outdir      <- "/Users/zisan/SGI_Paper/Full_paper_clean/Data/vfdb/vfdb_tax_plots"
dir.create(outdir, showWarnings = FALSE)

# ========================
# 🔹 Load and annotate
# ========================
vfdb <- read_csv(vfdb_path) |> rename(genome = Genome)
gtdb <- read_csv(gtdb_path)
genome_lengths <- read_tsv(genome_path) %>%
  mutate(genome = str_remove(basename(file), "\\.fasta$")) %>%
  select(genome, sum_len)

vfdb_annot <- vfdb %>%
  left_join(gtdb, by = "genome") %>%
  left_join(genome_lengths, by = "genome") %>%
  filter(!is.na(genus), !is.na(species), !is.na(sum_len))

# ========================
# 🔹 Summarize by genus
# ========================
vf_summary_genus <- vfdb_annot %>%
  group_by(genus, species) %>%
  mutate(n_species = n_distinct(species)) %>%
  group_by(genus) %>%
  summarise(
    total_vf       = n(),
    n_species      = unique(first(n_species)),
    norm_count     = total_vf / n_species,
    total_len      = sum(sum_len),
    vf_per_mb      = total_vf / (total_len / 1e6),
    .groups = "drop"
  ) %>%
  filter(n_species >= 1) %>%
  mutate(genus = fct_reorder(genus, norm_count))

# ========================
# 🔹 Dual Axis Bar Plot
# ========================
# ========================
# 🔹 Revised Plot: VF per Mb (primary), VF per species (secondary)
# ========================
p_dual_revised <- ggplot(vf_summary_genus, aes(x = genus)) +
  geom_col(aes(y = vf_per_mb), fill = "#1B9E77") +
  geom_line(aes(y = norm_count * max(vf_per_mb) / max(norm_count), group = 1), 
            color = "#377EB8", size = 1.2) +
  scale_y_continuous(
    name = "VF per Mb (green bars)",
    sec.axis = sec_axis(~ . * max(vf_summary_genus$norm_count) / max(vf_summary_genus$vf_per_mb),
                        name = "VF per Species (blue line)")
  ) +
  labs(
    title = "Virulence Factor Density and Normalized Abundance per Genus",
    x = "Genus"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y.left = element_text(color = "#1B9E77"),
    axis.title.y.right = element_text(color = "#377EB8")
  )

print(p_dual_revised)
ggsave(file.path(outdir, "VFDB_dual_axis_barplot_genus_REVISED.png"), p_dual_revised, width = 12, height = 7, dpi = 300)

# Assumes: vf_summary_genus has genus, vf_per_mb, total_count, n_genomes
vf_summary_genus <- vfdb_annotated %>%
  filter(!is.na(genus), !is.na(sum_len)) %>%
  group_by(genus, genome) %>%
  summarise(vf_count = n(), sum_len = unique(sum_len), .groups = "drop") %>%
  group_by(genus) %>%
  summarise(
    total_count = sum(vf_count),
    total_length = sum(sum_len),
    n_genomes = n_distinct(genome),
    vf_per_mb = total_count / (total_length / 1e6),
    vf_per_genome = total_count / n_genomes
  ) %>%
  filter(n_genomes >= 1) %>%
  mutate(genus = fct_reorder(genus, vf_per_mb, .desc = TRUE))

# Plot
p_dual_final <- ggplot(vf_summary_genus, aes(x = genus)) +
  geom_col(aes(y = vf_per_mb), fill = "#1B9E77", alpha=0.9) +
  geom_line(
    aes(y = vf_per_genome * max(vf_per_mb) / max(vf_per_genome), group = 1),
    color = "#377EB8", size = 1.2
  ) +
  geom_point(
    aes(y = vf_per_genome * max(vf_per_mb) / max(vf_per_genome)),
    color = "#377EB8", size = 2
  ) +
  scale_y_continuous(
    name = expression("VF Mb"^{-1}),
    sec.axis = sec_axis(
      ~ . * max(vf_summary_genus$vf_per_genome) / max(vf_summary_genus$vf_per_mb),
      name = expression("VF Prevalence")
    )
  ) +
  labs(
    #title = "Virulence Factor Density and Genome-Normalized Abundance per Genus",
    x = "Genus"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = 'italic'),
    axis.title.y.left = element_text(color = "#1B9E77"),
    axis.title.y.right = element_text(color = "#377EB8"),
    panel.grid = element_blank()
  )

print(p_dual_final)
ggsave(file.path(outdir, "VFDB_dual_axis_barplot_genus_final.png"), p_dual_final, width = 3.5, height = 3.5, dpi = 600)
