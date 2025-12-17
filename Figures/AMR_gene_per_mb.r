# ===== Load libraries =====
suppressPackageStartupMessages({
  library(tidyverse)
})

# ===== Input files =====
amr_file     <- "/Users/zisan/SGI_Paper/Master_tables/master_rgi_annotations.tsv"
length_file  <- "/Users/zisan/SGI_Paper/Full_paper_clean/Data/genome_stats_out/assembly_stats_all.tsv"
out_file     <- "/Users/zisan/Library/CloudStorage/Dropbox/SGI_Publication/Final Clean Paper/RPKM_new/feature_matrix_from_raw/amr_density_by_genome.tsv"
out_plot     <- "/Users/zisan/Library/CloudStorage/Dropbox/SGI_Publication/Final Clean Paper/RPKM_new/feature_matrix_from_raw/amr_density_distribution.png"

# ===== Load data =====
amr <- read_tsv(amr_file, show_col_types = FALSE)
head(amr)
# Count AMR genes per genome
amr_gene_counts <- amr %>%
  count(Genome_ID, name = "AMR_gene_count")
str(amr_gene_counts)


# Extract genome ID from file path
lengths <- lengths %>%
  mutate(genome = file %>%
           str_remove("^/Volumes/Expansion/New_transfer/SGI_genomes/Combined_SGI_PiBac/filtered_combined_fasta_file/") %>%
           str_remove("\\.fasta$"))
colnames(lengths)

amr_gene_counts <- amr %>%
  filter(!is.na(`Drug Class`)) %>%
  distinct(Genome_ID, `CARD_Hit_ARO`) %>%  # Or Gene_ID
  count(Genome_ID, name = "AMR_gene_count")
# ===== Merge and compute AMR gene density =====
amr_density <- amr_gene_counts %>%
  left_join(lengths, by = c("Genome_ID" = "genome")) %>%
  mutate(
    genome_length_mb = sum_len / 1e6,
    AMR_genes_per_mb = (AMR_gene_count/20) / genome_length_mb
  )
head(amr_density)
# ===== Save table =====
write_tsv(amr_density, out_file)

# ===== Plot density =====
p <- ggplot(amr_density, aes(x = AMR_genes_per_mb)) +
  geom_histogram(fill = "#1B9E77", color = "white", bins = 30, alpha = 1) +
  labs(
    x = "AMR Genes per Mb",
    y = "Number of Genomes"
  ) +
  theme_bw(base_size = 11) +
  theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
# Show plot
print(p)
ggsave(out_plot, p, width = 4, height = 4, dpi = 600)


