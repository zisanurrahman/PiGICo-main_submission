#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# ===== Paths =====
amr_file  <- "/Users/zisan/SGI_Paper/Master_tables/master_rgi_annotations.tsv"
gtdb_path <- "/Users/zisan/SGI_Paper/New_ML_Analysis_Trans/gtdb_taxonomy.csv"
out_dir   <- "/Users/zisan/Library/CloudStorage/Dropbox/SGI_Publication/Final Clean Paper/RPKM_new/feature_matrix_from_raw/amr_by_taxon_outputs"
dir.create(out_dir, showWarnings = FALSE)

# ===== Load data =====
amr  <- read_tsv(amr_file, show_col_types = FALSE)
# === Load GTDB and rename ===
gtdb <- read_csv(gtdb_path, show_col_types = FALSE) %>%
  rename(Genome_ID = genome)
head(amr)
head(gtdb)
# ===amr# ===== Merge AMR + taxonomy =====
amr_tax <- amr %>%
  filter(!is.na(`Drug Class`)) %>%
  left_join(gtdb, by = "Genome_ID")

# ===== Function: gene richness at any taxonomic level =====
summarize_richness <- function(df, tax_level) {
  df %>%
    distinct(.data[[tax_level]], `CARD_Hit_ARO`) %>%
    count(.data[[tax_level]], name = "AMR_gene_richness") %>%
    drop_na() %>%
    write_csv(file.path(out_dir, paste0("amr_richness_by_", tax_level, ".csv")))
}

tax_levels <- c("phylum", "class", "order", "family", "genus", "species")
walk(tax_levels, ~ summarize_richness(amr_tax, .x))

# ===== Mechanism by Genus =====
mech_summary <- amr_tax %>%
  count(genus, `Resistance Mechanism`, name = "n") %>%
  drop_na(genus) %>%
  write_csv(file.path(out_dir, "amr_mechanism_by_genus.csv"))

# ===== Drug Class by Genus =====
drugclass_summary <- amr_tax %>%
  count(genus, `Drug Class`, name = "n") %>%
  drop_na(genus) %>%
  write_csv(file.path(out_dir, "amr_drugclass_by_genus.csv"))

# ===== Drug Class by Family =====
drugclass_summary <- amr_tax %>%
  count(family, `Drug Class`, name = "n") %>%
  drop_na(family) %>%
  write_csv(file.path(out_dir, "amr_drugclass_by_family.csv"))

# ===== Keyword-Based Drug Class Summary (e.g., tetracycline, beta-lactam) =====
keywords <- c("tetracycline", "beta-lactam", "macrolide", "aminoglycoside", "quinolone")
for (kw in keywords) {
  amr_tax %>%
    filter(str_detect(`Drug Class`, regex(kw, ignore_case = TRUE))) %>%
    count(genus, name = paste0(kw, "_genes")) %>%
    drop_na() %>%
    write_csv(file.path(out_dir, paste0("keyword_", kw, "_by_genus.csv")))
}

# ===== Keyword-Based Drug Class Summary (e.g., tetracycline, beta-lactam) =====
keywords <- c("tetracycline", "beta-lactam", "macrolide", "aminoglycoside", "quinolone")
for (kw in keywords) {
  amr_tax %>%
    filter(str_detect(`Drug Class`, regex(kw, ignore_case = TRUE))) %>%
    count(family, name = paste0(kw, "_genes")) %>%
    drop_na() %>%
    write_csv(file.path(out_dir, paste0("keyword_", kw, "_by_family.csv")))
}
# === Gene Richness per Genus ===
richness_by_genus <- amr_tax %>%
  filter(!is.na(`CARD_Hit_ARO`)) %>%
  distinct(Genome_ID, `CARD_Hit_ARO`, genus) %>%
  group_by(genus) %>%
  summarise(AMR_gene_richness = n_distinct(`CARD_Hit_ARO`)) %>%
  arrange(desc(AMR_gene_richness))

# === Gene Richness per family ===
richness_by_family <- amr_tax %>%
  filter(!is.na(`CARD_Hit_ARO`)) %>%
  distinct(Genome_ID, `CARD_Hit_ARO`, family) %>%
  group_by(family) %>%
  summarise(AMR_gene_richness = n_distinct(`CARD_Hit_ARO`)) %>%
  arrange(desc(AMR_gene_richness))
# === Plot ===
p<-ggplot(richness_by_family, aes(x = reorder(family, -AMR_gene_richness), y = AMR_gene_richness)) +
  geom_bar(stat = "identity", fill = "#4C72B0") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Genus",
    y = "AMR Gene Richness",
    title = "AMR Gene Richness per Genus"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)

##### Top Antibiotics Found Across All Taxa ####

# === Extract First Word from 'Drug Class' and Count ===
# Extract and clean Drug Class keywords
# Process Drug Class
drug_keywords <- amr %>%
  filter(!is.na(`Drug Class`)) %>%
  mutate(Drug_Keyword = str_extract(`Drug Class`, "^[^ ;]+") %>% str_to_lower()) %>%
  filter(!Drug_Keyword %in% c("penam", "nucleoside", "disinfecting")) %>%
  mutate(
    Drug_Keyword = case_when(
      Drug_Keyword %in% c("glycopeptide", "peptide") ~ "Peptide",
      Drug_Keyword %in% c("phenicol", "chloramphenicol") ~ "Chloramphenicol",
      TRUE ~ str_to_title(Drug_Keyword)
    )
  ) %>%
  count(Drug_Keyword, sort = TRUE) %>%
  mutate(freq = n / sum(n)) %>%
  slice_max(freq, n = 20)

# Plot
ta<-ggplot(drug_keywords, aes(x = freq, y = reorder(Drug_Keyword, freq))) +
  geom_bar(stat = "identity", fill = "#1B9E77") +
  labs(
    x = "Proportion of AMR Hits",
    y = "Resistance vs Drug Classes"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 10))+
  scale_x_continuous(labels = scales::percent_format())
print(ta)
ggsave('/Users/zisan/Library/CloudStorage/Dropbox/SGI_Publication/Final Clean Paper/RPKM_new/feature_matrix_from_raw/resistance_vs_hits_propo.png', ta, width = 5, height = 4.5, dpi = 600)


### Plot Distribution of Resistance Mechanism by family #####
library(tidyverse)

# Summarize by family and resistance mechanism
mech_by_family <- amr_tax %>%
  filter(!is.na(`Resistance Mechanism`)) %>%
  mutate(`Resistance Mechanism` = str_remove_all(`Resistance Mechanism`, regex("(?i)to\\s+antibiotic|antibiotic "))) %>%
  count(family, `Resistance Mechanism`) %>%
  group_by(family) %>%
  mutate(prop = n / sum(n))


# Plot
# Use custom palette
custom_palette <- c(
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
  "#66A61E", "#E6AB02", "#A6761D", "#666666",
  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"
)

mech <- ggplot(mech_by_family, aes(x = family, y = prop, fill = `Resistance Mechanism`)) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Family",
    y = "AMR Mechanisms Proportion",
  ) +
  scale_fill_manual(values = custom_palette)+
  theme_bw(base_size = 09) +
  theme(axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 08, angle = 45, hjust = 1, face = 'italic'),
        legend.position = "top",
        legend.key.size = unit(0.3, "lines"),
        legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent_format())

print(mech)
ggsave('/Users/zisan/Library/CloudStorage/Dropbox/SGI_Publication/Final Clean Paper/RPKM_new/feature_matrix_from_raw/mechanism_vs_hits_propo.png', mech, width = 9, height = 3.5, dpi = 600)

