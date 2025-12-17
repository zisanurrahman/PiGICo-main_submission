#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

# ---------------- CONFIG ----------------
per_genome_tsv <- "/Users/zisan/Library/CloudStorage/Sync/SGI_guilds/SGI_guilds/Scripts/Final_scripts/AMP_density/amp_density_per_genome.tsv"
taxonomy_tsv   <- "/Volumes/Expansion/New_transfer/SGI_genomes/Combined_SGI_PiBac/GTDB-Tk_out_on_filtered_SGI_plus_PiBac/gtdbtk.bac120.summary.tsv"
outdir         <- "/Users/zisan/Library/CloudStorage/Sync/SGI_guilds/SGI_guilds/Scripts/Final_scripts/AMP_density/AMP_density_plots_R"
taxon_level    <- "Genus"   # choose: "Genus", "Phylum", "Family" or "Species"
min_genomes    <- 3         # groups need ≥ this many genomes for violin/box
topN_bubble    <- 10        # top N groups for bubble
# -----------------------------------------

dir.create(file.path(outdir, "figures"), showWarnings = FALSE, recursive = TRUE)

# ---------- helpers ----------
parse_rank <- function(class_str, prefix){
  m <- str_match(class_str %||% "", paste0(prefix, "([^;]+)"))
  ifelse(is.na(m[,2]), NA_character_, m[,2])
}
`%||%` <- function(x, y) if (is.null(x)) y else x

# ---------- load data ----------
dens <- read_tsv(per_genome_tsv, show_col_types = FALSE) %>%
  mutate(rho_AMP = as.numeric(rho_AMP),
         n_cAMPs = as.numeric(n_cAMPs),
         Length_bp = as.numeric(Length_bp))

tax <- read_tsv(taxonomy_tsv, show_col_types = FALSE)
ucol <- if ("user_genome" %in% names(tax)) "user_genome" else names(tax)[1]
if ("classification" %in% names(tax)){
  tax <- tax %>%
    mutate(
      Phylum  = parse_rank(classification, "p__"),
      Family  = parse_rank(classification, "f__"),
      Genus   = parse_rank(classification, "g__"),
      Species = parse_rank(classification, "s__")
    )
}
tax <- tax %>% select(Genome = all_of(ucol), any_of(c("Phylum","Family","Genus","Species")))

df <- dens %>% left_join(tax, by = "Genome")

# ---------- plots ----------
## Bubble
bubble_data <- df %>%
  filter(!is.na(.data[[taxon_level]]), !is.na(rho_AMP)) %>%
  group_by(.data[[taxon_level]]) %>%
  summarise(
    mean_rho = mean(rho_AMP, na.rm=TRUE),
    sd_rho = sd(rho_AMP, na.rm=TRUE),
    n_genomes = sum(!is.na(rho_AMP)),
    prevalence = mean(replace_na(n_cAMPs,0) >= 1),
    .groups="drop"
  ) %>%
  mutate(
    se_rho = sd_rho / sqrt(pmax(n_genomes,1)),
    ci95_lo = mean_rho - 1.96*se_rho,
    ci95_hi = mean_rho + 1.96*se_rho
  ) %>%
  arrange(desc(mean_rho)) %>%
  slice_head(n = topN_bubble)

p_bubble <- ggplot(bubble_data, aes(x=prevalence, y=mean_rho)) +
  geom_point(aes(size=n_genomes), alpha=0.85, stroke=0.5, colour="steelblue") +
  geom_errorbar(aes(ymin=ci95_lo, ymax=ci95_hi), width=0, alpha=0.7) +
  ggrepel::geom_text_repel(aes(label=.data[[taxon_level]]), max.overlaps=12, size=3) +
  scale_size_continuous(range=c(3,16)) +
  scale_x_continuous(limits=c(-0.02,1.02), labels=percent_format(accuracy=1)) +
  labs(x="Share of genomes with ≥1 cAMP (prevalence)", y="Mean AMP density (ρAMP)",
       title=paste("AMP density by", taxon_level, "— bubble plot")) +
  theme_bw(base_size=12)
print(p_bubble)

ggsave(file.path(outdir,"figures",paste0("bubble_",tolower(taxon_level),"_density.png")), p_bubble, width=10, height=7, dpi=300)
ggsave(file.path(outdir,"figures",paste0("bubble_",tolower(taxon_level),"_density.pdf")), p_bubble, width=10, height=7)

## Violin
keep <- df %>% count(.data[[taxon_level]]) %>% filter(n >= min_genomes) %>% pull(.data[[taxon_level]])

tax_col <- taxon_level      # "Genus", "Phylum", or "Species"
min_genomes <- 3            # keep groups with ≥ this many genomes

# 1) pick top 10 groups by median rho_AMP (after min_genomes filter)
top10_tbl <- df %>%
  filter(!is.na(.data[[tax_col]]), !is.na(rho_AMP)) %>%
  group_by(.data[[tax_col]]) %>%
  filter(n() >= min_genomes) %>%
  summarise(med = median(rho_AMP, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med)) %>%
  slice_head(n = 10)

top10 <- top10_tbl[[tax_col]]
# 2) subset data to top 10 and make clean labels (strip suffix after "_")
vdf_top10 <- df %>%
  filter(.data[[tax_col]] %in% top10) %>%
  mutate(clean_label = str_replace(.data[[tax_col]], "_.*", ""))

# 3) order by median again (ensures factor order matches top10 ranking)
order_top10 <- vdf_top10 %>%
  group_by(.data[[tax_col]]) %>%
  summarise(med = median(rho_AMP, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med)) %>%
  pull(.data[[tax_col]])

order_clean <- str_replace(order_top10, "_.*", "")

#For other taxon_levels
#vdf <- df %>% filter(.data[[taxon_level]] %in% keep)

#order <- vdf %>% group_by(.data[[taxon_level]]) %>% summarise(med=median(rho_AMP,na.rm=TRUE),.groups="drop") %>%
#  arrange(desc(med)) %>% pull(.data[[taxon_level]])
#head(order)

p_violin <- ggplot(vdf_top10, aes(x = factor(clean_label, levels = order),
                            y = rho_AMP)) +
  geom_violin(fill = "grey", color = "grey20", scale = "width", trim = TRUE) +
  stat_summary(fun = "median", geom = "point", size = 1.6, color = "darkred") +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1, color = "black") +
#  labs(x = paste0(taxon_level," (≥",min_genomes," genomes)"),
#       y = "AMP density (ρAMP)",
#       title = paste("Per-genome AMP density by", taxon_level,"— violin")) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )+
  ylab ("AMP Density (ρAMP)")+
  xlab ("Genus (at least 3 genomes)" )

print(p_violin)

ggsave(file.path(outdir,"figures",paste0("violin_",tolower(taxon_level),"_density.tiff")), p_violin, width=5, height=5, dpi=500, compression = 'lzw')
ggsave(file.path(outdir,"figures",paste0("violin_",tolower(taxon_level),"_density.pdf")), p_violin, width=max(10,0.35*length(order)), height=6)


####Plot top 10 #####
tax_col <- taxon_level      # "Genus", "Phylum", "Family" or "Species"
min_genomes <- 3            # keep groups with ≥ this many genomes
# 1) pick top 10 groups by median rho_AMP (after min_genomes filter)
top10_tbl <- df %>%
  filter(!is.na(.data[[tax_col]]), !is.na(rho_AMP)) %>%
  group_by(.data[[tax_col]]) %>%
  filter(n() >= min_genomes) %>%
  summarise(med = median(rho_AMP, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med)) %>%
  slice_head(n = 10)

top10 <- top10_tbl[[tax_col]]

# 2) subset data to top 10 and make clean labels (strip suffix after "_")
#vdf_top10 <- df %>%
#  filter(.data[[tax_col]] %in% top10) %>%
#  mutate(clean_label = str_replace(.data[[tax_col]], "_.*", ""))
vdf_top10 <- df %>%
  filter(.data[[tax_col]] %in% top10)

# 3) order by median again (ensures factor order matches top10 ranking)
order_top10 <- vdf_top10 %>%
  group_by(.data[[tax_col]]) %>%
  summarise(med = median(rho_AMP, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med)) %>%
  pull(.data[[tax_col]])


#order_clean <- str_replace(order_top10, "_.*", "")

# 4) violin plot: all violins dark blue, italic x labels, no y tick labels
#p_violin <- ggplot(vdf_top10,
#                   aes(x = factor(clean_label, levels = order_clean),
#                       y = rho_AMP)) +
p_violin <- ggplot(vdf_top10, aes(x=factor(.data[[taxon_level]], levels=order_top10), y=rho_AMP)) +
  geom_violin(fill = "darkblue", color = "grey20", scale = "width", trim = TRUE, alpha = 0.9) +
  stat_summary(fun = "median", geom = "point", size = 1.6, color = "white") +
  geom_jitter(width = 0.2, alpha = 0.35, size = 1, color = "black") +
  xlab(paste0(tax_col, " (Top 10 by median ρAMP; ≥", min_genomes, " genomes)")) +
  ylab("AMP density (ρAMP)") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    axis.text.y = element_blank(),     # hide y tick labels
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

print(p_violin)

## Box
p_box <- ggplot(vdf, aes(x=factor(.data[[taxon_level]], levels=order), y=rho_AMP)) +
  geom_boxplot(outlier.size=1.8, width=0.7) +
  geom_jitter(width=0.2, alpha=0.35, size=1) +
  labs(x=paste0(taxon_level," (≥",min_genomes," genomes)"), y="AMP density (ρAMP)",
       title=paste("Per-genome AMP density by", taxon_level,"— box")) +
  theme_bw(base_size=12) + theme(axis.text.x=element_text(angle=45,hjust=1))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
print(p_box)

ggsave(file.path(outdir,"figures",paste0("box_",tolower(taxon_level),"_density.png")), p_box, width=max(10,0.35*length(order)), height=6, dpi=300)
ggsave(file.path(outdir,"figures",paste0("box_",tolower(taxon_level),"_density.pdf")), p_box, width=max(10,0.35*length(order)), height=6)
