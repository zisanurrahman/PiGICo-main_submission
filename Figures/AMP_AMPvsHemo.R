#!/usr/bin/env Rscript
#install.packages('ggExtra')
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
  library(ggExtra)
})

# ------------------ inputs ------------------
AMP_FILE <- "/Users/zisan/Library/CloudStorage/Sync/SGI_guilds/SGI_guilds/AMP_physchem/AMP_physchem_features.tsv"   # <-- save your snippet as TSV with those 7 columns
OUTDIR   <- "/Users/zisan/SGI_Paper/Full_paper_clean/Figures/amp_meta_plots"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ------------------ read --------------------
df <- read_tsv(AMP_FILE, show_col_types = FALSE) %>%
  mutate(
    RowID = paste0(Access, "_", row_number()),    # unique label per row
    AMP_probability = as.numeric(AMP_probability),
    Hemolytic_probability = as.numeric(Hemolytic_probability),
    Hemolytic = factor(Hemolytic, levels = c("NonHemo","Hemo")),
    seq_len = nchar(Sequence)
  )



head(df)

# A small palette for hemolysis state
hemo_cols <- c(NonHemo = "#1B9E77",Hemo = "darkslateblue")

# ------------------ 1) AMP probability (bar/lollipop) ------------------
p1 <- df %>%
  arrange(desc(AMP_probability)) %>%
  mutate(Access = factor(Access, levels = RowID)) %>%
  ggplot(aes(x = AMP_probability, y = RowID)) +
  geom_segment(aes(x = 0, xend = AMP_probability, y = RowID, yend = RowID), color = "grey70") +
  geom_point(size = 3) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
  labs(x = "AMP probability", y = NULL, title = "AMP probability per peptide") +
  theme_classic(base_size = 12) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4))
print(p1)

ggsave(file.path(OUTDIR, "01_amp_probability_per_access.png"), p1, width = 7, height = 3.8, dpi = 300, bg = "white")

# ------------------ 2) Hemolytic probability (colored by label) --------
p2 <- df %>%
  arrange(desc(Hemolytic_probability)) %>%
  mutate(Access = factor(RowID, levels = RowID)) %>%
  ggplot(aes(x = Hemolytic_probability, y = RowID, color = Hemolytic)) +
  geom_segment(aes(x = 0, xend = Hemolytic_probability, y = RowID, yend = RowID), color = "grey70") +
  geom_point(size = 3) +
  scale_color_manual(values = hemo_cols, guide = guide_legend(title = "Hemolytic")) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
  labs(x = "Hemolytic probability", y = NULL, title = "Hemolysis probability per peptide") +
  theme_classic(base_size = 12) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4))
print(p2)
ggsave(file.path(OUTDIR, "02_hemolytic_probability_per_access.png"), p2, width = 7, height = 3.8, dpi = 300, bg = "white")

# ------------------ 3) AMP vs Hemolytic scatter ------------------------
head(df)
p3 <- ggplot(df, aes(x = AMP_probability, y = Hemolytic_probability,
                     color = Hemolytic, shape = AMP_family, label = RowID)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey60") +
  geom_point(size = 2, alpha = 0.6) +
  #ggrepel::geom_text_repel(size = 3, max.overlaps = 20, seed = 1, show.legend = FALSE) +
  scale_color_manual(values = hemo_cols) +
  labs(x = "Bioactivity Probability", y = "Hemolytic Probability",
       color = "Hemolytic", shape = "AMP Family")+
  theme_classic(base_size = 12) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4))+
  theme(legend.position = c(0.8,0.22),
        legend.direction = 'vertical',
        legend.key.size = unit(0.5, "lines"),      # shrink legend keys
        legend.title = element_text(size = 7),     # smaller title
        legend.text  = element_text(size = 6),     # smaller text
        legend.spacing.x = unit(0.2, "cm"),        # tighten horizontal gap
        legend.spacing.y = unit(0.2, "cm"),        # tighten vertical gap
        legend.box.margin = margin(0, 0, 0, 0))     # remove external margin)

print(p3)

# add marginal density plots
p3_density <- ggMarginal(
  p3,
  type = "density",         # density plots instead of histograms
  groupColour = TRUE,       # densities colored by group (Hemolytic)
  groupFill = TRUE,
  alpha = 0.4
)
print(p3_density)

ggsave(file.path(OUTDIR, "03_amp_vs_hemolytic_scatter.png"), p3_density, width = 4, height = 4, dpi = 600, bg = "white")

# ------------------ 3.2) AMP vs Hemolytic scatter with colored by collection------------------------

df <- df %>%
  mutate(Project = case_when(
    str_detect(Genome, "SGI")   ~ "SGI",
    str_detect(Genome, "PiBac") ~ "PiBac",
    TRUE                        ~ "Other"
  ))

proj_cols <- c(SGI = "#1B9E77", PiBac = "darkslateblue", Other = "grey50")

p3_genome <- ggplot(df, aes(x = AMP_probability, y = Hemolytic_probability,
                     color = Project, shape = AMP_family, label = RowID)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey60") +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = proj_cols) +
  labs(x = "Bioactivity Probability", y = "Hemolytic Probability",
       color = "PorciBiome", shape = "AMP Family") +
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    legend.position = c(0.8, 0.22),
    legend.direction = "vertical",
    legend.key.size = unit(0.5, "lines"),
    legend.title = element_text(size = 7),
    legend.text  = element_text(size = 6),
    legend.spacing.x = unit(0.2, "cm"),
    legend.spacing.y = unit(0.2, "cm"),
    legend.box.margin = margin(0, 0, 0, 0)
  )

print(p3_genome)


# add marginal density plots
p3_genome_density <- ggMarginal(
  p3_genome,
  type = "density",         # density plots instead of histograms
  groupColour = TRUE,       # densities colored by group (Hemolytic)
  groupFill = TRUE,
  alpha = 0.4
)
print(p3_genome_density)
ggsave(file.path(OUTDIR, "03_amp_vs_hemolytic_scatter.png"), p3_genome_density, width = 4, height = 4, dpi = 600, bg = "white")

library(ggpubr)

# Compare AMP probability
wilcox_amp <- wilcox.test(AMP_probability ~ Project, data = df)
wilcox_amp$p.value

# Compare Hemolytic probability
wilcox_hemo <- wilcox.test(Hemolytic_probability ~ Project, data = df)
wilcox_hemo$p.value

p3_genome_stats <- ggplot(df, aes(x = Project, y = AMP_probability, color = Project)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
  scale_color_manual(values = proj_cols) +
  stat_compare_means(method = "wilcox.test", label = "p.signif") + 
  labs(x = NULL, y = "AMP probability") +
  theme_classic(base_size = 12)

print(p3_genome_stats)

# ------------------ 4) AMP family counts -------------------------------
p4 <- df %>%
  count(AMP_family) %>%
  ggplot(aes(x = AMP_family, y = n)) +
  geom_col(fill = "grey70", color = "black", width = 0.75) +
  geom_text(aes(label = n), vjust = -0.3, size = 3.5) +
  labs(x = "AMP family", y = "Count", title = "AMP family composition") +
  theme_classic(base_size = 12) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4))
print(p4)
ggsave(file.path(OUTDIR, "04_amp_family_counts.png"), p4, width = 4.8, height = 3.6, dpi = 300, bg = "white")

# ------------------ 5) Sequence length per AMP -------------------------
p5 <- df %>%
  arrange(desc(seq_len)) %>%
  mutate(RowID = factor(RowID, levels = RowID)) %>%
  ggplot(aes(x = seq_len, y = Access)) +
  geom_col(fill = "grey70", color = "black", width = 0.7) +
  geom_text(aes(label = seq_len), hjust = -0.1, size = 3.2) +
  labs(x = "Length (aa)", y = NULL, title = "Sequence length") +
  xlim(0, max(df$seq_len) * 1.15) +
  theme_classic(base_size = 12) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4))
print(p5)

ggsave(file.path(OUTDIR, "05_sequence_length_per_access.png"), p5, width = 6.2, height = 3.8, dpi = 300, bg = "white")

# ------------------ 6) Mini heatmap of probabilities -------------------
prob_long <- df %>%
  select(RowID, AMP_probability, Hemolytic_probability) %>%
  pivot_longer(-RowID, names_to = "Metric", values_to = "Value")

p6 <- ggplot(prob_long, aes(x = Metric, y = RowID, fill = Value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "darkred") +
  labs(x = NULL, y = NULL, fill = "Prob.") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
print(p6)
ggsave(file.path(OUTDIR, "06_probability_heatmap.png"), p6, width = 4.6, height = 3.8, dpi = 300, bg = "white")

# ------------------ (Optional) sequence logos by family ----------------
# install.packages("ggseqlogo") if needed
if (requireNamespace("ggseqlogo", quietly = TRUE)) {
  library(ggseqlogo)
  # one panel per family
  logo_out <- file.path(OUTDIR, "07_seqlogo_by_family.png")
  p7 <- df %>%
    group_by(AMP_family) %>%
    summarise(seqs = list(Sequence), .groups = "drop") %>%
    mutate(plot = purrr::map(seqs, ~ ggseqlogo::ggseqlogo(.x) + ggtitle("")))
  
  # simple grid of logos
  g <- patchwork::wrap_plots(p7$plot, ncol = 2) +
    patchwork::plot_annotation(title = "Sequence logos by AMP family")
  ggsave(logo_out, g, width = 8, height = 6, dpi = 300, bg = "white")
}
