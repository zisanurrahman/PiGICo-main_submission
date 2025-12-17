#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(ggpubr)  # For stat_compare_means
})

# ======================= CONFIG =======================
per_genome_tsv <- "/Users/zisan/Library/CloudStorage/Sync/SGI_guilds/SGI_guilds/Scripts/Final_scripts/AMP_density/amp_density_per_genome.tsv"
taxonomy_tsv   <- "/Users/zisan/SGI_Paper/Full_paper_clean/Data/gtdbtk.bac120.summary.tsv"
outdir         <- "/Users/zisan/SGI_Paper/Full_paper_clean/Figures/amp_meta_plots"
taxon_level    <- "Genus"        # "Genus", "Phylum", "Family", "Species"
min_genomes    <- 3              # groups need ≥ this many genomes (for violin/box)
topN_bubble    <- 10             # top-N for bubble
dir.create(file.path(outdir, "figures"), showWarnings = FALSE, recursive = TRUE)

# ======================= HELPERS =======================
`%||%` <- function(x, y) if (is.null(x)) y else x
parse_rank <- function(class_str, prefix){
  m <- stringr::str_match(class_str %||% "", paste0(prefix, "([^;]+)"))
  ifelse(is.na(m[,2]), NA_character_, m[,2])
}

# ======================= READ =======================
dens <- read_tsv(per_genome_tsv, show_col_types = FALSE) %>%
  mutate(
    rho_AMP   = as.numeric(rho_AMP),
    n_cAMPs   = as.numeric(n_cAMPs),
    Length_bp = as.numeric(Length_bp)
  )

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

# ======================= BUBBLE PLOT =======================
bubble_data <- df %>%
  filter(!is.na(.data[[taxon_level]]), !is.na(rho_AMP)) %>%
  group_by(.data[[taxon_level]]) %>%
  summarise(
    mean_rho  = mean(rho_AMP, na.rm=TRUE),
    sd_rho    = sd(rho_AMP, na.rm=TRUE),
    n_genomes = dplyr::n(),
    prevalence= mean(replace_na(n_cAMPs,0) >= 1),
    .groups   = "drop"
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
  theme_bw(base_size=12) +
  theme(panel.border = element_rect(color="black", fill=NA, linewidth=0.4))
print(p_bubble)
ggsave(file.path(outdir,"figures",paste0("bubble_",tolower(taxon_level),"_density.png")), p_bubble, width=10, height=7, dpi=300, bg="white")
ggsave(file.path(outdir,"figures",paste0("bubble_",tolower(taxon_level),"_density.pdf")), p_bubble, width=10, height=7)

# ======================= VIOLIN (TOP 10 BY MEDIAN) WITH STATS =======================
# 1) top 10 groups by median rho_AMP (after n>=min_genomes filter)
top10_tbl <- df %>%
  filter(!is.na(.data[[taxon_level]]), !is.na(rho_AMP)) %>%
  group_by(.data[[taxon_level]]) %>%
  filter(dplyr::n() >= min_genomes) %>%
  summarise(med = median(rho_AMP, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med)) %>%
  slice_head(n = 10)

top10 <- top10_tbl[[taxon_level]]

# 2) subset data to top 10
vdf_top10 <- df %>%
  filter(.data[[taxon_level]] %in% top10, !is.na(rho_AMP))

# 3) order factor levels by median rho_AMP
order_top10 <- vdf_top10 %>%
  group_by(.data[[taxon_level]]) %>%
  summarise(med = median(rho_AMP, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med)) %>%
  pull(.data[[taxon_level]])

# 4) Perform Kruskal-Wallis test
# Create a formula using the taxon_level variable
formula_kw <- as.formula(paste("rho_AMP ~", taxon_level))
kw_test <- kruskal.test(formula_kw, data = vdf_top10)
kw_p <- kw_test$p.value

# 5) Prepare for pairwise comparisons
ref_group <- order_top10[1]
comparisons_vs_ref <- lapply(order_top10[-1], function(x) c(ref_group, x))

# After step 5, add this to calculate and display the p-values

# 5) Prepare for pairwise comparisons
ref_group <- order_top10[1]
comparisons_vs_ref <- lapply(order_top10[-1], function(x) c(ref_group, x))

# Calculate actual p-values for each comparison
pairwise_pvalues <- data.frame(
  comparison = character(),
  p_value = numeric(),
  p_adj = numeric(),
  significance = character(),
  stringsAsFactors = FALSE
)

for (comp in comparisons_vs_ref) {
  group1_data <- vdf_top10 %>% filter(.data[[taxon_level]] == comp[1]) %>% pull(rho_AMP)
  group2_data <- vdf_top10 %>% filter(.data[[taxon_level]] == comp[2]) %>% pull(rho_AMP)
  
  test_result <- wilcox.test(group1_data, group2_data)
  
  pairwise_pvalues <- rbind(pairwise_pvalues, data.frame(
    comparison = paste(comp[1], "vs", comp[2]),
    p_value = test_result$p.value,
    stringsAsFactors = FALSE
  ))
}

# Apply Benjamini-Hochberg correction
pairwise_pvalues$p_adj <- p.adjust(pairwise_pvalues$p_value, method = "BH")

# Add significance symbols
pairwise_pvalues$significance <- case_when(
  pairwise_pvalues$p_adj < 0.0001 ~ "****",
  pairwise_pvalues$p_adj < 0.001 ~ "***",
  pairwise_pvalues$p_adj < 0.01 ~ "**",
  pairwise_pvalues$p_adj < 0.05 ~ "*",
  TRUE ~ "ns"
)

# Print the table
print("Pairwise comparisons (all vs top genus):")
print(pairwise_pvalues)

# Save to file
write_tsv(pairwise_pvalues, 
          file.path(outdir, "figures", paste0("pairwise_wilcox_", tolower(taxon_level), ".tsv")))

# Also print in a nicer format
cat("\n=== P-value Legend ===\n")
cat("ns: p > 0.05\n")
cat("*: p ≤ 0.05\n")
cat("**: p ≤ 0.01\n")
cat("***: p ≤ 0.001\n")
cat("****: p < 0.0001\n\n")

cat("=== Comparisons to", ref_group, "===\n")
for (i in 1:nrow(pairwise_pvalues)) {
  cat(sprintf("%s: p = %.4e (adj. p = %.4e) %s\n", 
              pairwise_pvalues$comparison[i],
              pairwise_pvalues$p_value[i],
              pairwise_pvalues$p_adj[i],
              pairwise_pvalues$significance[i]))
}

# 6) violin plot with statistics
p_violin <- ggplot(vdf_top10,
                   aes(x=factor(.data[[taxon_level]], levels=order_top10), y=rho_AMP)) +
  geom_violin(fill = "#1B9E77", color = "grey20", scale = "width", trim = TRUE, alpha = 0.4) +
  stat_summary(fun = "median", geom = "point", size = 1.6, color = "darkred") +
  geom_jitter(width = 0.2, alpha = 0.8, size = 1, color = "black") +
  
  # Add pairwise comparisons (Wilcoxon test)
  stat_compare_means(
    aes(group = .data[[taxon_level]]),
    comparisons = comparisons_vs_ref,
    method = "wilcox.test",
    label = "p.signif",  # Use "p.format" for actual p-values
    hide.ns = TRUE,      # Hide non-significant comparisons
    size = 3
  ) +
  
  # Add overall Kruskal-Wallis p-value
  stat_compare_means(
    aes(group = .data[[taxon_level]]),
    method = "kruskal.test",
    label.y = max(vdf_top10$rho_AMP, na.rm = TRUE) * 1.15,
    size = 3.5
  ) +
  
  xlab(paste0(taxon_level, " (Top 10 ρAMP; ≥", min_genomes, " genomes)")) +
  ylab("AMP Density (ρAMP)") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )
print(p_violin)
ggsave(file.path(outdir,"figures",paste0("violin_",tolower(taxon_level),"_density_stats.tiff")),
       p_violin, width=5.5, height=6.5, dpi=600, compression='lzw', bg="white")
ggsave(file.path(outdir,"figures",paste0("violin_",tolower(taxon_level),"_density_stats.pdf")),
       p_violin, width=6.5, height=6, bg="white")

# Also save summary statistics table
stats_summary <- vdf_top10 %>%
  group_by(.data[[taxon_level]]) %>%
  summarise(
    n = n(),
    median = median(rho_AMP, na.rm = TRUE),
    mean = mean(rho_AMP, na.rm = TRUE),
    sd = sd(rho_AMP, na.rm = TRUE),
    q25 = quantile(rho_AMP, 0.25, na.rm = TRUE),
    q75 = quantile(rho_AMP, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median))

write_tsv(stats_summary, file.path(outdir, "figures", paste0("violin_", tolower(taxon_level), "_stats_summary.tsv")))

#Option 1: Compare Top Genus vs All Others
# Create binary comparison: top genus vs rest
top_genus <- order_top10[1]  # Highest median genus

vdf_top10_binary <- vdf_top10 %>%
  mutate(group_comparison = ifelse(.data[[taxon_level]] == top_genus, 
                                   top_genus, 
                                   "All other genera"))

# Statistical test
wilcox_result <- wilcox.test(rho_AMP ~ group_comparison, data = vdf_top10_binary)

# Plot
p_violin_binary <- ggplot(vdf_top10_binary,
                          aes(x = group_comparison, y = rho_AMP)) +
  geom_violin(aes(fill = group_comparison), color = "grey20", 
              scale = "width", trim = TRUE, alpha = 0.6) +
  scale_fill_manual(values = c("darkred", "grey70")) +
  stat_summary(fun = "median", geom = "point", size = 2, color = "black") +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
  stat_compare_means(method = "wilcox.test", 
                     label = "p.format",
                     size = 4.5,
                     label.y = max(vdf_top10_binary$rho_AMP) * 1.05) +
  labs(x = NULL, y = "AMP Density (ρAMP)",
       title = paste0(top_genus, " vs Other Top Genera")) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

print(p_violin_binary)
ggsave(file.path(outdir,"figures",paste0("violin_",tolower(taxon_level),"_binary_comparison.tiff")),
       p_violin_binary, width=5, height=5.5, dpi=500, compression='lzw', bg="white")


#Option 2: Show Top 3-5 with Pairwise Comparisons
# Focus on top 5 for cleaner visualization
top5 <- order_top10[1:5]

vdf_top5 <- vdf_top10 %>%
  filter(.data[[taxon_level]] %in% top5)

# Compare top genus to others
ref_group <- top5[1]
comparisons_vs_top <- lapply(top5[-1], function(x) c(ref_group, x))

p_violin_top5 <- ggplot(vdf_top5,
                        aes(x = factor(.data[[taxon_level]], levels = top5), 
                            y = rho_AMP)) +
  geom_violin(fill = "#1B9E77", color = "grey20", 
              scale = "width", trim = TRUE, alpha = 0.5) +
  stat_summary(fun = "median", geom = "point", size = 2, color = "darkred") +
  geom_jitter(width = 0.2, alpha = 0.7, size = 1.5, color = "black") +
  stat_compare_means(
    comparisons = comparisons_vs_top,
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = FALSE,  # Show all comparisons including ns
    size = 4,
    tip.length = 0.02,
    step.increase = 0.08
  ) +
  xlab(NULL) +
  ylab("AMP Density (ρAMP)") +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 11),
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

print(p_violin_top5)

#Option 3: Highlight Top Genus with Color

# Color code: top genus vs others
vdf_top10_colored <- vdf_top10 %>%
  mutate(highlight = ifelse(.data[[taxon_level]] == order_top10[1], 
                            "Highest", 
                            "Other"))

p_violin_highlight <- ggplot(vdf_top10_colored,
                             aes(x = factor(.data[[taxon_level]], levels = order_top10), 
                                 y = rho_AMP,
                                 fill = highlight)) +
  geom_violin(color = "grey20", scale = "width", trim = TRUE, alpha = 0.6) +
  scale_fill_manual(values = c("Highest" = "#D55E00", "Other" = "grey80")) +
  stat_summary(fun = "median", geom = "point", size = 1.8, color = "black") +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
  
  # Show Kruskal-Wallis
  stat_compare_means(
    method = "kruskal.test",
    label.y = max(vdf_top10$rho_AMP) * 1.1,
    size = 4
  ) +
  
  xlab(NULL) +
  ylab("AMP Density (ρAMP)") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

print(p_violin_highlight)


# 6) violin plot with statistics - only show significant comparisons
# Alternative approach with multcompView
# 6) violin plot with statistics - only show significant comparisons
p_violin <- ggplot(vdf_top10,
                   aes(x=factor(.data[[taxon_level]], levels=order_top10), y=rho_AMP)) +
  geom_violin(fill = "#1B9E77", color = "grey20", scale = "width", trim = TRUE, alpha = 0.4) +
  stat_summary(fun = "median", geom = "point", size = 1.6, color = "darkred") +
  geom_jitter(width = 0.2, alpha = 0.8, size = 1, color = "black") +
  
  # Add pairwise comparisons (Wilcoxon test) - only significant ones
  stat_compare_means(
    comparisons = comparisons_vs_ref,
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE,  # This hides non-significant
    size = 3.5,
    tip.length = 0.01,
    step.increase = 0.05,
    bracket.size = 0.4
  ) +
  
  # Add overall Kruskal-Wallis p-value
  annotate("text", 
           x = length(order_top10)/2 + 0.5, 
           y = max(vdf_top10$rho_AMP, na.rm = TRUE) * 1.18,
           label = paste0("Kruskal-Wallis, p = ", format.pval(kw_p, digits = 2)),
           size = 4) +
  
  xlab(paste0(taxon_level, " (Top 10 by median ρAMP; ≥", min_genomes, " genomes)")) +
  ylab("AMP Density (ρAMP)") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 11),
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )
print(p_violin)
ggsave(file.path(outdir,"figures",paste0("violin_",tolower(taxon_level),"_density_stats.tiff")),
       p_violin, width=7.5, height=6, dpi=500, compression='lzw', bg="white")

# ======================= BOX (ALL GROUPS WITH n>=min_genomes) =======================
vdf_all <- df %>%
  filter(!is.na(.data[[taxon_level]]), !is.na(rho_AMP)) %>%
  group_by(.data[[taxon_level]]) %>%
  filter(dplyr::n() >= min_genomes) %>%
  ungroup()

order_all <- vdf_all %>%
  group_by(.data[[taxon_level]]) %>%
  summarise(med = median(rho_AMP, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med)) %>%
  pull(.data[[taxon_level]])

p_box <- ggplot(vdf_all, aes(x=factor(.data[[taxon_level]], levels=order_all), y=rho_AMP)) +
  geom_boxplot(outlier.size=1.8, width=0.7, fill="grey90", color="black") +
  geom_jitter(width=0.2, alpha=0.35, size=1) +
  labs(x=paste0(taxon_level," (≥",min_genomes," genomes)"), y="AMP density (ρAMP)") +
  theme_bw(base_size=12) +
  theme(
    axis.text.x=element_text(angle=45,hjust=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
  )
w_box <- max(10, 0.35*length(order_all))
ggsave(file.path(outdir,"figures",paste0("box_",tolower(taxon_level),"_density.png")),
       p_box, width=w_box, height=6, dpi=300, bg="white")
ggsave(file.path(outdir,"figures",paste0("box_",tolower(taxon_level),"_density.pdf")),
       p_box, width=w_box, height=6, bg="white")