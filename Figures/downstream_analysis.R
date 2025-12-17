#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(circlize)
  library(grid)
})

# ==========================
# 0) Paths (EDIT IF NEEDED)
# ==========================
matrix_long_file <- "/Users/zisan/SGI_Paper/Gutsmash_analysis/gutsmash_pathway_matrix.genome_normalized_R.tsv"  # long: Genome, Pathway, Presence[, Category]
tax_file         <- "/Users/zisan/SGI_Paper/New_ML_Analysis_Trans/gtdb_taxonomy.csv"                            # genome, domain, phylum, class, order, family, genus, species
outdir           <- "/Users/zisan/SGI_Paper/Gutsmash_analysis/Downstream_Analyses"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# optional packages
have_vegan <- requireNamespace("vegan", quietly = TRUE)
have_ape   <- requireNamespace("ape",   quietly = TRUE)
if (!have_vegan) message("Note: 'vegan' not found -> using fallback Jaccard.")
if (!have_ape)   message("Note: 'ape'   not found -> PCoA will be skipped.")

# ==========================
# 1) Helpers
# ==========================
normalize_genome_id <- function(x) {
  x0 <- sub("\\..*$", "", x)  # strip any .gbk/.fa/.region*.gbk etc
  parts <- strsplit(x0, "_", fixed = TRUE)
  vapply(parts, function(p) if (length(p) >= 2) paste(p[1], p[2], sep = "_") else x0, FUN.VALUE = "")
}

canonize <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

# binary Jaccard (fallback if vegan missing)
jaccard_binary <- function(M) {
  M <- as.matrix(M > 0) * 1
  n <- nrow(M)
  D <- matrix(0, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      a <- sum(M[i,] & M[j,])
      b <- sum(M[i,] & !M[j,])
      c <- sum(!M[i,] & M[j,])
      D[i,j] <- 1 - (if ((a + b + c) == 0) 0 else a / (a + b + c))
    }
  }
  as.dist(D)
}

diversity_indices <- function(v) {
  p <- v / sum(v)
  p <- p[p > 0]
  shannon <- if (length(p)) -sum(p * log(p)) else 0
  simpson <- if (length(p)) 1 - sum(p^2)     else 0
  richness <- sum(v > 0)
  tibble(Shannon = shannon, Simpson = simpson, Richness = richness)
}

accumulation_curve <- function(bin_mat, n_perm = 100, seed = 1) {
  set.seed(seed)
  n <- nrow(bin_mat)
  acc <- matrix(0, nrow = n, ncol = n_perm)
  for (p in seq_len(n_perm)) {
    ord <- sample(seq_len(n))
    cum <- rep(0, ncol(bin_mat))
    uniq <- integer(n)
    for (k in seq_len(n)) {
      cum <- cum | (bin_mat[ord[k], ] > 0)
      uniq[k] <- sum(cum)
    }
    acc[, p] <- uniq
  }
  tibble(
    genomes = seq_len(n),
    mean_unique = rowMeans(acc),
    sd_unique   = apply(acc, 1, sd)
  )
}

# Generic Fisher enrichment (binary presence) for a matrix by a taxonomic level
fisher_enrichment_matrix <- function(mat_binary, lookup_df, level = c("Order","Family","Genus"), feature_label = "Feature", min_n = 10) {
  level <- match.arg(level)
  df <- as_tibble(mat_binary, rownames = "BaseGenome") |>
    left_join(lookup_df, by = "BaseGenome") |>
    pivot_longer(-c(BaseGenome, Order, Family, Genus, Source),
                 names_to = feature_label, values_to = "Presence")
  df <- df |>
    mutate(Group = .data[[level]]) |>
    filter(!is.na(Group))
  
  groups <- df |> distinct(Group) |> pull(Group)
  
  out <- purrr::map_dfr(groups, function(g) {
    sub <- df |>
      mutate(in_group = Group == g)
    
    n_in <- sub |> filter(in_group) |> distinct(BaseGenome) |> nrow()
    if (n_in < min_n) return(tibble(Group = g, !!feature_label := character(0),
                                    n_in = n_in, p = numeric(0), odds = numeric(0)))
    
    sub |>
      group_by(.data[[feature_label]]) |>
      summarise(
        a = sum( Presence > 0 &  in_group),
        b = sum( Presence == 0 &  in_group),
        c = sum( Presence > 0 & !in_group),
        d = sum( Presence == 0 & !in_group),
        .groups = "drop"
      ) |>
      rowwise() |>
      mutate(
        fisher = list(fisher.test(matrix(c(a,b,c,d), nrow = 2))),
        p = fisher$p.value,
        odds = fisher$estimate %||% NA_real_
      ) |>
      ungroup() |>
      mutate(Group = g, n_in = n_in) |>
      select(Group, all_of(feature_label), n_in, a:b, c:d, p, odds)
  }) |>
    mutate(p_adj = p.adjust(p, method = "fdr"),
           log2_odds = log2(odds))
  
  out
}

volcano_plot <- function(df, title, outprefix, pad_mult = 0.04) {
  library(ggplot2)
  df2 <- df |>
    dplyr::filter(is.finite(log2_odds), is.finite(p_adj), p_adj > 0) |>
    dplyr::mutate(sig = p_adj < 0.05,
                  neglog10 = -log10(p_adj))
  
  p <- ggplot(df2, aes(x = log2_odds, y = neglog10, color = sig)) +
    geom_point(alpha = 1, size = 1.2, shape = 16) +
    scale_color_manual(values = c(`TRUE` = "#1B9E77", `FALSE` = "grey65"),
                       guide = guide_legend(title = NULL, override.aes = list(size = 3))) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey50") +
    # Add padding so point circles are fully visible near the edges
    scale_x_continuous(expand = expansion(mult = pad_mult)) +
    scale_y_continuous(expand = expansion(mult = pad_mult)) +
    # Prevent clipping at the panel border
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 11) +
    theme(panel.border = element_rect(color='black', linewidth = 0.5, fill = NA),
          legend.position = c(0.2,0.85))+
    theme(plot.margin = grid::unit(c(6, 10, 6, 10), "pt")) +
    labs(title = title, x = "log2(odds ratio)", y = "-log10(FDR)")
  
  ggsave(file.path(outdir, paste0(outprefix, ".png")),
         p, width = 3, height = 3, dpi = 600, bg = "white")
  ggsave(file.path(outdir, paste0(outprefix, ".tiff")),
         p, width = 3, height = 3, dpi = 600, bg = "white", compression = "lzw")
}



# ==========================
# 2) Load inputs
# ==========================
long <- read.delim(matrix_long_file, sep = "\t", header = TRUE, check.names = FALSE) |>
  as_tibble()

stopifnot(all(c("Genome","Pathway","Presence") %in% names(long)))
long <- long |>
  mutate(
    Genome   = as.character(Genome),
    Pathway  = as.character(Pathway),
    Presence = as.numeric(Presence)
  )
long$Presence[is.na(long$Presence)] <- 0
long$Presence <- as.numeric(long$Presence > 0)

# Ensure Category column exists
if (!"Category" %in% names(long)) long$Category <- "Other"
long <- long |>
  mutate(Category = if_else(is.na(Category) | Category == "", "Other", Category))

# Map taxonomy (Order/Family/Genus) to base genome IDs
tax <- readr::read_csv(tax_file, show_col_types = FALSE) |>
  transmute(
    genome,
    Order  = if_else(is.na(order)  | order  == "", "Unknown", order),
    Family = if_else(is.na(family) | family == "", "Unknown", family),
    Genus  = if_else(is.na(genus)  | genus  == "", "Unknown", genus)
  )

long <- long |>
  mutate(BaseGenome = normalize_genome_id(Genome)) |>
  left_join(tax, by = c("BaseGenome" = "genome"))

long$Order[is.na(long$Order)]   <- "Unknown"
long$Family[is.na(long$Family)] <- "Unknown"
long$Genus[is.na(long$Genus)]   <- "Unknown"

# Source (SGI / PiBac / Other) from *base* genome name
long <- long |>
  mutate(Source = case_when(
    grepl("SGI",   BaseGenome) ~ "SGI",
    grepl("PiBac", BaseGenome) ~ "PiBac",
    TRUE ~ "Other"
  ))

# ==========================
# 3) Wide matrices
# ==========================
# (a) Pathway presence per genome (binary)
wide_path <- long |>
  group_by(BaseGenome, Pathway) |>
  summarise(Presence = as.numeric(any(Presence > 0)), .groups = "drop") |>
  pivot_wider(names_from = Pathway, values_from = Presence, values_fill = 0) |>
  arrange(BaseGenome)

mat_path <- as.matrix(wide_path[, -1, drop = FALSE])
rownames(mat_path) <- wide_path$BaseGenome

# (b) Category presence per genome (binary: any pathway in category)
wide_cat <- long |>
  group_by(BaseGenome, Category) |>
  summarise(Presence = as.numeric(any(Presence > 0)), .groups = "drop") |>
  pivot_wider(names_from = Category, values_from = Presence, values_fill = 0) |>
  arrange(BaseGenome)

mat_cat <- as.matrix(wide_cat[, -1, drop = FALSE])
rownames(mat_cat) <- wide_cat$BaseGenome

# lookup
lookup <- long |>
  distinct(BaseGenome, Order, Family, Genus, Source)

# ==========================
# 4) Prevalence (Pathway/Category × Order/Family/Genus)
# ==========================
prevalence_by <- function(mat, level = c("Order","Family","Genus"), value_label) {
  level <- match.arg(level)
  as_tibble(mat, rownames = "BaseGenome") |>
    left_join(lookup, by = "BaseGenome") |>
    pivot_longer(-c(BaseGenome, Order, Family, Genus, Source),
                 names_to = value_label, values_to = "Presence") |>
    group_by(.data[[level]], .data[[value_label]]) |>
    summarise(
      n_genomes = n(),
      prevalence = mean(Presence > 0),
      .groups = "drop"
    ) |>
    rename(Level = 1)
}

prev_path_order  <- prevalence_by(mat_path, "Order",  "Pathway")
prev_path_family <- prevalence_by(mat_path, "Family", "Pathway")
prev_path_genus  <- prevalence_by(mat_path, "Genus",  "Pathway")

prev_cat_order   <- prevalence_by(mat_cat, "Order",  "Category")
prev_cat_family  <- prevalence_by(mat_cat, "Family", "Category")
prev_cat_genus   <- prevalence_by(mat_cat, "Genus",  "Category")

write_tsv(prev_path_order,  file.path(outdir, "prevalence_pathway_by_order.tsv"))
write_tsv(prev_path_family, file.path(outdir, "prevalence_pathway_by_family.tsv"))
write_tsv(prev_path_genus,  file.path(outdir, "prevalence_pathway_by_genus.tsv"))
write_tsv(prev_cat_order,   file.path(outdir, "prevalence_category_by_order.tsv"))
write_tsv(prev_cat_family,  file.path(outdir, "prevalence_category_by_family.tsv"))
write_tsv(prev_cat_genus,   file.path(outdir, "prevalence_category_by_genus.tsv"))

# Publication-style heatmaps (Category × Order / Family top-K / Genus top-K)
plot_prev_heat <- function(prev_df, row_label, title, out_base, topK = NA_integer_) {
  df <- prev_df
  if (!is.na(topK)) {
    top_groups <- df |>
      group_by(Level) |>
      summarise(n = max(n_genomes), .groups = "drop") |>
      arrange(desc(n)) |>
      slice_head(n = topK) |>
      pull(Level)
    df <- df |> filter(Level %in% top_groups)
  }
  p <- df |>
    ggplot(aes(!!sym(names(df)[2]), reorder(Level, prevalence, FUN = median), fill = prevalence)) +
    geom_tile(color = "white", size = 0.2) +
    scale_fill_gradient(low = "white", high = "#1B9E77", limits = c(0,1)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(face = 'italic')
    ) +
    labs(x = NULL, y = row_label, fill = "Prevalence")
  print(p)
  ggsave(file.path(outdir, paste0(out_base, ".png")),  p, width = 3.9, height = 7.5, dpi = 300, bg = "white")
  ggsave(file.path(outdir, paste0(out_base, ".tiff")), p, width = 3.9, height = 7.5, dpi = 300, bg = "white", compression = "lzw")
}

plot_prev_heat(prev_cat_order,  "Order",  "Prevalence of categories by Order",                    "prevalence_category_by_order_heat")
plot_prev_heat(prev_cat_family, "Family", "Prevalence of categories by Family (top 50 families)", "prevalence_category_by_family_topK_heat", topK = 50)
plot_prev_heat(prev_cat_genus,  "Genus",  "Prevalence of categories by Genus (top 50 genera)",    "prevalence_category_by_genus_topK_heat", topK = 50)

# ==========================
# 5) Enrichment (Pathway) by Order/Family/Genus
# ==========================
enrich_path_order  <- fisher_enrichment_matrix(mat_path, lookup, "Order",  feature_label = "Pathway", min_n = 10)
enrich_path_family <- fisher_enrichment_matrix(mat_path, lookup, "Family", feature_label = "Pathway", min_n = 10)
enrich_path_genus  <- fisher_enrichment_matrix(mat_path, lookup, "Genus",  feature_label = "Pathway", min_n = 10)

write_tsv(enrich_path_order,  file.path(outdir, "enrichment_pathway_by_order.tsv"))
write_tsv(enrich_path_family, file.path(outdir, "enrichment_pathway_by_family.tsv"))
write_tsv(enrich_path_genus,  file.path(outdir, "enrichment_pathway_by_genus.tsv"))

volcano_plot(enrich_path_order,  "Pathway enrichment by Order",  "volcano_enrich_pathway_order")
volcano_plot(enrich_path_family, "Pathway enrichment by Family", "volcano_enrich_pathway_family")
volcano_plot(enrich_path_genus,  "Pathway enrichment by Genus",  "volcano_enrich_pathway_genus")

# ==========================
# 6) Diversity (Shannon/Simpson/Richness)
# ==========================
# Per-genome: counts of pathways per category
genome_cat_counts <- long |>
  group_by(BaseGenome, Category) |>
  summarise(n_pathways = n_distinct(Pathway[Presence > 0]), .groups = "drop") |>
  pivot_wider(names_from = Category, values_from = n_pathways, values_fill = 0)

div_genome <- genome_cat_counts |>
  rowwise() |>
  mutate(tmp = list(diversity_indices(c_across(where(is.numeric))))) |>
  unnest(tmp) |>
  left_join(lookup, by = c("BaseGenome")) |>
  relocate(BaseGenome, Order, Family, Genus, Source)

write_tsv(div_genome, file.path(outdir, "diversity_indices_per_genome.tsv"))

# Per Order/Family/Genus (sum category counts across genomes)
agg_diversity <- function(df_long, level = c("Order","Family","Genus")) {
  level <- match.arg(level)
  df_long |>
    group_by(.data[[level]], BaseGenome, Category) |>
    summarise(n_pathways = n_distinct(Pathway[Presence > 0]), .groups = "drop") |>
    group_by(.data[[level]], Category) |>
    summarise(n_pathways = sum(n_pathways), .groups = "drop") |>
    pivot_wider(names_from = Category, values_from = n_pathways, values_fill = 0) |>
    rowwise() |>
    mutate(tmp = list(diversity_indices(c_across(where(is.numeric))))) |>
    unnest(tmp) |>
    rename(Level = 1)
}

div_order  <- agg_diversity(long, "Order")
div_family <- agg_diversity(long, "Family")
div_genus  <- agg_diversity(long, "Genus")

write_tsv(div_order,  file.path(outdir, "diversity_indices_by_order.tsv"))
write_tsv(div_family, file.path(outdir, "diversity_indices_by_family.tsv"))
write_tsv(div_genus,  file.path(outdir, "diversity_indices_by_genus.tsv"))

# quick plots
p1 <- ggplot(div_order, aes(x = reorder(Level, Shannon), y = Shannon)) +
  geom_col() + coord_flip() + theme_classic() +
  labs(title = "Shannon diversity of categories by Order", x = NULL, y = "Shannon")
ggsave(file.path(outdir, "diversity_shannon_by_order.png"), p1, width = 7, height = 8, dpi = 300, bg = "white")

p2 <- ggplot(div_family, aes(x = Shannon)) +
  geom_histogram(bins = 40) + theme_classic() +
  labs(title = "Shannon diversity of categories across Families", x = "Shannon", y = "Count")
ggsave(file.path(outdir, "diversity_shannon_family_hist.png"), p2, width = 6, height = 4, dpi = 300, bg = "white")

# ==========================
# 7) PCoA on Jaccard (Genome × Pathway)
# ==========================
if (nrow(mat_path) > 2 && ncol(mat_path) > 1 && have_ape) {
  D <- if (have_vegan) vegan::vegdist(mat_path, method = "jaccard", binary = TRUE) else jaccard_binary(mat_path)
  pcoa <- ape::pcoa(as.matrix(D))
  ord  <- as_tibble(pcoa$vectors[, 1:2], .name_repair = "minimal") |>
    mutate(BaseGenome = rownames(mat_path)) |>
    rename(PCoA1 = 1, PCoA2 = 2) |>
    left_join(lookup, by = "BaseGenome")
  
  ggplot(ord, aes(PCoA1, PCoA2, color = Order)) +
    geom_point(alpha = 0.7, size = 1.5) + theme_classic() +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(title = "PCoA (Jaccard) on pathway presence — colored by Order") +
    ggsave(file.path(outdir, "pcoa_jaccard_by_order.png"), width = 6, height = 5, dpi = 300, bg = "white")
  
  ggplot(ord, aes(PCoA1, PCoA2, color = Family)) +
    geom_point(alpha = 0.7, size = 1.5) + theme_classic() +
    labs(title = "PCoA (Jaccard) — colored by Family") +
    ggsave(file.path(outdir, "pcoa_jaccard_by_family.png"), width = 6, height = 5, dpi = 300, bg = "white")
  
  ggplot(ord, aes(PCoA1, PCoA2, color = Genus)) +
    geom_point(alpha = 0.7, size = 1.5) + theme_classic() +
    labs(title = "PCoA (Jaccard) — colored by Genus") +
    ggsave(file.path(outdir, "pcoa_jaccard_by_genus.png"), width = 6, height = 5, dpi = 300, bg = "white")
  
  ggplot(ord, aes(PCoA1, PCoA2, color = Source)) +
    geom_point(alpha = 0.8, size = 1.8) + theme_classic() +
    scale_color_manual(values = c(SGI = "#E41A1C", PiBac = "#377EB8", Other = "grey60")) +
    labs(title = "PCoA (Jaccard) — colored by Source (SGI / PiBac / Other)") +
    ggsave(file.path(outdir, "pcoa_jaccard_by_source.png"), width = 6, height = 5, dpi = 300, bg = "white")
} else {
  message("PCoA skipped (need >2 genomes & >1 pathway, and package 'ape').")
}

# ==========================
# 8) SGI vs PiBac comparisons
# ==========================
library(ggpubr)
library(rstatix)
#install.packages('effectsize')
library(effectsize)
sgi_pibac_prev <- as_tibble(mat_cat, rownames = "BaseGenome") |>
  left_join(lookup, by = "BaseGenome") |>
  filter(Source %in% c("SGI","PiBac")) |>
  pivot_longer(-c(BaseGenome, Order, Family, Genus, Source), names_to = "Category", values_to = "Presence") |>
  group_by(Source, Category) |>
  summarise(n = n(), prev = mean(Presence > 0), .groups = "drop") |>
  pivot_wider(names_from = Source, values_from = c(n, prev))
write_tsv(sgi_pibac_prev, file.path(outdir, "sgi_vs_pibac_prevalence_by_category.tsv"))

per_genome_catcount <- as_tibble(mat_cat, rownames = "BaseGenome") |>
  mutate(n_cat_present = rowSums(across(-BaseGenome) > 0)) |>
  left_join(lookup, by = "BaseGenome")
wil_out <- per_genome_catcount |>
  filter(Source %in% c("SGI","PiBac")) |>
  summarise(wilcox_p = wilcox.test(n_cat_present ~ Source)$p.value)
write_tsv(per_genome_catcount, file.path(outdir, "per_genome_category_counts.tsv"))
write_tsv(wil_out,            file.path(outdir, "sgi_vs_pibac_wilcox.tsv"))

df_box <- per_genome_catcount |> filter(Source %in% c("SGI","PiBac"))

# Wilcoxon + Cliff’s delta
wil_tbl <- df_box |>
  wilcox_test(n_cat_present ~ Source) |>
  adjust_pvalue(method = "BH") |>
  add_significance()
eff_tbl <- effectsize::cliffs_delta(
  x = df_box$n_cat_present,
  y = df_box$Source
)
print(eff_tbl)

# For plotting annotation with ggpubr
comp <- list(c("SGI","PiBac"))
library(scales)   # for pvalue()

df_box <- per_genome_catcount |> dplyr::filter(Source %in% c("SGI","PiBac"))
wil    <- wilcox.test(n_cat_present ~ Source, data = df_box, exact = FALSE)
p_lbl  <- paste0("Wilcoxon, p = ", scales::pvalue(wil$p.value))

y_top <- max(df_box$n_cat_present, na.rm = TRUE)

p_breadth <- ggplot(df_box, aes(Source, n_cat_present, fill = Source)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(values = c(SGI = "#1B9E77", PiBac = "#377EB8")) +
  # make room above; keep label off the panel clip
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  coord_cartesian(clip = "off") +
  annotate("text", x = 1.5, y = y_top * 1.3, label = p_lbl, size = 3.6) +
  theme_classic() +
  theme(panel.border = element_rect(color='black', linewidth = 0.5, fill = NA),
        legend.position = 'none',
        legend.background = element_rect(fill = NA)) +
  labs(y = "# categories present per genome", x = NULL)
print(p_breadth)
ggsave(file.path(outdir, "sgi_vs_pibac_category_breadth.png"),
       p_breadth, width = 3, height = 3, dpi = 600, bg = "white")


# Save stats tables too
readr::write_tsv(wil_tbl, file.path(outdir, "sgi_vs_pibac_wilcox_rstatix.tsv"))
readr::write_tsv(eff_tbl, file.path(outdir, "sgi_vs_pibac_cliffs_delta.tsv"))

ggsave(file.path(outdir, "sgi_vs_pibac_category_breadth.png"),
       p_breadth, width = 3, height = 3, dpi = 600, bg = "white")

eff_tbl <- effectsize::cliffs_delta(
  df_box$n_cat_present ~ df_box$Source
)
readr::write_tsv(as.data.frame(eff_tbl),
                 file.path(outdir, "sgi_vs_pibac_cliffs_delta.tsv"))# ==========================
# 9) Category co-occurrence (Jaccard) across genomes
# ==========================
if (ncol(mat_cat) >= 2) {
  # Distance between categories (columns), using binary Jaccard
  D_cat <- if (have_vegan) vegan::vegdist(t(mat_cat), method = "jaccard", binary = TRUE) else jaccard_binary(t(mat_cat))
  S_cat <- 1 - as.matrix(D_cat)
  dimnames(S_cat) <- list(colnames(mat_cat), colnames(mat_cat))
  S_cat[lower.tri(S_cat, diag = TRUE)] <- NA
  
  cooc_edges <- as.data.frame(S_cat, stringsAsFactors = FALSE) |>
    tibble::rownames_to_column("Category1") |>
    pivot_longer(cols = -Category1, names_to = "Category2", values_to = "Jaccard") |>
    filter(!is.na(Jaccard))
  
  readr::write_tsv(cooc_edges, file.path(outdir, "category_cooccurrence_jaccard.tsv"))
  
  thr <- 0.25
  readr::write_tsv(
    dplyr::filter(cooc_edges, Jaccard >= thr),
    file.path(outdir, sprintf("category_cooccurrence_edges_jaccard_ge%.2f.tsv", thr))
  )
  
  p_cooc_rank <- cooc_edges |>
    arrange(desc(Jaccard)) |>
    mutate(rank = row_number()) |>
    ggplot(aes(x = rank, y = Jaccard)) +
    geom_point(size = 0.8, alpha = 0.6) +
    geom_hline(yintercept = thr, linetype = 2, color = "grey50") +
    theme_classic(base_size = 11) +
    labs(title = "Category co-occurrence (Jaccard similarity)",
         x = "Pair rank (descending similarity)", y = "Jaccard")
  ggsave(file.path(outdir, "cooccurrence_category_ranked.png"),
         p_cooc_rank, width = 5, height = 3.5, dpi = 300, bg = "white")
  ggsave(file.path(outdir, "cooccurrence_category_ranked.tiff"),
         p_cooc_rank, width = 5, height = 3.5, dpi = 300, bg = "white", compression = "lzw")
  
  topN <- 50
  top_pairs <- cooc_edges |>
    arrange(desc(Jaccard)) |>
    slice_head(n = topN)
  
  p_cooc_heat <- top_pairs |>
    ggplot(aes(Category1, Category2, fill = Jaccard)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "#1B9E77") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = sprintf("Top %d category co-occurrences", topN), x = NULL, y = NULL)
  ggsave(file.path(outdir, "cooccurrence_category_topN_heat.png"),
         p_cooc_heat, width = 6, height = 5, dpi = 300, bg = "white")
  ggsave(file.path(outdir, "cooccurrence_category_topN_heat.tiff"),
         p_cooc_heat, width = 6, height = 5, dpi = 300, bg = "white", compression = "lzw")
}

# ==========================
# 10) Accumulation (category discovery) per Order/Family/Genus
# ==========================
make_accum_for_level <- function(level = c("Order","Family","Genus"), n_perm = 200) {
  level <- match.arg(level)
  df <- as_tibble(mat_cat, rownames = "BaseGenome") |>
    left_join(lookup, by = "BaseGenome")
  groups <- df |> distinct(.data[[level]]) |> pull(1)
  
  out <- map_dfr(groups, function(g) {
    sub <- df |> filter(.data[[level]] == g)
    if (nrow(sub) < 3) return(tibble(Level = level, Group = g, genomes = integer(0), mean_unique = numeric(0), sd_unique = numeric(0)))
    bin <- as.matrix(sub[, colnames(mat_cat), drop = FALSE])
    acc <- accumulation_curve(bin, n_perm = n_perm, seed = 42)
    acc |> mutate(Level = level, Group = g, .before = 1)
  })
  out
}

acc_order  <- make_accum_for_level("Order",  n_perm = 200)
acc_family <- make_accum_for_level("Family", n_perm = 200)
acc_genus  <- make_accum_for_level("Genus",  n_perm = 200)

write_tsv(acc_order,  file.path(outdir, "accumulation_categories_by_order.tsv"))
write_tsv(acc_family, file.path(outdir, "accumulation_categories_by_family.tsv"))
write_tsv(acc_genus,  file.path(outdir, "accumulation_categories_by_genus.tsv"))

if (nrow(acc_order) > 0) {
  p_acc <- ggplot(acc_order, aes(genomes, mean_unique, color = Group)) +
    geom_line() +
    theme_classic() +
    labs(x = "# genomes sampled", y = "Mean categories discovered",
         title = "Accumulation of categories with sampling — by Order")
  ggsave(file.path(outdir, "accumulation_categories_order.png"),
         p_acc, width = 7, height = 5, dpi = 300, bg = "white")
}

cat("✓ Done. Outputs written to: ", outdir, "\n")

# ---- Accumulation curves: categories discovered vs # genomes sampled (by Order) ----
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(readr); library(purrr)
})

# Inputs expected:
# mat_cat : matrix/data.frame with rows = BaseGenome (rownames), cols = categories (counts)
# lookup  : data.frame with columns BaseGenome, Order (and possibly others)
# outdir  : directory to save outputs

accum_by_group <- function(mat_cat, lookup, group_col = "Order",
                           reps = 200, seed = 42, min_n = 3) {
  
  stopifnot("BaseGenome" %in% c(colnames(lookup)))
  stopifnot(!is.null(rownames(mat_cat)))
  
  set.seed(seed)
  
  # make presence/absence
  pa <- as.matrix(mat_cat > 0) * 1L
  genomes <- rownames(pa)
  
  # attach group
  meta <- lookup %>%
    dplyr::select(BaseGenome, !!rlang::sym(group_col)) %>%
    dplyr::filter(BaseGenome %in% genomes)
  
  # split genomes by group
  groups <- split(meta$BaseGenome, meta[[group_col]])
  
  # keep groups with enough genomes
  groups <- groups[vapply(groups, length, 1L) >= min_n]
  
  # helper: mean categories discovered for each k in one group
  one_group_curve <- function(gvec) {
    g_pa <- pa[gvec, , drop = FALSE]
    n_g  <- nrow(g_pa)
    ks   <- 1:n_g
    
    res <- map_dfr(ks, function(k) {
      # sample k genomes, union categories, count; repeat reps
      counts <- replicate(reps, {
        sel <- sample.int(n_g, k, replace = FALSE)
        sum(colSums(g_pa[sel, , drop = FALSE]) > 0)
      })
      tibble(k = k,
             mean = mean(counts),
             sd   = sd(counts),
             se   = sd(counts)/sqrt(reps))
    })
    res
  }
  
  # compute for all groups
  df <- imap_dfr(groups, function(gvec, gname) {
    one_group_curve(gvec) %>% mutate(Group = gname, .before = 1)
  })
  
  df
}

# ==== Run ====
acc_df <- accum_by_group(mat_cat, lookup, group_col = "Order",
                         reps = 300, seed = 1, min_n = 3)

# Save the table (useful for supplement)
write_tsv(acc_df, file.path(outdir, "accumulation_by_order.tsv"))

# ==== Plot ====
p_acc <- ggplot(acc_df, aes(x = k, y = mean, color = Group)) +
  geom_line(linewidth = 0.7, alpha = 0.9) +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se), alpha = 0.12, color = NA) +
  labs(title = "Accumulation of categories with sampling — by Order",
       x = "# genomes sampled",
       y = "Mean categories discovered") +
  theme_classic(base_size = 12) +
  geom_line(linewidth = 0.7) +
  guides(fill = "none")+
  theme(legend.position = "right",
        legend.title = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6))
print(p_acc)
keep <- lookup %>% count(Order, name="n") %>% arrange(desc(n)) %>% slice_head(n=12) %>% pull(Order)
p_acc <- p_acc + scale_color_discrete(limits = keep, drop = TRUE)
print(p_acc)
ggsave(file.path(outdir, "accumulation_by_order.png"),
       p_acc, width = 9, height = 6, dpi = 300, bg = "white")
ggsave(file.path(outdir, "accumulation_by_order.tiff"),
       p_acc, width = 9, height = 6, dpi = 300, bg = "white", compression = "lzw")
