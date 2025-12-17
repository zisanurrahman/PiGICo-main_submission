#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
})

# =========================
# ====== INPUTS (edit) ====
# =========================
seqkit_tsv <- "/Users/zisan/SGI_Paper/Full_paper_clean/Data/genome_stats_out/assembly_stats_with_group.tsv"
gtdb_path  <- "/Users/zisan/SGI_Paper/New_ML_Analysis_Trans/gtdb_taxonomy.csv"

# EITHER provide a precomputed cluster map (from earlier FastANI pipeline)...
annot_clusters_csv <- "/Users/zisan/SGI_Paper/Full_paper_clean/Data/genome_stats_out/fastani_diversity/genome_cluster_map_species95_strain99.csv"
# ...OR provide the FastANI all-vs-all TSV and we will compute clusters here:
fastani_tsv <- "/Users/zisan/SGI_Paper/Full_paper_clean/Data/genome_stats_out/fastani_allpairs.tsv"

# Outputs
out_taxa    <- "/Users/zisan/SGI_Paper/Full_paper_clean/Data/genome_stats_out/diversity_taxa"
out_fastani <- "/Users/zisan/SGI_Paper/Full_paper_clean/Data/genome_stats_out/fastani_diversity"
dir.create(out_taxa,    showWarnings = FALSE, recursive = TRUE)
dir.create(out_fastani, showWarnings = FALSE, recursive = TRUE)

# Try ComplexUpset if installed (for UpSet plots)
have_upset <- requireNamespace("ComplexUpset", quietly = TRUE)

# Colors
col_map <- c("SGI"="#1B9E77", "PiBac"="darkslateblue", "Union"="#666666")

# =========================
# ====== Utilities ========
# =========================
strip_ext <- function(x) sub("\\.[^.]+$", "", basename(x))

# Common writer for a simple richness bar (SGI, PiBac, Union)
plot_richness_bar <- function(n_SGI, n_PiBac, n_union, title, outfile) {
  rich_df <- tibble(group=c("SGI","PiBac","Union"), n=c(n_SGI,n_PiBac,n_union))
  p <- ggplot(rich_df, aes(group, n, fill = group)) +
    geom_col(color = "black") +
    scale_fill_manual(values = col_map) +
    labs(x = NULL, y = "Unique count", title = title) +
    theme_classic(base_size = 12) +
    theme(panel.border = element_rect(color="black", fill=NA), legend.position = "none")
  ggsave(outfile, p, width = 4.8, height = 3.2, dpi = 300)
}

# =========================
# == A) TAXONOMIC (GTDB) ==
# =========================
# Load SGI/PiBac mapping
map <- read_tsv(seqkit_tsv, show_col_types = FALSE) %>%
  transmute(genome_file = file,
            genome_stem = strip_ext(file),
            Group = factor(Group, levels = c("SGI","PiBac","Other"))) %>%
  filter(Group %in% c("SGI","PiBac")) %>%
  select(genome_stem, Group) %>%
  distinct()

# Load GTDB taxonomy
gtdb <- read_csv(gtdb_path, show_col_types = FALSE) %>%
  mutate(genome = as.character(genome),
         genome_stem = strip_ext(genome))

tax <- gtdb %>%
  inner_join(map, by = "genome_stem") %>%
  select(genome, genome_stem, Group, genus, species) %>%
  mutate(
    genus   = na_if(trimws(genus),   ""),
    species = na_if(trimws(species), "")
  )

# Generic analyzer for a taxonomic level (genus/species)
analyze_tax_level <- function(df, level = c("genus","species"), label, outdir) {
  level <- match.arg(level)
  sets <- df %>%
    filter(!is.na(.data[[level]])) %>%
    distinct(.data[[level]], Group) %>%
    mutate(value = TRUE) %>%
    pivot_wider(names_from = Group, values_from = value, values_fill = FALSE) %>%
    mutate(SGI = ifelse(is.na(SGI), FALSE, SGI),
           PiBac = ifelse(is.na(PiBac), FALSE, PiBac))
  
  unique_SGI    <- sets %>% filter(SGI & !PiBac)
  unique_PiBac  <- sets %>% filter(!SGI & PiBac)
  shared        <- sets %>% filter(SGI & PiBac)
  
  n_SGI   <- nrow(sets %>% filter(SGI))
  n_PiBac <- nrow(sets %>% filter(PiBac))
  n_union <- nrow(sets %>% filter(SGI | PiBac))
  
  # CSVs
  write_csv(unique_SGI,   file.path(outdir, paste0("unique_", label, "_SGI.csv")))
  write_csv(unique_PiBac, file.path(outdir, paste0("unique_", label, "_PiBac.csv")))
  write_csv(shared,       file.path(outdir, paste0("shared_", label, ".csv")))
  
  # UpSet (if available). NOTE: ComplexUpset>=1.3 uses `intersect=`.
  if (have_upset) {
    p_up <- ComplexUpset::upset(
      sets %>% select(any_of(c(level, "SGI","PiBac"))),
      intersect = c("SGI","PiBac"),
      base_annotations = list(
        "Intersection size" = ComplexUpset::intersection_size(text = list(size = 3))
      ),
      width_ratio = 0.2
    ) + labs(title = paste0("UpSet: ", tools::toTitleCase(label), " overlap (SGI vs PiBac)"))
    ggsave(file.path(outdir, paste0("upset_", label, ".png")), p_up, width = 6, height = 4, dpi = 300)
  }
  
  # Richness bar
  plot_richness_bar(n_SGI, n_PiBac, n_union,
                    paste0(tools::toTitleCase(label), " richness"),
                    file.path(outdir, paste0("richness_", label, ".png")))
  
  list(
    text = paste0(
      "=== ", toupper(label), " ===\n",
      "Unique to SGI   : ", nrow(unique_SGI),   "\n",
      "Unique to PiBac : ", nrow(unique_PiBac), "\n",
      "Shared          : ", nrow(shared),       "\n",
      "Total SGI       : ", n_SGI, " | Total PiBac: ", n_PiBac, " | Union: ", n_union, "\n",
      "Lift (Union - SGI)   : +", n_union - n_SGI,   "\n",
      "Lift (Union - PiBac) : +", n_union - n_PiBac, "\n\n"
    )
  )
}

gen_res <- analyze_tax_level(tax, "genus",   "genera",  out_taxa)
sp_res  <- analyze_tax_level(tax, "species", "species", out_taxa)
writeLines(paste0(gen_res$text, sp_res$text),
           file.path(out_taxa, "uniques_shared_summary.txt"))
message("Taxonomy uniques/shared written to: ", out_taxa)

# =========================================
# == B) FASTANI species/strain CLUSTERS  ==
# =========================================
# If a precomputed cluster map exists, use it; otherwise compute from FastANI TSV.
if (file.exists(annot_clusters_csv)) {
  annot <- read_csv(annot_clusters_csv, show_col_types = FALSE) %>%
    mutate(Group = factor(Group, levels = c("SGI","PiBac")))
} else {
  # Build clusters from FastANI pairs
  fastani <- read_tsv(
    fastani_tsv,
    col_names = c("query","reference","ANI","frags_mapped","frags_total"),
    show_col_types = FALSE
  ) %>%
    mutate(q = basename(query), r = basename(reference)) %>%
    filter(q != r, is.finite(ANI))
  
  grp_of <- function(x) ifelse(str_detect(x, "^SGI"), "SGI", "PiBac")
  all_genomes <- sort(unique(c(fastani$q, fastani$r)))
  
  components_at <- function(edges_tbl, threshold, all_nodes) {
    ed <- edges_tbl %>%
      filter(ANI >= threshold) %>%
      transmute(from = q, to = r) %>%
      distinct()
    g <- graph_from_data_frame(ed, directed = FALSE, vertices = tibble(name = all_nodes))
    comp <- components(g)
    tibble(genome = names(comp$membership),
           cluster_id = as.integer(comp$membership))
  }
  
  species_map <- components_at(fastani, 95, all_genomes) %>%
    rename(species_cluster = cluster_id)
  strain_map  <- components_at(fastani, 99, all_genomes) %>%
    rename(strain_cluster  = cluster_id)
  
  annot <- tibble(genome = all_genomes,
                  Group = factor(grp_of(all_genomes), levels = c("SGI","PiBac"))) %>%
    left_join(species_map, by = "genome") %>%
    left_join(strain_map,  by = "genome")
  
  write_csv(annot, file.path(out_fastani, "genome_cluster_map_species95_strain99.csv"))
}

# Generic uniques/shared over a cluster column (species_cluster/strain_cluster)
analyze_cluster_level <- function(df, cl_col, label, outdir) {
  df2 <- df %>% select(genome, Group, cluster = {{cl_col}}) %>% distinct()
  
  sets <- df2 %>%
    filter(!is.na(cluster)) %>%
    distinct(cluster, Group) %>%
    mutate(value = TRUE) %>%
    pivot_wider(names_from = Group, values_from = value, values_fill = FALSE) %>%
    mutate(SGI = ifelse(is.na(SGI), FALSE, SGI),
           PiBac = ifelse(is.na(PiBac), FALSE, PiBac))
  
  unique_SGI    <- sets %>% filter(SGI & !PiBac)
  unique_PiBac  <- sets %>% filter(!SGI & PiBac)
  shared        <- sets %>% filter(SGI & PiBac)
  
  n_SGI   <- nrow(sets %>% filter(SGI))
  n_PiBac <- nrow(sets %>% filter(PiBac))
  n_union <- nrow(sets %>% filter(SGI | PiBac))
  
  # CSVs
  write_csv(unique_SGI,   file.path(outdir, paste0("unique_", label, "_clusters_SGI.csv")))
  write_csv(unique_PiBac, file.path(outdir, paste0("unique_", label, "_clusters_PiBac.csv")))
  write_csv(shared,       file.path(outdir, paste0("shared_", label, "_clusters.csv")))
  
  # UpSet (if available)
  if (have_upset) {
    p_up <- ComplexUpset::upset(
      sets %>% select(cluster, SGI, PiBac),
      intersect = c("SGI","PiBac"),
      base_annotations = list(
        "Intersection size" = ComplexUpset::intersection_size(text = list(size = 3))
      ),
      width_ratio = 0.2
    ) + labs(title = paste0("UpSet: ", tools::toTitleCase(label), " clusters (SGI vs PiBac)"))
    ggsave(file.path(outdir, paste0("upset_", label, "_clusters.png")), p_up, width = 6, height = 4, dpi = 300)
  }
  
  # Richness bar
  plot_richness_bar(n_SGI, n_PiBac, n_union,
                    paste0(tools::toTitleCase(label), " cluster richness"),
                    file.path(outdir, paste0("richness_", label, "_clusters.png")))
  
  paste0(
    "=== ", toupper(label), " CLUSTERS ===\n",
    "Unique to SGI   : ", nrow(unique_SGI),   "\n",
    "Unique to PiBac : ", nrow(unique_PiBac), "\n",
    "Shared          : ", nrow(shared),       "\n",
    "Total SGI       : ", n_SGI, " | Total PiBac: ", n_PiBac, " | Union: ", n_union, "\n",
    "Lift (Union - SGI)   : +", n_union - n_SGI,   "\n",
    "Lift (Union - PiBac) : +", n_union - n_PiBac, "\n\n"
  )
}

txt_species <- analyze_cluster_level(annot, species_cluster, "species", out_fastani)
txt_strain  <- analyze_cluster_level(annot, strain_cluster,  "strain",  out_fastani)
writeLines(paste0(txt_species, txt_strain),
           file.path(out_fastani, "unique_shared_cluster_summary.txt"))

message("Done.\n- Taxonomy outputs: ", out_taxa,
        "\n- Cluster outputs: ", out_fastani,
        "\n(UpSet plots require ComplexUpset; install.packages('ComplexUpset') if missing.)")
