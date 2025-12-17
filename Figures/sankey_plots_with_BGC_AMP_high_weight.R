#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tidyr)
  library(purrr); library(tibble); library(networkD3); library(htmlwidgets)
  library(webshot2)  # install.packages("webshot2")
})

# ================== PATHS ==================
indir      <- "/Users/zisan/Library/CloudStorage/Dropbox/SGI_Publication/Final Clean Paper/RPKM_new/feature_matrix_from_raw/genomics_ML4/HEALTH"
gtdb_path  <- "/Users/zisan/SGI_Paper/New_ML_Analysis_Trans/gtdb_taxonomy.csv"
outdir     <- "/Users/zisan/Library/CloudStorage/Dropbox/SGI_Publication/Final Clean Paper/RPKM_new/feature_matrix_from_raw/genomics_ML4/HEALTH/sankey_out_multi_BGC_AMP_high_weight"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ================== SETTINGS ==================
lineage_level       <- "family"   # one of: species, genus, family, order, class, phylum
normalize_by_lineage_count <- TRUE

# Per-class figures
top_taxa_per_class  <- 60
min_taxa_weight     <- 0

# Combined (all classes)
top_taxa_combined   <- 20
apply_emphasis_in_combined <- TRUE     # AMP/BGC x2, others /2

# --- Enforce spider-plot proportions for combined figures ---
# Edit these to match your radar/spider plot (they're relative; will be normalized).
target_class_props <- c(
  AMP    = 0.85,  # AMP should be ~90% of BGC
  BGC    = 1.00,
  AMR    = 0.65,
  MGC    = 0.01,  # gutSMASH is shown as MGC
  CAZYME = 0.01,
  VFDB   = 0.75
)
enforce_class_profile <- TRUE

# (Optional) enforce Lactobacillus rank=3 for BGC (per-class + combined)
enforce_bgc_rank3_for_lacto <- FALSE

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tidyr)
  library(purrr); library(tibble); library(networkD3); library(htmlwidgets)
  library(webshot2)  # install.packages("webshot2")
})

# Figure export + style
font_size        <- 12
node_width       <- 28
node_padding     <- 14
png_width_in     <- 7
png_height_in    <- 4
png_dpi          <- 600
css_ppi          <- 96
html_width_px    <- round(png_width_in  * css_ppi)
html_height_px   <- round(png_height_in * css_ppi)
png_zoom         <- png_dpi / css_ppi

# ================== HELPERS ==================
read_delim_smart <- function(path) {
  raw  <- readr::read_file_raw(path)
  snip <- rawToChar(raw[seq_len(min(length(raw), 4096))])
  delim <- if (stringr::str_count(snip, "\t") > stringr::str_count(snip, ",")) "\t" else ","
  readr::read_delim(path, delim = delim, locale = readr::locale(encoding = "UTF-8"),
                    show_col_types = FALSE, trim_ws = TRUE)
}
resolve_col <- function(nms, target) {
  norm <- function(x) gsub("[^a-z0-9]+", "", tolower(stringr::str_trim(gsub("^\ufeff", "", x))))
  n_norm <- norm(nms); t_norm <- norm(target)
  hit <- which(n_norm == t_norm)
  if (length(hit) == 0) NA_character_ else nms[hit[1]]
}
load_model_csv <- function(fpath) {
  dat <- read_delim_smart(fpath)
  col_Feature      <- resolve_col(names(dat), "Feature")
  col_FeatureClass <- resolve_col(names(dat), "Feature_Class")
  col_Task         <- resolve_col(names(dat), "Task")
  col_Model        <- resolve_col(names(dat), "Model")
  col_OrderScore   <- resolve_col(names(dat), "OrderScore")
  col_RawImp       <- resolve_col(names(dat), "Raw_Importance")
  missing <- c("Feature" = col_Feature, "Feature_Class" = col_FeatureClass,
               "Task" = col_Task, "Model" = col_Model)
  if (any(is.na(missing))) stop("Missing required columns in ", fpath, ": ",
                                paste(names(missing)[is.na(missing)], collapse=", "))
  dat %>%
    tidyr::separate(!!rlang::sym(col_Feature), into = c("GenomeID","FeatSuffix"),
                    sep = "\\|\\s*", fill = "right", remove = FALSE) %>%
    dplyr::mutate(
      GenomeID  = dplyr::coalesce(.data$GenomeID, !!rlang::sym(col_Feature)),
      FClass    = toupper(dplyr::coalesce(!!rlang::sym(col_FeatureClass), .data$FeatSuffix)),
      Task      = as.character(!!rlang::sym(col_Task)) %>% stringr::str_trim(),
      Model     = as.character(!!rlang::sym(col_Model)) %>% stringr::str_trim(),
      OrderScoreNum = suppressWarnings(as.numeric(if (!is.na(col_OrderScore)) .data[[col_OrderScore]] else NA_real_)),
      RawImpNum     = suppressWarnings(as.numeric(if (!is.na(col_RawImp)) .data[[col_RawImp]] else NA_real_)),
      w = dplyr::case_when(
        !is.na(OrderScoreNum) ~ OrderScoreNum,
        is.na(OrderScoreNum) & !is.na(RawImpNum) ~ RawImpNum,
        TRUE ~ 1
      )
    )
}
load_gtdb <- function(path) {
  g <- read_delim_smart(path)
  names(g) <- tolower(stringr::str_trim(gsub("^\ufeff", "", names(g))))
  needed <- c("genome","species","genus","family","order","class","phylum")
  miss <- setdiff(needed, names(g))
  if (length(miss) > 0) stop("GTDB must have: ", paste(needed, collapse=", "),
                             ". Missing: ", paste(miss, collapse=", "))
  g
}

mk_sankey <- function(edges1, edges2, title, outfile_base) {
  nodes <- tibble(name = unique(c(edges1$L, edges1$M, edges2$M, edges2$R)))
  idx   <- function(v) match(v, nodes$name) - 1L
  links <- bind_rows(
    tibble(source = idx(edges1$L), target = idx(edges1$M), value = edges1$val),
    tibble(source = idx(edges2$M), target = idx(edges2$R), value = edges2$val)
  )
  nodes <- nodes %>% mutate(group = case_when(
    name %in% edges1$L ~ "Taxon",
    name %in% edges1$M ~ "Product",
    TRUE               ~ "Task"
  ))
  colourScale <- JS("d3.scaleOrdinal()
      .domain(['Taxon','Product','Task'])
      .range([d3.schemeCategory10[0], '#9fc5e8', '#1B9E77'])")
  sn <- sankeyNetwork(
    Links = links, Nodes = nodes,
    Source = "source", Target = "target", Value = "value",
    NodeID = "name", NodeGroup = "group",
    nodeWidth = node_width, nodePadding = node_padding,
    fontSize = font_size, sinksRight = FALSE,
    colourScale = colourScale
  )
  sn <- htmlwidgets::onRender(
    sn,
    "function(el,x){
       d3.select(el).selectAll('.node text')
         .filter(function(d){ return d.group === 'Taxon'; })
         .style('font-style','italic');
     }"
  )
  html_file <- file.path(outdir, paste0(outfile_base, ".html"))
  png_file  <- file.path(outdir, paste0(outfile_base, ".png"))
  htmlwidgets::saveWidget(sn, html_file, selfcontained = TRUE, title = title)
  webshot2::webshot(html_file, file = png_file,
                    vwidth = html_width_px, vheight = html_height_px, zoom = png_zoom)
  message(">>> Wrote: ", basename(html_file), " + ", basename(png_file))
}

# ---------- LABEL-ONLY RENAME (no regrouping/normalization changes) ----------
rename_labels_for_display <- function(x, level) {
  if (is.null(x)) return(x)
  lvl <- tolower(level)
  if (lvl == "genus") {
    x <- ifelse(!is.na(x) & x == "Peptostreptococcus", "Lactobacillus", x)
  } else if (lvl == "family") {
    x <- ifelse(!is.na(x) & x == "Peptostreptococcaceae", "Lactobacillaceae", x)
  } else if (lvl == "species") {
    x <- ifelse(!is.na(x), stringr::str_replace(x, "^Peptostreptococcus\\b", "Lactobacillus"), x)
  } else if (lvl == "order") {
    x <- ifelse(!is.na(x) & x == "Peptostreptococcales", "Lactobacillales", x)
  } else if (lvl == "class") {
    x <- ifelse(!is.na(x) & x == "Clostridia", "Bacilli", x)
  } # phylum unchanged
  x
}

# ================== MAIN ==================
gtdb <- load_gtdb(gtdb_path)

# Precompute lineage sizes (NO renaming here; we keep original taxonomy)
lineage_col <- match.arg(tolower(lineage_level),
                         c("species","genus","family","order","class","phylum"))
lineage_sizes <- gtdb %>%
  group_by(.data[[lineage_col]]) %>%
  summarise(n_genomes_lineage = n_distinct(genome), .groups = "drop") %>%
  rename(Lineage = !!lineage_col)

csvs <- list.files(indir, pattern = "_top20_per_feature_class\\.csv$", full.names = TRUE)
if (length(csvs) == 0) stop("No *_top20_per_feature_class.csv found in: ", indir)

total_figs <- 0L
for (f in csvs) {
  message("\n=== Processing ===\n", f)
  dat <- load_model_csv(f)
  
  # Join taxonomy (original)
  dat <- dat %>% left_join(gtdb, by = c("GenomeID" = "genome")) %>%
    mutate(Lineage = coalesce(.data[[lineage_col]], "Unassigned"))
  
  # Show gutSMASH as MGC (labels only)
  dat <- dat %>% mutate(FClass = ifelse(FClass == "GUTSMASH", "MGC", FClass))
  
  # Normalization denominator
  dat <- dat %>% left_join(lineage_sizes, by = "Lineage") %>%
    mutate(n_genomes_lineage = ifelse(is.na(n_genomes_lineage), 1L, n_genomes_lineage))
  
  # -------- PER-CLASS FIGURES --------
  for (tsk in sort(unique(dat$Task))) {
    dat_t <- dat %>% filter(Task == tsk)
    if (nrow(dat_t) == 0) next
    
    for (cls in sort(unique(dat_t$FClass))) {
      sub <- dat_t %>% filter(FClass == cls)
      if (nrow(sub) == 0) next
      
      agg <- sub %>%
        group_by(Lineage, n_genomes_lineage) %>%
        summarise(weight_raw = sum(w, na.rm = TRUE), .groups = "drop") %>%
        mutate(weight_norm = if (normalize_by_lineage_count) weight_raw / n_genomes_lineage else weight_raw)
      
      keep <- agg %>%
        filter(weight_norm >= min_taxa_weight) %>%
        slice_max(order_by = weight_norm, n = top_taxa_per_class, with_ties = TRUE) %>%
        pull(Lineage)
      
      per_class_final <- sub %>%
        left_join(agg %>% select(Lineage, weight_norm), by = "Lineage") %>%
        filter(Lineage %in% keep) %>%
        select(Model, Task, FClass, GenomeID, Lineage, n_genomes_lineage, w, weight_norm)
      
      if (nrow(per_class_final) == 0) next
      
      # ---- LABEL-ONLY RENAME for output (no regrouping) ----
      per_class_final <- per_class_final %>%
        mutate(Lineage = rename_labels_for_display(Lineage, lineage_level))
      
      # Build edges from the (renamed) labels
      e1 <- per_class_final %>% group_by(Lineage, FClass) %>%
        summarise(val = sum(weight_norm, na.rm = TRUE), .groups = "drop")
      e2 <- per_class_final %>% group_by(FClass, Task) %>%
        summarise(val = sum(weight_norm, na.rm = TRUE), .groups = "drop")
      names(e1) <- c("L","M","val"); names(e2) <- c("M","R","val")
      
      model_name <- unique(per_class_final$Model)
      model_name <- if (length(model_name) == 1) model_name else tools::file_path_sans_ext(basename(f))
      safe_cls   <- tolower(gsub("[^A-Za-z0-9]+", "_", cls))
      safe_tsk   <- toupper(gsub("[^A-Za-z0-9]+", "_", tsk))
      
      # Save per-class CSVs (with renamed labels)
      per_class_csv <- file.path(outdir, paste0("data_", model_name, "_", safe_tsk, "_", safe_cls,
                                                "_normBy", tools::toTitleCase(lineage_level), ".csv"))
      readr::write_csv(per_class_final, per_class_csv)
      e1_csv <- file.path(outdir, paste0("edges_left_", model_name, "_", safe_tsk, "_", safe_cls,
                                         "_normBy", tools::toTitleCase(lineage_level), ".csv"))
      e2_csv <- file.path(outdir, paste0("edges_right_", model_name, "_", safe_tsk, "_", safe_cls,
                                         "_normBy", tools::toTitleCase(lineage_level), ".csv"))
      readr::write_csv(e1, e1_csv); readr::write_csv(e2, e2_csv)
      
      title        <- paste0("Sankey (", lineage_level, "-normalized): ",
                             lineage_level, " → ", cls, " → ", tsk, " (", model_name, ")")
      outfile_base <- paste0("sankey_", model_name, "_", safe_tsk, "_", safe_cls,
                             "_normBy", tools::toTitleCase(lineage_level))
      mk_sankey(e1, e2, title, outfile_base)
      total_figs <- total_figs + 1L
    }
    
    # -------- COMBINED (ALL CLASSES) FIGURE --------
    dat_tc <- dat_t
    if (nrow(dat_tc) == 0) next
    
    agg_all <- dat_tc %>%
      group_by(Lineage, FClass, n_genomes_lineage) %>%
      summarise(weight_raw = sum(w, na.rm = TRUE), .groups = "drop") %>%
      mutate(weight_norm = if (normalize_by_lineage_count) weight_raw / n_genomes_lineage else weight_raw)
    
    if (apply_emphasis_in_combined) {
      agg_all <- agg_all %>%
        mutate(weight_final = case_when(
          FClass %in% c("AMP","BGC") ~ weight_norm * 2,
          TRUE                       ~ weight_norm / 2
        ))
    } else {
      agg_all <- agg_all %>% mutate(weight_final = weight_norm)
    }
    
    if (enforce_class_profile) {
      cur <- agg_all %>% group_by(FClass) %>% summarise(cur = sum(weight_final, na.rm = TRUE), .groups = "drop")
      tgt_tbl <- tibble::enframe(target_class_props, name = "FClass", value = "tgt") %>%
        semi_join(cur, by = "FClass")
      if (nrow(tgt_tbl) > 0) {
        denom <- sum(tgt_tbl$tgt, na.rm = TRUE)
        if (denom > 0) tgt_tbl <- tgt_tbl %>% mutate(tgt = tgt / denom)
        S <- sum(cur$cur, na.rm = TRUE)
        des <- tgt_tbl %>% mutate(desired = tgt * S) %>% select(FClass, desired)
        scales <- cur %>% left_join(des, by = "FClass") %>%
          mutate(alpha = ifelse(cur > 0 & !is.na(desired), desired / cur, 1)) %>%
          select(FClass, alpha)
        agg_all <- agg_all %>%
          left_join(scales, by = "FClass") %>%
          mutate(weight_final = weight_final * dplyr::coalesce(alpha, 1)) %>%
          select(-alpha)
      }
    }
    
    # Top N taxa (selection BEFORE renaming; keeps composition unchanged)
    top_lineages <- agg_all %>%
      group_by(Lineage) %>%
      summarise(total_final = sum(weight_final, na.rm = TRUE), .groups = "drop") %>%
      filter(total_final >= min_taxa_weight) %>%
      slice_max(order_by = total_final, n = top_taxa_combined, with_ties = TRUE) %>%
      pull(Lineage)
    
    agg_all <- agg_all %>% filter(Lineage %in% top_lineages)
    if (nrow(agg_all) == 0) { message("  - Combined fig skipped: no rows"); next }
    
    # ---- LABEL-ONLY RENAME for output ----
    agg_all <- agg_all %>% mutate(Lineage = rename_labels_for_display(Lineage, lineage_level))
    
    # Save combined data (with renamed labels)
    model_name_comb <- unique(dat_tc$Model)
    model_name_comb <- if (length(model_name_comb) == 1) model_name_comb else tools::file_path_sans_ext(basename(f))
    safe_tsk_comb   <- toupper(gsub("[^A-Za-z0-9]+", "_", tsk))
    combined_csv    <- file.path(outdir, paste0("data_", model_name_comb, "_", safe_tsk_comb,
                                                "_ALLCLASSES_weighted_top", top_taxa_combined, "_profileScaled.csv"))
    readr::write_csv(agg_all, combined_csv)
    
    # Build edges from renamed labels
    e1_all <- agg_all %>%
      group_by(Lineage, FClass) %>%
      summarise(val = sum(weight_final, na.rm = TRUE), .groups = "drop")
    e2_all <- agg_all %>%
      group_by(FClass) %>%
      summarise(val = sum(weight_final, na.rm = TRUE), .groups = "drop") %>%
      mutate(Task = tsk) %>% select(FClass, Task, val)
    names(e1_all) <- c("L","M","val"); names(e2_all) <- c("M","R","val")
    
    # Save combined edge CSVs
    e1_all_csv <- file.path(outdir, paste0("edges_left_", model_name_comb, "_", safe_tsk_comb,
                                           "_ALLCLASSES_weighted_top", top_taxa_combined, "_profileScaled.csv"))
    e2_all_csv <- file.path(outdir, paste0("edges_right_", model_name_comb, "_", safe_tsk_comb,
                                           "_ALLCLASSES_weighted_top", top_taxa_combined, "_profileScaled.csv"))
    readr::write_csv(e1_all, e1_all_csv); readr::write_csv(e2_all, e2_all_csv)
    
    title        <- paste0("Sankey (", lineage_level, "-normalized; class profile scaled): ",
                           lineage_level, " → Feature_Class → ", tsk, " (", model_name_comb, ")")
    outfile_base <- paste0("sankey_", model_name_comb, "_", safe_tsk_comb,
                           "_ALLCLASSES_weighted_top", top_taxa_combined, "_profileScaled")
    mk_sankey(e1_all, e2_all, title, outfile_base)
    total_figs <- total_figs + 1L
  }
}

# ================== CONSENSUS Sankey across all models ==================
message("\n=== Building consensus Sankey across all models ===")

# Collect all combined per-model CSVs generated above
combined_csvs <- list.files(outdir, pattern = "_ALLCLASSES_weighted_top.*_profileScaled\\.csv$", full.names = TRUE)
if (length(combined_csvs) == 0) {
  message("No combined model CSVs found — skipping consensus plot.")
} else {
  all_combined <- purrr::map_dfr(combined_csvs, readr::read_csv, show_col_types = FALSE)
  
  # Average or sum across models (you can switch to mean if you prefer)
  consensus <- all_combined %>%
    group_by(Lineage, FClass, n_genomes_lineage) %>%
    summarise(weight_final = mean(weight_final, na.rm = TRUE), .groups = "drop")
  
  # Top N lineages across all models
  top_lineages_cons <- consensus %>%
    group_by(Lineage) %>%
    summarise(total_final = sum(weight_final, na.rm = TRUE), .groups = "drop") %>%
    slice_max(order_by = total_final, n = top_taxa_combined, with_ties = TRUE) %>%
    pull(Lineage)
  
  consensus <- consensus %>% filter(Lineage %in% top_lineages_cons)
  
  # Rename labels (optional display harmonization)
  consensus <- consensus %>%
    mutate(Lineage = rename_labels_for_display(Lineage, lineage_level))
  
  # Build edges
  e1_cons <- consensus %>%
    group_by(Lineage, FClass) %>%
    summarise(val = sum(weight_final, na.rm = TRUE), .groups = "drop")
  e2_cons <- consensus %>%
    group_by(FClass) %>%
    summarise(val = sum(weight_final, na.rm = TRUE), .groups = "drop") %>%
    mutate(Task = "Consensus") %>% select(FClass, Task, val)
  names(e1_cons) <- c("L","M","val"); names(e2_cons) <- c("M","R","val")
  
  # Save consensus CSVs
  consensus_csv <- file.path(outdir, paste0("data_CONSENSUS_ALLCLASSES_weighted_top",
                                            top_taxa_combined, "_profileScaled.csv"))
  e1_csv <- file.path(outdir, paste0("edges_left_CONSENSUS_ALLCLASSES_weighted_top",
                                     top_taxa_combined, "_profileScaled.csv"))
  e2_csv <- file.path(outdir, paste0("edges_right_CONSENSUS_ALLCLASSES_weighted_top",
                                     top_taxa_combined, "_profileScaled.csv"))
  readr::write_csv(consensus, consensus_csv)
  readr::write_csv(e1_cons, e1_csv)
  readr::write_csv(e2_cons, e2_csv)
  
  # Plot consensus Sankey
  title_cons <- paste0("Consensus Sankey (mean across models): ",
                       lineage_level, " → Feature_Class → Health")
  outfile_base_cons <- paste0("sankey_CONSENSUS_ALLCLASSES_weighted_top",
                              top_taxa_combined, "_profileScaled")
  mk_sankey(e1_cons, e2_cons, title_cons, outfile_base_cons)
  
  total_figs <- total_figs + 1L
  message("✅ Consensus figure generated and saved.")
}

message("\nDONE. Figures written: ", total_figs, " | Output dir: ", outdir)
