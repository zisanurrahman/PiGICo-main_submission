#!/usr/bin/env Rscript

# Hard-coded version: reads two TSVs (health + FE), runs k-means on UMAP and t-SNE,
# and saves plots with four hulls (group hulls filled; k-means hulls dashed).

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(rlang)
})

# =============================
# ---- USER-HARDCODED PART ----
# =============================
HEALTH_FILE <- "/Users/zisan/SGI_Paper/tSNE/health_coords_v2.tsv"   # must contain columns: UMAP1, UMAP2, TSNE1, TSNE2, Health
FE_FILE     <- "/Users/zisan/SGI_Paper/tSNE/fe_coords_v2.tsv"       # must contain columns: UMAP1, UMAP2, TSNE1, TSNE2, FE
OUTDIR      <- "/Users/zisan/SGI_Paper/tSNE/kmeans_plots"                   # output folder
K           <- 2                                 # k-means k
SEED        <- 42                                # random seed
POINT_ALPHA <- 0.7
POINT_SIZE  <- 1.3

# =============================
# ---------- helpers ----------
require_cols <- function(df, cols) {
  miss <- setdiff(cols, colnames(df))
  if (length(miss) > 0) {
    stop(sprintf("Missing required columns: %s", paste(miss, collapse = ", ")))
  }
}

hull_df <- function(df, x, y, group) {
  # Robust convex-hull builder without tidy-eval trickery
  gx <- rlang::as_name(ensym(x))
  gy <- rlang::as_name(ensym(y))
  gg <- rlang::as_name(ensym(group))
  
  # Keep only needed cols and complete cases
  cols_needed <- c(gx, gy, gg)
  d <- df[, cols_needed, drop = FALSE]
  d <- d[stats::complete.cases(d), , drop = FALSE]
  if (nrow(d) == 0) return(d[0, , drop = FALSE])
  
  # Split by group and compute hull per group
  pieces <- split(d, d[[gg]], drop = TRUE)
  res <- lapply(pieces, function(s) {
    if (nrow(s) < 3) return(s)  # not enough points for a hull; return as-is
    idx <- tryCatch(grDevices::chull(s[[gx]], s[[gy]]), error = function(e) integer())
    if (length(idx) == 0) return(s)
    s[idx, , drop = FALSE]
  })
  out <- do.call(rbind, res)
  rownames(out) <- NULL
  out
}

plot_embedding <- function(df, xcol, ycol, group_col, cluster_col, title, outpath,
                           point_alpha = 0.7, point_size = 1.3) {
  group_sym   <- rlang::sym(group_col)
  cluster_sym <- rlang::sym(cluster_col)
  x_sym <- rlang::sym(xcol)
  y_sym <- rlang::sym(ycol)
  
  # dynamic axis titles
  x_title <- if (grepl("^UMAP", xcol, ignore.case = TRUE)) "UMAP 1" else
    if (grepl("^TSNE", xcol, ignore.case = TRUE)) "t-SNE 1" else xcol
  y_title <- if (grepl("^UMAP", ycol, ignore.case = TRUE)) "UMAP 2" else
    if (grepl("^TSNE", ycol, ignore.case = TRUE)) "t-SNE 2" else ycol
  
  # legend titles
  group_title <- if (identical(group_col, "Health")) "Health" else group_col
  
  # hulls
  hull_group   <- hull_df(df, !!x_sym, !!y_sym, !!group_sym)
  hull_cluster <- hull_df(df, !!x_sym, !!y_sym, !!cluster_sym)
  
  p <- ggplot(df, aes(x = !!x_sym, y = !!y_sym)) +
    # points (Health legend)
    geom_point(aes(color = !!group_sym), alpha = point_alpha, size = point_size) +
    # filled group hulls (Health legend)
    #geom_polygon(
    #  data = hull_group,
    #  aes(fill = !!group_sym, group = !!group_sym),
    #  alpha = 0.10, color = NA
    #) +
    # dashed cluster outlines (Cluster legend)
    #geom_polygon(
    #  data = hull_cluster,
    #  aes(fill = !!cluster_sym, group = !!cluster_sym),
    #  alpha = 0.15, color = "black", linewidth = 0.5, linetype = "solid"
    #) +
    # hulls shaded by the same group color (linked to fill)
    geom_polygon(
      data = hull_group,
      aes(fill = !!group_sym, group = !!group_sym),
      alpha = 0.15, color = "black", linewidth = 0.5
    ) +
    # custom color scales
    scale_color_manual(
      name = group_title,
      values = c("HEALTHY" = "darkslateblue", "PWD" = "cornflowerblue",
                 "LFE" = "coral3", "HFE" = "darkolivegreen4")
    ) +
    scale_fill_manual(
      name = group_title,
      values = c("HEALTHY" = "darkslateblue", "PWD" = "cornflowerblue",
                 "LFE" = "coral3", "HFE" = "darkolivegreen4")
    ) +
    # separate legend headers
  # separate legend headers and ordering
   #scale_color_discrete(name = "Health") +
   #scale_fill_discrete(name  = "Cluster") +
   guides(
     color = guide_legend(order = 1),
     fill  = guide_legend(order = 2)
   ) +
    theme_classic(base_size = 12) +
    theme(
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.title.x = element_text(size = 7, face = "bold"),
      axis.title.y = element_text(size = 7, face = "bold"),
      legend.position = "right",
      legend.box = "vertical",
      legend.key.size = unit(1, "lines"),
      legend.title     = element_text(size = 8),  # legend title size/style
      legend.text      = element_text(size = 6),  # legend text size
      legend.spacing.x = unit(0.5, "cm")             # horizontal gap between the two boxes
    ) +
    theme(legend.position = 'none')+
    labs(
      x = x_title,
      y = y_title
    )
  
  ggsave(outpath, p, width = 2, height = 1.8, dpi = 600,units = "in")
  message("✓ Saved: ", outpath)
}

run_one <- function(input, group_col, outdir, base_prefix) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  message("[LOAD] ", input)
  df <- suppressMessages(readr::read_tsv(input, guess_max = 100000))
  num_cols <- intersect(c("UMAP1","UMAP2","TSNE1","TSNE2"), colnames(df))
  if (length(num_cols) > 0) df[num_cols] <- lapply(df[num_cols], as.numeric)
  
  require_cols(df, c("UMAP1", "UMAP2", "TSNE1", "TSNE2", group_col))
  
  # ---- UMAP kmeans ----
  set.seed(SEED)
  umap_mat <- as.matrix(df[, c("UMAP1", "UMAP2")])
  km_u <- kmeans(umap_mat, centers = K, nstart = 50)
  df$Cluster_UMAP <- paste0("C", km_u$cluster)
  
  # ---- t-SNE kmeans ----
  set.seed(SEED)
  tsne_mat <- as.matrix(df[, c("TSNE1", "TSNE2")])
  km_t <- kmeans(tsne_mat, centers = K, nstart = 50)
  df$Cluster_TSNE <- paste0("C", km_t$cluster)
  
  # ---- plots ----
  plot_embedding(
    df, "UMAP1", "UMAP2", group_col, "Cluster_UMAP",
    #title = sprintf("UMAP + k-means (k=%d) — %s", K, base_prefix),
    outpath = file.path(outdir, sprintf("%s_umap_kmeans_hulls_k%d.png", base_prefix, K)),
    point_alpha = POINT_ALPHA, point_size = POINT_SIZE
  )
  
  
  plot_embedding(
    df, "TSNE1", "TSNE2", group_col, "Cluster_TSNE",
    #title = sprintf("t-SNE + k-means (k=%d) — %s", K, base_prefix),
    outpath = file.path(outdir, sprintf("%s_tsne_kmeans_hulls_k%d.png", base_prefix, K)),
    point_alpha = POINT_ALPHA, point_size = POINT_SIZE
  )
  
  
  # also save clustered tables
  clustered_out <- file.path(outdir, sprintf("%s_with_clusters_k%d.tsv", base_prefix, K))
  readr::write_tsv(df, clustered_out)
  message("✓ Saved table: ", clustered_out)
}

# =============================
# Run for HEALTH (if present)
# =============================
if (file.exists(HEALTH_FILE)) {
  run_one(HEALTH_FILE, group_col = "Health", outdir = OUTDIR, base_prefix = "health_coords")
} else {
  message("[SKIP] HEALTH file not found: ", HEALTH_FILE)
}

# =============================
# Run for FE (if present)
# =============================
if (file.exists(FE_FILE)) {
  run_one(FE_FILE, group_col = "FE", outdir = OUTDIR, base_prefix = "fe_coords")
} else {
  message("[SKIP] FE file not found: ", FE_FILE)
}

