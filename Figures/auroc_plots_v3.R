# =================== CONFIG ===================
root_dir <- "/Users/zisan/Library/CloudStorage/Dropbox/SGI_Publication/Final Clean Paper/RPKM_new/feature_matrix_from_raw/genomics_ML4"
tasks    <- c("HEALTH", "FE")
models_order <- c("LR", "RF", "XGB")

# Consistent colors per model
model_cols <- c(LR = "#1F77B4", RF = "darkslateblue", XGB = "#1B9E77")

# Save sizes
single_w <- 3; single_h <- 3; single_dpi <- 600
combo_w  <- 4.5; combo_h <- 4.5; combo_dpi <- 600
# ==============================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(purrr); library(ggplot2)
})

theme_roc <- function(base_size = 11){
  theme_bw(base_size = base_size) +
    theme(panel.grid = element_blank(),
          legend.background = element_rect(color = "grey85", fill = "white"))
}

for (task in tasks) {
  task_dir <- file.path(root_dir, task)
  if (!dir.exists(task_dir)) {
    message("Skipping (missing dir): ", task_dir)
    next
  }
  
  # ----- Summary (for AUROC text & legend labels) -----
  sum_path <- file.path(task_dir, "auroc_summary.csv")
  if (!file.exists(sum_path)) {
    message("No auroc_summary.csv in ", task_dir, " — did the Python step run?")
    next
  }
  sum_df <- read_csv(sum_path, show_col_types = FALSE) %>%
    mutate(Model = toupper(Model)) %>%
    filter(Model %in% models_order)
  
  # To build legend labels like "LR AUROC = 0.739 ± 0.191"
  legend_labels <- sum_df %>%
    mutate(label = sprintf("%s AUROC = %.3f \u00B1 %.3f", Model, Mean_AUROC, SD_AUROC)) %>%
    select(Model, label)
  
  # For combined overlay data
  overlay_df <- NULL
  
  # ================= Per-model single plots =================
  for (m in models_order) {
    roc_path <- file.path(task_dir, sprintf("%s_mean_roc.csv", m))
    if (!file.exists(roc_path)) next
    
    df <- read_csv(roc_path, show_col_types = FALSE)
    val <- sum_df %>% filter(Model == m)
    if (nrow(val) == 0) next
    
    # In-plot text label (mean only)
    auc_label <- sprintf("AUROC = %.3f", val$Mean_AUROC)
    
    p_single <-
      ggplot(df, aes(FPR, TPR_mean)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.4) +
      geom_ribbon(
        aes(ymin = pmax(TPR_mean - TPR_sd, 0), ymax = pmin(TPR_mean + TPR_sd, 1)),
        fill = scales::alpha(model_cols[m], 0.18), color = NA, inherit.aes = TRUE
      ) +
      geom_line(linewidth = 0.9, color = model_cols[m]) +
      geom_text(
        aes(x = 0.50, y = 0.05), label = auc_label,
        inherit.aes = FALSE, size = 3.8, hjust = 0, color = "black"
      ) +
      labs(x = "False Positive Rate", y = "True Positive Rate") +
      theme_roc() +
      theme(legend.position = "none") # ← no legend on single plots
    
    ggsave(file.path(task_dir, sprintf("roc_%s.png", m)),  p_single,
           width = single_w, height = single_h, dpi = single_dpi)
    ggsave(file.path(task_dir, sprintf("roc_%s.tiff", m)), p_single,
           width = single_w, height = single_h, dpi = single_dpi, compression = "lzw")
    
    # add model tag for overlay
    df$Model <- m
    overlay_df <- bind_rows(overlay_df, df)
  }
  
  # ================= Combined overlay (all models) =================
  if (!is.null(overlay_df) && nrow(overlay_df)) {
    overlay_df <- overlay_df %>%
      mutate(Model = toupper(Model)) %>%
      left_join(legend_labels, by = "Model")
    
    # Build named palettes keyed by "label" so legend shows desired strings
    label_levels <- unique(overlay_df$label)
    label_to_color <- setNames(unname(model_cols[unique(overlay_df$Model)]), label_levels)
    label_to_fill  <- setNames(scales::alpha(unname(model_cols[unique(overlay_df$Model)]), 0.18), label_levels)
    
    p_overlay <-
      ggplot(overlay_df, aes(FPR, TPR_mean, color = label, fill = label)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.4,
                  inherit.aes = FALSE) +
      geom_ribbon(aes(ymin = pmax(TPR_mean - TPR_sd, 0),
                      ymax = pmin(TPR_mean + TPR_sd, 1)),
                  alpha = 0.14, color = NA) +
      geom_line(linewidth = 0.9) +
      scale_color_manual(values = label_to_color) +
      scale_fill_manual(values  = label_to_fill) +
      labs(x = "False Positive Rate", y = "True Positive Rate")+
           #title = sprintf("%s – Mean ROC (±1 SD)", task), color = NULL, fill = NULL) +
      theme_roc() +
      theme(
        legend.position  = c(0.715, 0.13),     # inside plot, bottom-center
        legend.direction = "vertical",
        legend.title     = element_blank()
      )
    
    ggsave(file.path(task_dir, "roc_models.png"),  p_overlay,
           width = combo_w, height = combo_h, dpi = combo_dpi)
    ggsave(file.path(task_dir, "roc_models.tiff"), p_overlay,
           width = combo_w, height = combo_h, dpi = combo_dpi, compression = "lzw")
  } else {
    message("No model ROC CSVs found in ", task_dir)
  }
}

message("Done. Per-model PNG/TIFF and combined plots saved in: ",
        file.path(root_dir, "HEALTH"), " and ", file.path(root_dir, "FE"))

