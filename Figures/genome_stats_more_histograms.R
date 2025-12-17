#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggtern)
  library(plotly)
  library(patchwork)
  library(htmlwidgets)
})

# ====== INPUTS ======
seqkit_tsv <- "/Users/zisan/SGI_Paper/Full_paper_clean/Data/genome_stats_out/assembly_stats_with_group.tsv"
checkm_tsv <- "/Volumes/Expansion/New_transfer/SGI_genomes/Combined_SGI_PiBac/drep_SGI_plus_PiBac_genomes/data/checkM/checkM_outdir/results.tsv"
outdir     <- "/Users/zisan/SGI_Paper/Full_paper_clean/Data/genome_stats_out/plots_genome_with_stats"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Find a GC column and coerce to numeric safely (works for "GC(%)", "GC %", etc.)
get_gc_numeric <- function(d) {
  # candidate column names containing 'GC' but not Q20/Q30
  cand <- grep("(^|[^A-Za-z])GC([^A-Za-z]|$)", names(d), ignore.case = TRUE, value = TRUE)
  cand <- setdiff(cand, grep("^Q(20|30)", names(d), value = TRUE))
  if (length(cand) == 0) stop("Could not find a GC column in seqkit table. Names: ", paste(names(d), collapse=", "))
  col <- cand[1]
  x <- d[[col]]
  if (is.numeric(x)) as.numeric(x) else readr::parse_number(as.character(x))
}


# QC thresholds (edit as needed)
max_contam_for_filter <- 8     # filter out >8% contamination from all plots
qc_min_completeness   <- 90    # pass if completeness >= 90
qc_max_contamination  <- 8     # pass if contamination <= 5

# ====== LOAD & JOIN ======
# Seqkit
df_base <- read_tsv(seqkit_tsv, show_col_types = FALSE) %>%
  mutate(
    basename  = basename(file),
    genome_mb = sum_len / 1e6,
    GC_num    = get_gc_numeric(cur_data_all()),     # <-- add GC here
    Group     = factor(Group, levels = c("SGI","PiBac","Other"))
  ) %>%
  select(basename, Group, genome_mb, GC_num)

# CheckM (sanitize headers)
chk <- read_tsv(checkm_tsv, show_col_types = FALSE)
names(chk) <- gsub("[^A-Za-z0-9]+","_",names(chk))
chk <- chk %>%
  transmute(
    basename       = Bin_Id,
    completeness   = as.numeric(Completeness),
    contamination  = as.numeric(Contamination)
  )

dfj <- df_base %>%
  left_join(chk, by = "basename") %>%
  filter(Group %in% c("SGI","PiBac"),
         is.finite(genome_mb), is.finite(GC_num),
         is.finite(completeness), is.finite(contamination)) %>%
  # Global quality filter: keep only <= 8% contamination
  filter(contamination <= max_contam_for_filter) %>%
  droplevels()


# ====== SCALE to simplex (A+B+C=1) ======
scale01 <- function(x){
  r <- range(x, na.rm = TRUE)
  if (diff(r) == 0) return(rep(0.5, length(x)))
  (x - r[1]) / (r[2] - r[1])
}

dfx <- dfj %>%
  mutate(
    A_sz = scale01(genome_mb),       # Size
    B_cp = scale01(completeness),    # Completeness
    C_ct = scale01(contamination)    # Contamination (higher = worse; as requested)
  ) %>%
  mutate(S = A_sz + B_cp + C_ct,
         A = A_sz / S,
         B = B_cp / S,
         C = C_ct / S) %>%
  filter(is.finite(A) & is.finite(B) & is.finite(C))

col_map <- c("SGI" = "#1B9E77", "PiBac" = "darkslateblue")

# ====== 1) INTERACTIVE TERNARY (hover shows RAW values) ======
p_tern <- ggtern(
  dfx,
  aes(x = A, y = B, z = C, color = Group,
      text = paste0(
        basename, "\n",
        "Size: ", round(genome_mb, 2), " Mb\n",
        "Completeness: ", round(completeness, 1), "%\n",
        "Contamination: ", round(contamination, 2), "%"
      ))
) +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_manual(values = col_map, guide = guide_legend(title = NULL)) +
  Tlab("Size (Scaled)") + Llab("Completeness (Scaled) ") + Rlab("Contamination (Scaled)") +
  theme_bw(base_size = 12) + theme_showarrows() +
  theme(legend.title = element_blank(), legend.key.size = unit(0.45,"cm")) +
  labs(title = "Ternary: Size vs Completeness vs Contamination")

p_tern_int <- plotly::ggplotly(p_tern, tooltip = "text")
htmlwidgets::saveWidget(p_tern_int, file.path(outdir, "ternary_interactive.html"))

# ====== 2) CORNER ANNOTATIONS + RAW MEDIANS in caption ======
corner_labels <- tibble(
  x = c(1, 0, 0.02),
  y = c(0, 1, 0.02),
  z = c(0, 0, 0.96),
  lab = c("\u2191 Size (Mb)", "\u2191 Completeness (%)", "\u2191 Contamination (%)")
)

summ <- dfx %>%
  group_by(Group) %>%
  summarise(
    med_size = median(genome_mb, na.rm=TRUE),
    med_comp = median(completeness, na.rm=TRUE),
    med_cont = median(contamination, na.rm=TRUE),
    .groups="drop"
  )

cap_txt <- paste(
  paste0("SGI median: ", round(summ$med_size[summ$Group=="SGI"],2)," Mb, ",
         round(summ$med_comp[summ$Group=="SGI"],1),"%, ",
         round(summ$med_cont[summ$Group=="SGI"],2),"% contam."),
  paste0("PiBac median: ", round(summ$med_size[summ$Group=="PiBac"],2)," Mb, ",
         round(summ$med_comp[summ$Group=="PiBac"],1),"%, ",
         round(summ$med_cont[summ$Group=="PiBac"],2),"% contam."),
  sep = "\n"
)

p_annot <- p_tern +
  geom_text(data = corner_labels,
            aes(x = x, y = y, z = z, label = lab),
            inherit.aes = FALSE, size = 3.5, fontface = "italic") +
  labs(caption = cap_txt)
print(p_annot)
ggsave(file.path(outdir, "ternary_annotated.pdf"), p_annot, width = 8, height = 8.0, dpi = 600)

# ====== 3) QC THRESHOLD VIEW (pass/fail) ======
dfx <- dfx %>%
  mutate(pass_qc = completeness >= qc_min_completeness & contamination <= qc_max_contamination)

p_thresh <- ggtern(dfx, aes(A, B, C)) +
  geom_point(aes(color = pass_qc, shape = Group), size = 2, alpha = 0.95) +
  scale_color_manual(values = c(`TRUE` = "#2ca02c", `FALSE` = "#d62728"),
                     labels = c("Fail", "Pass")) +
  Tlab("Size (scaled)") + Llab("Completeness (scaled)") + Rlab("Contamination (scaled)") +
  theme_bw(base_size = 12) + theme_showarrows() +
  theme(legend.title = element_blank(), legend.key.size = unit(0.45, "cm")) +
  labs(
    title = "Ternary with QC thresholds",
    subtitle = paste0("Pass: Completeness ≥ ", qc_min_completeness,
                      "% & Contamination ≤ ", qc_max_contamination, "%")
  )
print(p_thresh)

ggsave(file.path(outdir, "ternary_qc_thresholds.png"), p_thresh, width = 6.8, height = 6.0, dpi = 300)

# ====== 4) PAIR WITH RAW HISTOGRAMS ======
# ====== 4) PAIR WITH RAW HISTOGRAMS ======
# dfx already exists; ensure Group is just SGI/PiBac and carry GC_num over
dfx <- dfj %>%
  mutate(Group = factor(Group, levels = c("SGI","PiBac")))

col_map <- c("SGI" = "#1B9E77", "PiBac" = "darkslateblue")

# GC histogram (combined)
# ===== Add medians =====
add_median_line <- function(p, data, var, col_map) {
  meds <- data %>%
    group_by(Group) %>%
    summarise(med = median(.data[[var]], na.rm = TRUE), .groups="drop")
  
  p + geom_vline(data = meds, aes(xintercept = med, color = Group),
                 linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
    geom_text(data = meds, aes(x = med, y = 0,
                               label = paste0("Median = ", round(med, 2))),
              angle = 90, vjust = -0.5, hjust = -0.2, size = 2.5, inherit.aes = FALSE) +
    scale_color_manual(values = col_map)
}

p_gc <- ggplot(dfx, aes(x = GC_num, fill = Group)) +
  geom_histogram(bins = 30, alpha = 0.8, position = "identity") +
  scale_fill_manual(values = col_map) +
  labs(x = "GC Content (%)", y = "Genomes (Count)") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA))
p_gc <- add_median_line(p_gc, dfx, "GC_num", col_map)
print(p_gc)
ggsave(file.path(outdir, "p_gc.png"), p_gc, width = 2.5, height = 2.5, dpi = 600)

# Size / Completeness / Contamination histograms (combined)
p_size <- ggplot(dfx, aes(genome_mb, fill = Group)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = col_map) +
  labs(x = "Size (Mb)", y = "Genomes (Count)") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none", panel.border = element_rect(color="black", fill=NA))
p_size  <- add_median_line(p_size,  dfx, "genome_mb",    col_map)
print(p_size)
ggsave(file.path(outdir, "p_size.png"), p_size, width = 2.5, height = 2.5, dpi = 600)

p_comp <- ggplot(dfx, aes(completeness, fill = Group)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = col_map) +
  labs(x = "Completeness (%)", y = "Genomes (Count)") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none", panel.border = element_rect(color="black", fill=NA))
p_comp  <- add_median_line(p_comp,  dfx, "completeness", col_map)
print(p_comp)
ggsave(file.path(outdir, "p_comp.png"), p_comp, width = 2.5, height = 2.5, dpi = 600)
ggsave(file.path(outdir, "p_comp.pdf"), p_comp, width = 2.5, height = 2.5, dpi = 600)

p_cont <- ggplot(dfx, aes(contamination, fill = Group)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = col_map) +
  labs(x = "Contamination (%)", y = "Genomes (Count)") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none", panel.border = element_rect(color="black", fill=NA))
p_cont  <- add_median_line(p_cont,  dfx, "contamination", col_map)
print(p_cont)
ggsave(file.path(outdir, "p_cont.png"), p_cont, width = 2.5, height = 2.5, dpi = 600)
ggsave(file.path(outdir, "p_cont.pdf"), p_cont, width = 2.5, height = 2.5, dpi = 600)
# Panel with ternary + three histograms (use GC instead of completeness if you prefer)
#layout_final <- (p_annot | (p_size / p_comp / p_cont)) + patchwork::plot_layout(widths = c(2, 1))
#ggsave(file.path(outdir, "ternary_with_raw_dists.png"), layout_final, width = 10, height = 6.2, dpi = 300)

# ---------------------------------------------------------------
# Save per-group histograms with median line + numeric label
# bins=30, alpha=0.7, position="identity", median color = group color
# ---------------------------------------------------------------
save_hist_per_group <- function(d, outdir, col_map, varcol, xlab,
                                file_prefix, height = 2.2, width = 2.1, dpi = 600) {
  stopifnot(is.character(varcol), length(varcol) == 1, varcol %in% names(d))
  grps <- levels(d$Group)
  
  for (g in grps) {
    sub <- d %>% filter(Group == g)
    if (nrow(sub) == 0) next
    
    # Copy target column to a fixed name to avoid tidy-eval issues
    sub <- sub %>% mutate(.val = .data[[varcol]])
    med_val <- median(sub$.val, na.rm = TRUE)
    
    # Build histogram (bins=30 etc.)
    p <- ggplot(sub, aes(x = .val)) +
      geom_histogram(bins = 30, alpha = 0.7, position = "identity",
                     fill = col_map[g]) +
      labs(title = paste0(g, " genomes"), x = xlab, y = "Genomes (Count)") +
      theme_classic(base_size = 10) +
      theme(panel.border = element_rect(color = "black", fill = NA))
    
    # Find ymax from the built plot so we can place the label above bars
    gb <- ggplot_build(p)
    ys <- unlist(lapply(gb$data, function(dd) if (!is.null(dd$ymax)) dd$ymax else dd$count))
    ymax <- max(ys, na.rm = TRUE); if (!is.finite(ymax)) ymax <- 1
    
    # Add median line (group color) + numeric label near top
    p <- p +
      geom_vline(xintercept = med_val, color = col_map[g],
                 linetype = "dashed", linewidth = 0.7) +
      annotate("text",
               x = med_val, y = ymax * 0.2,
               label = paste0("Median = ", round(med_val, 2)),
               angle = 90, vjust = -0.5, hjust = -0.2, size = 2.7)
    
    # Save file
    ggsave(file.path(outdir, paste0(file_prefix, "_", g, "_median.png")),
           p, width = width, height = height, dpi = dpi)
  }
}

col_map <- c("SGI" = "#1B9E77", "PiBac" = "darkslateblue")

save_hist_per_group(dfx, outdir, col_map, "genome_mb",    "Size (Mb)",         "hist_size")
save_hist_per_group(dfx, outdir, col_map, "completeness", "Completeness (%)",  "hist_completeness")
save_hist_per_group(dfx, outdir, col_map, "contamination","Contamination (%)", "hist_contamination")
save_hist_per_group(dfx, outdir, col_map, "GC_num",       "GC Content (%)",    "hist_gc")


message("Saved:")
message("  - ", file.path(outdir, "ternary_interactive.html"), "  (open in browser)")
message("  - ", file.path(outdir, "ternary_annotated.png"))
message("  - ", file.path(outdir, "ternary_qc_thresholds.png"))
message("  - ", file.path(outdir, "ternary_with_raw_dists.png"))
