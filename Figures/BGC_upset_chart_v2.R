#!/usr/bin/env Rscript
# ==========================================================
# UpSetR: NRPS / PKS / TERPENE / RIPP / HYBRID / OTHER (ALL CAPS)
# - Multi-label aware (hybrids → multiple sets)
# - HYBRID is separate bucket (if Category contains a dot)
# - OTHER only if none of the other 5 are present
# ==========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(UpSetR)
})

# ---------- user params ----------
input_fp   <- "/Users/zisan/SGI_Paper/BGC_analysis/BiG-SCAPE/v2/output_files/2025-09-05_22-07-12_c0.3/record_annotations.tsv"
out_dir    <- "/Users/zisan/SGI_Paper/Full_paper_clean/Figures"
nintersects <- 40
pdf_w <- 10.5; pdf_h <- 7
tif_w <- 2600; tif_h <- 3050; tif_res <- 600  # pixels/res for crisp TIFF

# ---------- buckets ----------
BUCKETS <- c("NRPS","PKS","TERPENE","RIPP","HYBRID","OTHER")

# ---------- helpers ----------
`%||%` <- function(x, y) if (is.null(x)) y else x
tokenize <- function(s) {
  s <- tolower(trimws(s %||% ""))
  unlist(strsplit(s, "[/;,+:&|\\-_. ]+"), use.names = FALSE)
}
has_token <- function(toks, key) {
  key <- tolower(key)
  any(grepl(paste0("^", key, "$"), toks, ignore.case = TRUE)) ||
    any(grepl(key, toks, fixed = TRUE))
}

# ---------- load ----------
bgc <- read_tsv(input_fp, show_col_types = FALSE) %>%
  mutate(id = Description,
         class = Category %||% "Unknown",
         class = as.character(class))

# ---------- build membership sets ----------
det <- lapply(seq_len(nrow(bgc)), function(i) {
  cls <- bgc$class[i]
  toks <- tokenize(cls)
  present <- set_names(rep(FALSE, length(BUCKETS)), BUCKETS)
  
  # Check main classes
  for (k in c("NRPS","PKS","TERPENE","RIPP")) {
    present[[k]] <- has_token(toks, k)
  }
  
  # HYBRID = if Category contains dot
  present[["HYBRID"]] <- grepl("\\.", cls)
  
  # OTHER = only if none of the other 5 present
  present[["OTHER"]] <- !any(present[BUCKETS[1:5]])
  
  present
})

mem <- do.call(rbind, det) %>% as.data.frame()
mem <- mem %>% mutate(id = bgc$id)
mem[BUCKETS] <- lapply(mem[BUCKETS], as.logical)

# ---------- UpSet input ----------
class_list <- lapply(setNames(nm = BUCKETS), function(k) {
  mem$id[mem[[k]]]
})
up_df <- UpSetR::fromList(class_list) %>% as.data.frame()

# ---------- plot & save ----------
cutoff <- "record_annotations"

tiff(file.path(out_dir, sprintf("fig_classes_upset_%s.tiff", cutoff)),
     width = tif_w, height = tif_h, res = tif_res, compression = "lzw")
upset(
  up_df,
  sets = BUCKETS,
  order.by = "freq",
  nintersects = nintersects,
  mainbar.y.label = "BGC Intersections",
  sets.x.label    = "BGC Class Sizes",
  keep.order = TRUE
)
dev.off()

pdf(file.path(out_dir, sprintf("fig_classes_upset_%s.pdf", cutoff)),
    width = pdf_w, height = pdf_h, useDingbats = FALSE)
class_cols <- c(
  "NRPS"    = "cornflowerblue",
  "PKS"     = "#E6AB02",
  "TERPENE" = "#1B9E77",
  "RIPP"    = "darkslateblue",
  "HYBRID"  = "#D95F02",
  "OTHER"   = "#BEBEBE"
)

upset(
  up_df,
  sets = BUCKETS,
  order.by = "freq",
  nintersects = nintersects,
  mainbar.y.label = "BGC Intersections",
  sets.x.label    = "BGC Class Sizes",
  keep.order = TRUE,
  sets.bar.color = class_cols[BUCKETS],
  text.scale = c(3, 3, 3, 2, 2, 2)
)

dev.off()

message("✅ Saved UpSet plots to: ", out_dir)
