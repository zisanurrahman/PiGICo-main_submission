#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# ========== INPUTS ==========
ann_file <- "/Users/zisan/SGI_Paper/BGC_analysis/BiG-SCAPE/v2/output_files/2025-09-05_22-07-12_c0.3/record_annotations.tsv"
tax_file <- "/Users/zisan/SGI_Paper/New_ML_Analysis_Trans/gtdb_taxonomy.csv"
out_csv  <- "/Users/zisan/SGI_Paper/BGC_analysis/BiG-SCAPE/v2/output_files/2025-09-05_22-07-12_c0.3/Category_Count_per_Order.csv"

# ========== BUCKET LOGIC ==========
BUCKETS <- c("NRPS", "PKS", "TERPENE", "RIPP", "HYBRID", "OTHER")

to_bucket <- function(x) {
  x <- tolower(trimws(x %||% ""))
  if (grepl("\\.", x)) return("HYBRID")
  if (x %in% c("nrps", "pks", "terpene", "ripp")) return(toupper(x))
  return("OTHER")
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# ========== LOAD FILES ==========
bgc <- read_tsv(ann_file, show_col_types = FALSE) %>%
  select(Description, Category)

# Extract normalized genome: first two parts of Description
bgc <- bgc %>%
  mutate(
    genome = str_extract(Description, "^([^_]+_[^_]+)")
  )

# Load GTDB taxonomy
tax <- read_csv(tax_file, show_col_types = FALSE) %>%
  select(genome, order)

# ========== JOIN + BUCKET ==========
bgc_tax <- bgc %>%
  left_join(tax, by = "genome") %>%
  filter(!is.na(order)) %>%
  mutate(
    Bucket = vapply(Category, to_bucket, character(1)),
    Bucket = factor(Bucket, levels = BUCKETS)
  )

# ========== COUNT PER ORDER × BUCKET ==========
count_df <- bgc_tax %>%
  count(order, Bucket, name = "Avg_Cat") %>%
  complete(order, Bucket = factor(BUCKETS, levels = BUCKETS), fill = list(Avg_Cat = 0)) %>%
  rename(Order = order, Category = Bucket)

# Add total per Order
count_df <- count_df %>%
  group_by(Order) %>%
  mutate(Avg_total = sum(Avg_Cat)) %>%
  ungroup()

# ========== SAVE ==========
write_csv(count_df, out_csv)
message("✅ Saved: ", out_csv)
