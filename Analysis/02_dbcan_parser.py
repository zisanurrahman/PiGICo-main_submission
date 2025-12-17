#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 15:21:09 2025

@author: zisan
"""

import os
import pandas as pd
import re
from collections import defaultdict

DBCAN_DIR = "/Volumes/Expansion/New_transfer/SGI_genomes/Combined_SGI_PiBac/dbcan_out_filtered_combined"
output_file = "/Users/zisan/Library/CloudStorage/Dropbox/SGI_Publication/Analysis/Data/genome_cazyme_feature_matrix.tsv"

def extract_families(cazy_str):
    if pd.isna(cazy_str) or cazy_str.strip() == "-":
        return []
    cazy_str = re.sub(r'\(.*?\)', '', cazy_str)
    matches = re.findall(r'([A-Z]{2,4}\d+)', cazy_str)
    return matches

def parse_dbcan_genome(genome, filepath):
    df = pd.read_csv(filepath, sep="\t", header=0)
    df.columns = [c.strip() for c in df.columns]

    # Require at least 2 tools to support the prediction
    if not set(["HMMER", "dbCAN_sub", "DIAMOND"]).issubset(df.columns):
        return pd.Series(name=genome)

    df["Support"] = df[["HMMER", "dbCAN_sub", "DIAMOND"]].apply(lambda x: sum(pd.notna(x)), axis=1)
    df = df[df["Support"] >= 2]

    all_fams = []
    for col in ["HMMER", "dbCAN_sub", "DIAMOND"]:
        for val in df[col].dropna():
            all_fams += extract_families(val)

    # Count by prefix (e.g., GH, GT)
    prefix_counts = defaultdict(int)
    for fam in all_fams:
        prefix = re.match(r"([A-Z]{2,4})", fam)
        if prefix:
            prefix_counts[f"Num_{prefix.group(1)}"] += 1

    return pd.Series(prefix_counts, name=genome)

# Process all subdirectories
all_genomes = []
for subdir in os.listdir(DBCAN_DIR):
    sub_path = os.path.join(DBCAN_DIR, subdir)
    if os.path.isdir(sub_path):
        overview_file = os.path.join(sub_path, "overview.txt")
        if os.path.exists(overview_file):
            row = parse_dbcan_genome(subdir, overview_file)
            all_genomes.append(row)

# Concatenate to dataframe
df = pd.DataFrame(all_genomes).fillna(0).astype(int)
df.index.name = "Genome"
df.to_csv(output_file, sep="\t")
print(f"✅ CAZyme feature matrix saved: {output_file}")
