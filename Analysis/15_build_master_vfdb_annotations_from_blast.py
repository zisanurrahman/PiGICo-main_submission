#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 12:19:34 2025

@author: zisan
"""

import os
import glob
import pandas as pd

# Input directory containing per-genome VFDB BLAST results
VFDB_DIR = "/Users/zisan/Library/CloudStorage/OneDrive-UniversityofManitoba/Animal Science/SGI/SGI_PiBac_analysis/Data/vfdb_filtered"

records = []

for vfdb_file in glob.glob(os.path.join(VFDB_DIR, "*.tsv")):
    genome_id = os.path.basename(vfdb_file).split(".")[0]  # assumes filename = SGI_870.tsv

    try:
        df = pd.read_csv(vfdb_file, sep="\t")
        if "qseqid" not in df.columns:
            print(f"⚠️ Skipping {vfdb_file}: missing 'qseqid'")
            continue

        df["Contig_ID"] = df["qseqid"]
        df["Genome_ID"] = genome_id
        df["Feature_Class"] = "VFDB"

        records.append(df[["Contig_ID", "Genome_ID", "Feature_Class"]].drop_duplicates())

    except Exception as e:
        print(f"❌ Failed to process {vfdb_file}: {e}")

# Combine and export
if records:
    master_vfdb_df = pd.concat(records, ignore_index=True)
    master_vfdb_df.to_csv("/Users/zisan/Library/CloudStorage/Dropbox/SGI_Publication/Final_Analysis/Results/Abundance_Data/master_vfdb_annotations.tsv", sep="\t", index=False)
    print("✅ Saved: master_vfdb_annotations.tsv")
else:
    print("❌ No valid VFDB annotation files processed.")
