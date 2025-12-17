#!/usr/bin/env python3

import os
import pandas as pd
from pathlib import Path

# Input and output paths
base_dir = Path("/Volumes/Expansion/New_transfer/SGI_genomes/Combined_SGI_PiBac/virulencefinder_output")
output_csv = "/Users/zisan/SGI_Paper/Full_paper_clean/Data/vfdb/combined_virulencefinder_results.csv"

# List to collect DataFrames
combined_df_list = []

# Traverse each genome subdirectory
for genome_dir in base_dir.iterdir():
    if genome_dir.is_dir():
        genome_id = genome_dir.name
        result_file = genome_dir / "results_tab.tsv"
        
        if result_file.exists():
            try:
                df = pd.read_csv(result_file, sep="\t")
                df["Genome"] = genome_id
                combined_df_list.append(df)
            except Exception as e:
                print(f"Failed to read {result_file}: {e}")

# Concatenate all
if combined_df_list:
    combined_df = pd.concat(combined_df_list, ignore_index=True)
    combined_df.to_csv(output_csv, index=False)
    print(f"✅ Combined file saved to: {output_csv}")
else:
    print("❌ No result_tab.tsv files found.")

