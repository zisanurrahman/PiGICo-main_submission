#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 12:34:15 2025
@author: zisan
"""


import pandas as pd
import glob
import os
import sys

def parse_rgi_tsv(filepath, genome_id):
    try:
        df = pd.read_csv(filepath, sep='\t', encoding='utf-8', engine='python')

        # Rename columns based on actual RGI output
        df = df.rename(columns={
            'Contig': 'Contig_ID',
            'Start': 'Start',
            'Stop': 'Stop',
            'Best_Hit_ARO': 'CARD_Hit_ARO',
            'Drug Class': 'Drug Class',
            'Resistance Mechanism': 'Resistance Mechanism',
            'AMR Gene Family':'AMR Gene Family',
            'Antibiotic': 'Antibiotic'
        })

        required_cols = ['Contig_ID','Start', 'Stop', 'CARD_Hit_ARO', 'Drug Class','Resistance Mechanism', 'AMR Gene Family', 'Antibiotic']
        missing = [col for col in required_cols if col not in df.columns]
        if missing:
            print(f"⚠️ Missing expected columns in {filepath}: {missing}")
            return pd.DataFrame()

        df_subset = df[required_cols].copy()
        df_subset['Genome_ID'] = genome_id
        return df_subset

    except Exception as e:
        print(f"❌ Error parsing {filepath}: {e}")
        return pd.DataFrame()

# --- Define input/output paths ---
rgi_output_dir = "/Volumes/Expansion/New_transfer/SGI_genomes/Combined_SGI_PiBac/resistance_gene_identifier_SGI"
output_path = "/Users/zisan/Library/CloudStorage/Dropbox/SGI_Publication/Final_Analysis/Results/Adundance_Data/master_rgi_annotations.tsv"

# --- Parse all .txt.txt files ---
all_rgi_annotations = []
for rgi_file in glob.glob(os.path.join(rgi_output_dir, "*.txt.txt")):
    genome_id = os.path.basename(rgi_file).replace(".txt.txt", "")
    print(f"📄 Parsing rgi for: {genome_id}")
    parsed_df = parse_rgi_tsv(rgi_file, genome_id)
    if not parsed_df.empty:
        all_rgi_annotations.append(parsed_df)

# --- Export master table ---
if all_rgi_annotations:
    master_rgi_df = pd.concat(all_rgi_annotations, ignore_index=True)
    master_rgi_df.to_csv(output_path, sep='\t', index=False)
    print(f"\n✅ Master rgi table created: {output_path}")
    print(master_rgi_df.head())
else:
    print("\n⚠️ No valid rgi files parsed.")
    sys.exit()
