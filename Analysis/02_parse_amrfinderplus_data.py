#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 12:20:26 2025
@author: zisan
"""

import sys
import pandas as pd
import glob
import os
import matplotlib.pyplot as plt
import seaborn as sns

def parse_amrfinderplus_tsv(filepath, genome_id):
    try:
        df = pd.read_csv(filepath, sep='\t', encoding='utf-8', engine='python')

        # Updated column renaming based on your file's actual header
        column_renames = {
            'Protein id': 'Protein_ID',
            'Contig id': 'Contig_ID',
            'Start': 'Gene_Start',
            'Stop': 'Gene_Stop',
            'Strand': 'Gene_Strand',
            'Element symbol': 'AMR_Gene',
            'Element name': 'AMR_Sequence_Name',
            'Class': 'AMR_Class',
            'Subclass': 'AMR_Subclass',
            'Type': 'Type',
            'Subtype': 'Subtype',
        }
        df = df.rename(columns={k: v for k, v in column_renames.items() if k in df.columns})

        required_cols = ['Protein_ID', 'Contig_ID', 'Gene_Start', 'Gene_Stop', 'Gene_Strand',
                         'AMR_Gene', 'AMR_Sequence_Name', 'AMR_Class', 'AMR_Subclass', 'Type', 'Subtype']
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


# --- Set Paths ---
input_dir = "/Volumes/Expansion/New_transfer/SGI_genomes/Combined_SGI_PiBac/amrfinderplus_output"  # <-- UPDATE THIS with the path to cleaned .tsv files
output_table = "/Users/zisan/Library/CloudStorage/Dropbox/SGI_Publication/Final_Analysis/Results/Adundance_Data/master_amrfinderplus_annotations.tsv"
output_plot_class = "/Users/zisan/Library/CloudStorage/Dropbox/SGI_Publication/Final_Analysis/Plots/amr_class_summary.png"
output_plot_subclass = "/Users/zisan/Library/CloudStorage/Dropbox/SGI_Publication/Final_Analysis/Plots/amr_subclass_summary.png"

# --- Parse All .tsv Files ---
all_amrfinder_annotations = []
for filepath in glob.glob(os.path.join(input_dir, "*.tsv")):
    genome_id = os.path.basename(filepath).split(".")[0]
    print(f"📄 Parsing: {genome_id}")
    parsed_df = parse_amrfinderplus_tsv(filepath, genome_id)
    if not parsed_df.empty:
        all_amrfinder_annotations.append(parsed_df)

# --- Combine & Save ---
if all_amrfinder_annotations:
    combined_df = pd.concat(all_amrfinder_annotations, ignore_index=True)
    combined_df.to_csv(output_table, sep='\t', index=False)
    print(f"\n✅ Master AMRFinderPlus table saved to: {output_table}")
else:
    print("\n⚠️ No valid input found.")
    sys.exit()

# --- Plotting AMR Class Summary ---
plt.figure(figsize=(10, 6))
sns.countplot(data=combined_df, y='AMR_Class', order=combined_df['AMR_Class'].value_counts().index)
plt.title("AMR Class Distribution")
plt.xlabel("Count")
plt.ylabel("AMR Class")
plt.tight_layout()
plt.savefig(output_plot_class)
print(f"📊 AMR class plot saved: {output_plot_class}")
plt.close()

# --- Plotting AMR Subclass Summary ---
plt.figure(figsize=(12, 8))
top_subclasses = combined_df['AMR_Subclass'].value_counts().nlargest(20)
sns.barplot(x=top_subclasses.values, y=top_subclasses.index)
plt.title("Top 20 AMR Subclasses")
plt.xlabel("Count")
plt.ylabel("AMR Subclass")
plt.tight_layout()
plt.savefig(output_plot_subclass)
print(f"📊 AMR subclass plot saved: {output_plot_subclass}")
plt.close()
