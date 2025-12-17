
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 14 12:35:01 2025

@author: zisan
"""

from Bio import SeqIO
import os
import glob
import pandas as pd

input_dir = "/Volumes/Expansion/New_transfer/SGI_genomes/Combined_SGI_PiBac/gut_pathway_analysis/gbk_files"
output_file = "/Users/zisan/Library/CloudStorage/Dropbox/SGI_Publication/Final_Analysis/Results/Abundance_Data/master_gutsmash_table.tsv"

records = []

for gbk_file in glob.glob(os.path.join(input_dir, "*.gbk")):
    genome_id = os.path.basename(gbk_file).split(".")[0]
    for record in SeqIO.parse(gbk_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS" and 'locus_tag' in feature.qualifiers:
                protein_id = feature.qualifiers['locus_tag'][0]
                records.append({"Protein_ID": protein_id, "Genome_ID": genome_id})

if records:
    df = pd.DataFrame(records)
    df.to_csv(output_file, sep='\t', index=False)
    print(f"✅ gutSMASH mapping saved to {output_file}")
else:
    print("❌ No valid CDS features with locus_tag found.")
