#!/usr/bin/env python3
"""
tsv_to_parquet_batch.py
Author: Zisanur Rahman

Description:
    Convert all TSV files in a directory into individual Parquet files.
Usage:
    python tsv_to_parquet_batch.py --indir /path/to/tsv_dir --outdir /path/to/parquet_dir
"""

import os
import argparse
import pandas as pd
from tqdm import tqdm

def main():
    parser = argparse.ArgumentParser(description="Convert all TSV files in a directory to individual Parquet files.")
    parser.add_argument("--indir", required=True, help="Input directory containing TSV files")
    parser.add_argument("--outdir", required=True, help="Output directory for Parquet files")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    tsv_files = [f for f in os.listdir(args.indir) if f.endswith(".tsv")]
    if not tsv_files:
        print("❌ No TSV files found in the provided directory.")
        return

    print(f"🔍 Found {len(tsv_files)} TSV files. Converting to Parquet...")

    for f in tqdm(tsv_files, desc="Converting files"):
        tsv_path = os.path.join(args.indir, f)
        parquet_path = os.path.join(args.outdir, os.path.splitext(f)[0] + ".parquet")

        try:
            df = pd.read_csv(tsv_path, sep="\t")
            df.to_parquet(parquet_path, index=False)
        except Exception as e:
            print(f"⚠️ Skipping {f} due to error: {e}")

    print(f"✅ Conversion complete. Parquet files saved to: {args.outdir}")

if __name__ == "__main__":
    main()

