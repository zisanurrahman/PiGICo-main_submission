#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build genome-level abundance matrix from RPKM files
Consolidates AMR_rgi + AMR_AmrPlus → AMR
Enforces 6 feature rows (AMP, AMR, BGC, CAZyme, VFDB, gutSMASH)
per Genome_ID × Sample_ID (missing → 0).
"""

import os, re, glob
import pandas as pd
from tqdm import tqdm

# -------- Paths --------
RPKM_DIR = "/Users/zisan/SGI_Paper/results"
OUT_DIR  = "/Users/zisan/SGI_Paper/results/feature_matrix_from_raw"
os.makedirs(OUT_DIR, exist_ok=True)

OUT_TSV      = os.path.join(OUT_DIR, "genome_feature_abundance_matrix.tsv")
SUMMARY_TSV  = os.path.join(OUT_DIR, "summary_feature_class.tsv")
FE_META      = "/Users/zisan/SGI_Paper/results/Raw_input/fe_Metadata.csv"

# regex to pull genome IDs from contigs
RGX_GENOME = re.compile(r"(PiBac_\d+|SGI_\d+)", re.I)
RGX_SRR    = re.compile(r"(SRR\d+)", re.I)

# consolidated feature set
FEATURE_CLASSES = ["AMP", "AMR", "BGC", "CAZyme", "VFDB", "gutSMASH"]

# ------------------------
def load_fe_meta(path):
    if not os.path.exists(path):
        return {}
    df = pd.read_csv(path, sep=None, engine="python")
    sra = next((c for c in df.columns if "sra" in c.lower()), None)
    fe  = next((c for c in df.columns if "feed" in c.lower()), None)
    if not sra or not fe: 
        return {}
    return dict(zip(df[sra].astype(str).str.upper(),
                    df[fe].astype(str).str.lower()))

def infer_group(fp, fe_map, sample_id):
    low = fp.lower()
    if "guelph" in low: 
        return "PWD"
    elif "pigcat" in low: 
        return "Healthy"
    elif "tan_and_quan_fe1" in low or "fe1" in low:
        text = fp + "_" + sample_id
        m = RGX_SRR.search(text)
        if m:
            fe = fe_map.get(m.group(1).upper(),"")
            if fe.startswith("high"): return "HFE"
            if fe.startswith("low"):  return "LFE"
        return "Unknown"
    else:
        return "Unknown"

def detect_cols(fp):
    hdr = pd.read_csv(fp, sep="\t", nrows=0, engine="c")
    cols = [str(c).strip() for c in hdr.columns]
    contig_col = next((c for c in cols if c.lower().startswith("contig")), None)
    rpkm_col   = next((c for c in cols if c.lower()=="rpkm"), None)
    if not contig_col or not rpkm_col:
        raise KeyError(f"{fp}: no Contig/RPKM columns (got {cols})")
    return contig_col, rpkm_col

def parse_one(fp):
    contig_col, rpkm_col = detect_cols(fp)
    df = pd.read_csv(fp, sep="\t", engine="c", usecols=[contig_col, rpkm_col])
    df.columns = ["Contig","RPKM"]
    df["RPKM"] = pd.to_numeric(df["RPKM"], errors="coerce")
    df = df.dropna(subset=["RPKM"])
    df["Genome_ID"] = df["Contig"].astype(str).str.extract(RGX_GENOME, expand=False)
    df = df.dropna(subset=["Genome_ID"])
    return df

def infer_feature_class(fp):
    low = fp.lower()
    if "_caz" in low:      return "CAZyme"
    if "_gut" in low:      return "gutSMASH"
    if "_rgi" in low:      return "AMR"
    if "_amrplus" in low:  return "AMR"
    if "_vfdb" in low:     return "VFDB"
    if "_amp" in low:      return "AMP"
    return "BGC"   # default

# ------------------------
def main():
    fe_map = load_fe_meta(FE_META)
    files = glob.glob(os.path.join(RPKM_DIR, "**", "*.rpkm.tsv"), recursive=True)
    print(f"Found {len(files)} RPKM files")

    out_rows = []
    for fp in tqdm(files, desc="Aggregating genome×feature"):
        try:
            df = parse_one(fp)
        except Exception as e:
            print(f"❌ {fp}: {e}")
            continue
        sample_id = os.path.basename(fp).replace(".rpkm.tsv","")
        group     = infer_group(fp, fe_map, sample_id)
        feature   = infer_feature_class(fp)

        # sum per genome
        agg = df.groupby("Genome_ID", as_index=False)["RPKM"].sum()
        agg["Sample_ID"] = sample_id
        agg["Group"] = group
        agg["Feature_Class"] = feature
        agg.rename(columns={"RPKM":"Abundance"}, inplace=True)
        out_rows.append(agg)

    if not out_rows:
        print("⚠️ Nothing processed")
        return

    final = pd.concat(out_rows, ignore_index=True)

    # ---- Enforce full 6-feature rows per genome×sample ----
    combos = final[["Genome_ID","Sample_ID","Group"]].drop_duplicates()
    cartesian = (
        combos.assign(key=1)
        .merge(pd.DataFrame({"Feature_Class":FEATURE_CLASSES,"key":1}), on="key")
        .drop("key", axis=1)
    )
    final = cartesian.merge(final,
                            on=["Genome_ID","Sample_ID","Group","Feature_Class"],
                            how="left")
    final["Abundance"] = final["Abundance"].fillna(0)

    final.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"✅ Saved main matrix: {OUT_TSV}")

    # ---- Summary ----
    summary = final.groupby("Feature_Class").agg(
        rows=("Genome_ID","size"),
        unique_genomes=("Genome_ID","nunique"),
        unique_samples=("Sample_ID","nunique"),
        total_abundance=("Abundance","sum"),
        mean_abundance=("Abundance","mean")
    ).reset_index()
    summary.to_csv(SUMMARY_TSV, sep="\t", index=False)
    print(f"📊 Saved summary: {SUMMARY_TSV}")
    print(summary)

if __name__ == "__main__":
    main()
