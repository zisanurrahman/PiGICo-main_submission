#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  2 12:32:37 2025

@author: zisan
"""

import pandas as pd
import os, re, sys
from pathlib import Path

# ---- CONFIG: edit paths ----
DIR_IN  = Path("/Volumes/Expansion/New_transfer/SGI_genomes/Combined_SGI_PiBac/dbcan_out_on_Combined/Renamed_overview_files")
OUT_DIR = Path("/Users/zisan/SGI_Paper/Full_paper_clean/Data/Cazyme_analysis")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_LONG = OUT_DIR / "dbcan_cazy_long.tsv"                # Genome, Family, Count
OUT_WIDE = OUT_DIR / "dbcan_cazy_count_wide.tsv"          # wide matrix Genome x Family
OUT_CLASS = OUT_DIR / "dbcan_cazy_class_counts_by_genome.tsv"  # per-class counts

# ---- helpers ----
ENCODINGS = ["utf-8", "utf-8-sig", "latin-1", "cp1252"]
DELIMS = [",", "\t", ";"]

def try_read_table(path):
    last_err = None
    for enc in ENCODINGS:
        for delim in DELIMS:
            try:
                df = pd.read_csv(path, sep=delim, encoding=enc, on_bad_lines="skip", engine="python")
                # Heuristic: require at least ~5 columns to consider parse "ok"
                if df.shape[1] >= 3:
                    return df
            except Exception as e:
                last_err = e
                continue
    raise RuntimeError(f"Failed to read {path} with common encodings/delims. Last error:\n{last_err}")

def normalize_sample_name(fname: str) -> str:
    # Strip trailing _overview.csv/.tsv, keep the base sample id
    base = os.path.basename(fname)
    name = re.sub(r'_overview\.(csv|tsv)$', '', base, flags=re.IGNORECASE)
    return name

def find_col(df, patterns):
    """Find first column whose name matches any regex in patterns (case-insensitive)."""
    cols = list(df.columns)
    for pat in patterns:
        regex = re.compile(pat, re.IGNORECASE)
        for c in cols:
            if regex.search(str(c)):
                return c
    return None

CAZY_REGEX = re.compile(r'\b(GH|GT|CE|PL|AA|CBM)\d+\b', re.IGNORECASE)

def extract_cazy_families(cell):
    """Return list of CAZy families (e.g., GH43, GT2) from a text cell."""
    if pd.isna(cell):
        return []
    s = str(cell)
    # Remove any parentheses content to avoid clutter like "(subfamily x)"
    s = re.sub(r'\(.*?\)', '', s)
    # Find all CAZy family tokens
    hits = CAZY_REGEX.findall(s)  # returns tuples because of groups
    # Re-run to get full tokens (group0)
    full = re.findall(CAZY_REGEX, s)  # still tuples; fallback below
    # Simple fallback: use re.findall with non-capturing group
    full2 = re.findall(r'\b(?:GH|GT|CE|PL|AA|CBM)\d+\b', s, flags=re.IGNORECASE)
    # Normalize to upper-case prefix + digits
    return [h.upper() for h in full2]

def to_long_counts(df, sample_name, hmm_col, dia_col):
    """Produce long rows (Genome, Family, Count=1) per hit row."""
    # Work on a copy to avoid SettingWithCopy warnings
    x = df.copy()

    # Extract families from HMMER and DIAMOND columns
    x["_HMM_FAMS"] = x[hmm_col].apply(extract_cazy_families) if hmm_col else [[]]*len(x)
    x["_DIA_FAMS"] = x[dia_col].apply(extract_cazy_families) if dia_col else [[]]*len(x)

    # Determine if this sample has any HMMER hits at all
    hmm_any = any(len(v) > 0 for v in x["_HMM_FAMS"])

    # Prefer HMMER; if none present anywhere for this genome, use DIAMOND
    if hmm_any:
        x["_FAMS"] = x["_HMM_FAMS"]
    else:
        x["_FAMS"] = x["_DIA_FAMS"]

    # Explode families to long
    y = x[["_FAMS"]].explode("_FAMS").dropna()
    if y.empty:
        return pd.DataFrame(columns=["Genome","Family","Count"])

    y = y.rename(columns={"_FAMS":"Family"})
    y["Genome"] = sample_name
    y["Count"]  = 1
    return y[["Genome","Family","Count"]]

def cazy_class(fam: str) -> str:
    m = re.match(r'^(GH|GT|CE|PL|AA|CBM)', fam, flags=re.IGNORECASE)
    return m.group(1).upper() if m else None

# ---- main ----
def main():
    if not DIR_IN.exists():
        print(f"Error: {DIR_IN} does not exist.", file=sys.stderr)
        sys.exit(1)

    long_rows = []

    files = [p for p in DIR_IN.iterdir() if p.is_file() and p.suffix.lower() in (".csv",".tsv")]
    files.sort()
    if not files:
        print(f"No .csv/.tsv overview files found in {DIR_IN}", file=sys.stderr)
        sys.exit(1)

    for fp in files:
        try:
            df = try_read_table(fp)
        except Exception as e:
            print(f"[WARN] Skipping {fp.name}: {e}", file=sys.stderr)
            continue

        sample = normalize_sample_name(fp.name)

        # Identify likely HMMER/DIAMOND columns
        hmm_col = find_col(df, patterns=[
            r'^dbcan[_\s-]*hmm$', r'^hmm$', r'hmmer', r'hmm\s*hits?', r'hmm\s*result'
        ])
        dia_col = find_col(df, patterns=[
            r'^diamond$', r'diamond\s*hits?', r'blastp', r'diamond\s*result'
        ])

        # If neither present, try to salvage by looking for any column with CAZy tokens
        if hmm_col is None and dia_col is None:
            # find any column containing CAZy tokens in at least 1 row
            for c in df.columns:
                vals = df[c].astype(str).head(200)
                if any(CAZY_REGEX.search(v) for v in vals):
                    # pick first as HMMER stand-in
                    hmm_col = c
                    break

        # Still nothing? skip the file
        if hmm_col is None and dia_col is None:
            print(f"[WARN] {fp.name}: no HMMER/DIAMOND-like columns found; skipping.", file=sys.stderr)
            continue

        # Ensure missing columns exist for uniform code path
        if hmm_col is None:
            df["__HMM_EMPTY__"] = ""
            hmm_col = "__HMM_EMPTY__"
        if dia_col is None:
            df["__DIA_EMPTY__"] = ""
            dia_col = "__DIA_EMPTY__"

        # Build long
        part = to_long_counts(df, sample, hmm_col, dia_col)
        if not part.empty:
            long_rows.append(part)

    if not long_rows:
        print("No CAZy families extracted from any file.", file=sys.stderr)
        sys.exit(1)

    long_df = pd.concat(long_rows, ignore_index=True)

    # Collapse to counts per Genome×Family
    long_df = (long_df
               .groupby(["Genome","Family"], as_index=False)["Count"]
               .sum()
               .sort_values(["Genome","Family"]))

    # Save long
    long_df.to_csv(OUT_LONG, sep="\t", index=False)

    # Wide matrix
    wide = long_df.pivot_table(index="Genome", columns="Family", values="Count",
                               aggfunc="sum", fill_value=0).reset_index()
    wide.to_csv(OUT_WIDE, sep="\t", index=False)

    # Per-class counts by genome (useful sanity check + optional plot)
    tmp = long_df.copy()
    tmp["Class"] = tmp["Family"].apply(cazy_class)
    class_counts = (tmp.groupby(["Genome","Class"], as_index=False)["Count"]
                    .sum()
                    .sort_values(["Genome","Class"]))
    class_counts.to_csv(OUT_CLASS, sep="\t", index=False)

    print(f"OK ✓  Wrote:\n  - {OUT_LONG}\n  - {OUT_WIDE}\n  - {OUT_CLASS}")

if __name__ == "__main__":
    main()

