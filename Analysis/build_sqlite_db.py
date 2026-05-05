#!/usr/bin/env python3
"""
build_sqlite_db.py

Convert SGI project TSV/CSV files into a single SQLite database.
Large files (> CHUNK_THRESH rows) are loaded in chunks to avoid memory issues.

Output: SGI_ML_final/report/sgi_features.db
"""

import os
import sqlite3
import time

import pandas as pd

# ============================================================
# CONFIG
# ============================================================
BASE = (
    "/Users/zisanurrahman/Library/CloudStorage/"
    "GoogleDrive-zisan.rr@gmail.com/Other computers/My MacBook Air/"
    "Desktop/Research/SGI/SGI_paper_revision_microbiome/ML_analysis"
)
DATA      = os.path.join(BASE, "data")
FINAL_DIR = os.path.join(BASE, "SGI_ML_final")

DB_PATH    = os.path.join(FINAL_DIR, "report", "sgi_features.db")
CHUNK_THRESH = 500_000   # rows; files larger than this are loaded in chunks
CHUNK_SIZE   = 200_000

# ============================================================
# FILES TO LOAD
# table_name : { path, sep, indexes }
# ============================================================
FILES = {
    # --- core feature matrix (large) ---
    "features": {
        "path"    : os.path.join(DATA, "final_abundance_matrix", "final_master_with_KO_SGI.tsv"),
        "sep"     : "\t",
        "indexes" : ["Sample_ID", "Genome_ID", "Feature_Class"],
    },
    # --- genome-level aggregated matrices ---
    "genome_relabund_genecnt": {
        "path"    : os.path.join(DATA, "sgi_new_features", "sgi_relabund_genecnt.tsv"),
        "sep"     : "\t",
        "indexes" : ["Sample_ID", "Genome_ID", "Feature_Class"],
    },
    "cpm_genome_class": {
        "path"    : os.path.join(DATA, "sgi_new_features", "cpm_sgi_pa_plus_abund.tsv"),
        "sep"     : "\t",
        "indexes" : ["Sample_ID", "Genome_ID", "Feature_Class"],
    },
    "genome_presence_absence": {
        "path"    : os.path.join(DATA, "sgi_new_features", "sgi_presence_absence.tsv"),
        "sep"     : "\t",
        "indexes" : ["Sample_ID", "Genome_ID", "Feature_Class"],
    },
    "genome_functional_gene_counts": {
        "path"    : os.path.join(DATA, "sgi_new_features", "genome_functional_gene_counts.tsv"),
        "sep"     : "\t",
        "indexes" : ["Genome_ID", "Feature_Class"],
    },
    # --- genome coverage (very large) ---
    "genome_coverage_per_sample": {
        "path"    : os.path.join(DATA, "final_abundance_matrix", "genome_coverage_per_sample.tsv"),
        "sep"     : "\t",
        "indexes" : ["Sample_ID", "Genome_ID"],
    },
    "genome_coverage_aggregated": {
        "path"    : os.path.join(DATA, "final_abundance_matrix", "genome_coverage_aggregated.tsv"),
        "sep"     : "\t",
        "indexes" : ["Genome_ID"],
    },
    # --- GCF abundance ---
    "gcf_abundance_cpm": {
        "path"    : os.path.join(DATA, "final_abundance_matrix", "gcf_sample_abundance_long_CPM.csv"),
        "sep"     : ",",
        "indexes" : ["sample", "novel_number", "group"],
    },
    "gcf_abundance_rpkm": {
        "path"    : os.path.join(DATA, "final_abundance_matrix", "gcf_sample_abundance_long_RPKM.csv"),
        "sep"     : ",",
        "indexes" : ["sample", "novel_number", "group"],
    },
    "gcf_stats_cpm": {
        "path"    : os.path.join(DATA, "final_abundance_matrix", "gcf_direction_summary_CPM.csv"),
        "sep"     : ",",
        "indexes" : ["novel_number"],
    },
    "gcf_stats_rpkm": {
        "path"    : os.path.join(DATA, "final_abundance_matrix", "gcf_direction_summary_RPKM.csv"),
        "sep"     : ",",
        "indexes" : ["novel_number"],
    },
    # --- AMP features ---
    "amp_physchem": {
        "path"    : os.path.join(DATA, "final_abundance_matrix", "AMP_physchem_features.tsv"),
        "sep"     : "\t",
        "indexes" : [],
    },
    # --- genome gene counts ---
    "genome_gene_counts_bakta": {
        "path"    : os.path.join(DATA, "final_abundance_matrix", "genome_total_gene_counts_bakta.tsv"),
        "sep"     : "\t",
        "indexes" : ["Genome_ID"],
    },
    # --- metadata ---
    "sample_farm_mapping": {
        "path"    : os.path.join(DATA, "metadata", "sample_farm_mapping.tsv"),
        "sep"     : "\t",
        "indexes" : ["Sample_ID"],
    },
    "sample_province_mapping": {
        "path"    : os.path.join(DATA, "metadata", "sample_province_mapping.tsv"),
        "sep"     : "\t",
        "indexes" : [],
    },
    "taxonomy": {
        "path"    : os.path.join(DATA, "metadata", "gtdbtk.bac120.summary.tsv"),
        "sep"     : "\t",
        "indexes" : [],
    },
}

# ============================================================
# HELPERS
# ============================================================
def n_rows(path, sep):
    """Quick row count without loading the file."""
    with open(path, "r") as f:
        return sum(1 for _ in f) - 1  # subtract header


def load_chunked(path, sep, con, table, chunk_size=CHUNK_SIZE):
    """Load a large file in chunks using if_exists='append'."""
    first = True
    total = 0
    for chunk in pd.read_csv(path, sep=sep, chunksize=chunk_size,
                              low_memory=False, on_bad_lines="warn"):
        chunk.columns = [c.strip() for c in chunk.columns]
        chunk.to_sql(table, con, if_exists="replace" if first else "append",
                     index=False)
        first = False
        total += len(chunk)
        print(f"    loaded {total:,} rows...", end="\r")
    print()
    return total


def load_small(path, sep, con, table):
    df = pd.read_csv(path, sep=sep, low_memory=False, on_bad_lines="warn")
    df.columns = [c.strip() for c in df.columns]
    df.to_sql(table, con, if_exists="replace", index=False)
    return len(df)


def add_indexes(con, table, cols):
    for col in cols:
        idx_name = f"idx_{table}_{col}"
        try:
            con.execute(f"CREATE INDEX IF NOT EXISTS {idx_name} ON {table} ({col})")
        except sqlite3.OperationalError as e:
            print(f"    [warn] index {idx_name}: {e}")
    con.commit()


# ============================================================
# MAIN
# ============================================================
def main():
    os.makedirs(os.path.dirname(DB_PATH), exist_ok=True)

    print(f"Database: {DB_PATH}\n")
    con = sqlite3.connect(DB_PATH)
    con.execute("PRAGMA journal_mode=WAL")
    con.execute("PRAGMA synchronous=NORMAL")
    con.execute("PRAGMA cache_size=-1000000")   # ~1 GB cache

    for table, cfg in FILES.items():
        path = cfg["path"]
        sep  = cfg["sep"]
        idxs = cfg["indexes"]

        if not os.path.exists(path):
            print(f"[SKIP] {table} — file not found: {path}")
            continue

        print(f"[{table}]  {os.path.basename(path)}")
        t0 = time.time()

        nrows = n_rows(path, sep)
        print(f"  rows in file: {nrows:,}")

        if nrows > CHUNK_THRESH:
            n = load_chunked(path, sep, con, table)
        else:
            n = load_small(path, sep, con, table)

        print(f"  loaded {n:,} rows in {time.time()-t0:.1f}s")

        if idxs:
            print(f"  indexing: {idxs}")
            add_indexes(con, table, idxs)

        print()

    # Summary
    print("=" * 50)
    print("Tables in database:")
    for (tbl,) in con.execute("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"):
        (cnt,) = con.execute(f"SELECT COUNT(*) FROM {tbl}").fetchone()
        print(f"  {tbl:<35} {cnt:>12,} rows")

    con.close()
    size_mb = os.path.getsize(DB_PATH) / 1e6
    print(f"\nDatabase size: {size_mb:.1f} MB")
    print(f"Done: {DB_PATH}")


if __name__ == "__main__":
    main()
