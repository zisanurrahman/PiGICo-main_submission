#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ETEC enrichment (v1.2) — supports genome-level abundance matrices.

Abundance example:
  Genome_ID  Sample_ID  Group  Feature_Class  Abundance
  PiBac_015  SRR...     LFE    AMP            1537600.7
  PiBac_015  SRR...     LFE    BGC            48.07
  ...

Assumptions:
- Each row is a genome's abundance for a given Feature_Class in a Sample.
- We label a GENOME as ETEC-active if:
   * BGC: its row in the BGC predictions has max(anti_gram_negative_*) ≥ cutoff
   * AMP: its row in the AMP table has AMP_probability ≥ cutoff  (proxy)
- For each sample & class, we count/sum ONLY genomes flagged as active.
- Group can come from the abundance file's "Group" column (Healthy/PWD/LFE/HFE).
  If you also provide a groups file (Healthy/PWD), it will override the group in abundance.

Run example:
python etec_enrichment_from_structures_v1_2.py \
  --abundance "/path/genome_feature_abundance_matrix.tsv" \
  --bgc "/path/combined_probabilities.csv" \
  --amp "/path/AMP_physchem_features.tsv" \
  --prob-cutoff 0.7 \
  --outdir "results_etec_structured" \
  --abund-sample-col Sample_ID \
  --abund-genome-col Genome_ID \
  --abund-class-col Feature_Class \
  --abund-value-col Abundance
"""

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu
import matplotlib.pyplot as plt


# -------------------------
# Helpers
# -------------------------
def std_cols(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [c.strip() for c in df.columns]
    return df

def norm_id(s: pd.Series) -> pd.Series:
    return s.astype(str).str.replace(r'[^A-Za-z0-9]+', '', regex=True).str.upper()

def _find_col(df, candidates):
    low = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in low:
            return low[cand.lower()]
    return None

def detect(df, names, required=True):
    col = _find_col(df, names)
    if col is None and required:
        raise ValueError(f"Missing any of columns: {names}\nAvailable: {list(df.columns)}")
    return col

def load_groups(path: str|None) -> pd.DataFrame|None:
    if path is None:
        return None
    g = pd.read_csv(path, sep=None, engine="python")
    g = std_cols(g)
    s = detect(g, ["sample","sample_id"])
    t = detect(g, ["group","status","health_status"])
    out = g[[s,t]].rename(columns={s:"sample", t:"group"})
    out["group"] = out["group"].astype(str).str.strip().str.lower().map({"healthy":"Healthy","pwd":"PWD"})
    if out["group"].isna().any():
        raise ValueError("Groups file must have Healthy/PWD labels.")
    return out

def load_abundance(path: str|Path,
                   sample_col: str|None,
                   genome_col: str|None,
                   class_col: str|None,
                   value_col: str|None) -> pd.DataFrame:
    df = pd.read_csv(path, sep=None, engine="python")
    df = std_cols(df)
    s = sample_col or detect(df, ["Sample_ID","sample","sample_id"])
    g = genome_col or detect(df, ["Genome_ID","genome","genome_id","filename","Access"])
    c = class_col or detect(df, ["Feature_Class","feature_class","class"])
    v = value_col or detect(df, ["Abundance","abundance","value","rpkm","normalized_abundance"])
    out = df[[s,g,c,v]].rename(columns={s:"sample", g:"genome_id", c:"feature_class", v:"abundance"})
    out["feature_class"] = out["feature_class"].astype(str).str.upper()
    # if group present in abundance, keep it
    grp_col = _find_col(df, ["Group","group"])
    if grp_col:
        grp = df[[s, grp_col]].drop_duplicates().rename(columns={s:"sample", grp_col:"group"})
        out = out.merge(grp, on="sample", how="left")
    out["genome_norm"] = norm_id(out["genome_id"])
    return out

def load_bgc_preds(path: str|Path) -> tuple[pd.DataFrame, list[str]]:
    df = pd.read_csv(path, sep=None, engine="python")
    df = std_cols(df)
    idc = detect(df, ["filename","Genome_ID","genome","id"])
    df["genome_id"] = df[idc].astype(str)
    df["genome_norm"] = norm_id(df["genome_id"])
    gn_cols = [c for c in df.columns if c.lower().startswith("anti_gram_negative_")]
    if not gn_cols:
        raise ValueError("No 'anti_gram_negative_*' columns in BGC table.")
    df["prob_gn_max"] = df[gn_cols].apply(pd.to_numeric, errors="coerce").max(axis=1)
    # allow antibacterial_* as backup
    ab_cols = [c for c in df.columns if c.lower().startswith("antibacterial_")]
    if ab_cols:
        df["prob_ab_max"] = df[ab_cols].apply(pd.to_numeric, errors="coerce").max(axis=1)
        df["prob_any_neg"] = df[["prob_gn_max","prob_ab_max"]].max(axis=1)
    else:
        df["prob_any_neg"] = df["prob_gn_max"]
    return df, gn_cols

def load_amp_table(path: str|Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=None, engine="python")
    df = std_cols(df)
    g = detect(df, ["Genome","genome","Genome_ID","filename","Access","id"])
    p = detect(df, ["AMP_probability","amp_prob","probability"])
    out = df.rename(columns={g:"Genome", p:"AMP_probability"})
    out["genome_id"] = out["Genome"].astype(str)
    out["genome_norm"] = norm_id(out["genome_id"])
    return out

def add_presence(x: pd.Series, thresh: float) -> pd.Series:
    return (x > thresh).astype(int)

def summarize(ps_all: pd.DataFrame, outdir: Path):
    gsum = ps_all.groupby(["what","group"]).agg(
        prop_any=("any_present","mean"),
        mean_count=("n_items","mean"),
        median_count=("n_items","median"),
        mean_abund=("total_abundance","mean"),
        median_abund=("total_abundance","median"),
        n=("sample","nunique"),
    ).reset_index()
    gsum.to_csv(outdir/"group_summary.csv", index=False)

    # stats
    def fisher_any(df, label):
        sub = df[df["what"]==label]
        tab = pd.crosstab(sub["group"], sub["any_present"])
        for g in ["Healthy","PWD"]:
            if g not in tab.index: tab.loc[g] = 0
        for c in [0,1]:
            if c not in tab.columns: tab[c] = 0
        tab = tab.loc[["Healthy","PWD"], [0,1]]
        a,b = tab.loc["Healthy",1], tab.loc["Healthy",0]
        c,d = tab.loc["PWD",1],     tab.loc["PWD",0]
        OR = ((a+0.5)*(d+0.5))/((b+0.5)*(c+0.5))
        _, p = fisher_exact(tab.values, alternative="two-sided")
        return OR, p, tab

    def mwu(df, label, col):
        sub = df[df["what"]==label]
        x = sub.loc[sub["group"]=="Healthy", col].values
        y = sub.loc[sub["group"]=="PWD", col].values
        if len(x)==0 or len(y)==0:
            return np.nan, np.nan
        stat, p = mannwhitneyu(x, y, alternative="two-sided")
        return stat, p

    lines = []
    for label in ["BGC","AMP_proxy","BOTH"]:
        OR, p_f, tab = fisher_any(ps_all, label)
        Uc, pc = mwu(ps_all, label, "n_items")
        Ua, pa = mwu(ps_all, label, "total_abundance")
        lines.append(f"[{label}] Fisher any-present: OR={OR:.3g}, p={p_f:.3g}\n{tab}\n")
        lines.append(f"[{label}] Mann–Whitney counts: U={Uc:.3g}, p={pc:.3g}")
        lines.append(f"[{label}] Mann–Whitney total_abundance: U={Ua:.3g}, p={pa:.3g}\n")
    (outdir/"stats_tests.txt").write_text("\n".join(lines))


# -------------------------
# Main
# -------------------------
def main():
    ap = argparse.ArgumentParser(description="ETEC enrichment using genome-level abundance matrix (v1.2)")
    ap.add_argument("--abundance", required=True)
    ap.add_argument("--bgc", required=True)
    ap.add_argument("--amp", required=True)
    ap.add_argument("--groups", default=None, help="Optional groups file (Healthy/PWD). If omitted, uses Group in abundance if present.")
    ap.add_argument("--prob-cutoff", type=float, default=0.7)
    ap.add_argument("--presence-threshold", type=float, default=0.0)
    ap.add_argument("--outdir", default="results_etec_structured")

    # Abundance column overrides
    ap.add_argument("--abund-sample-col", default=None)
    ap.add_argument("--abund-genome-col", default=None)
    ap.add_argument("--abund-class-col", default=None)
    ap.add_argument("--abund-value-col", default=None)

    args = ap.parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    abund = load_abundance(args.abundance, args.abund_sample_col, args.abund_genome_col, args.abund_class_col, args.abund_value_col)
    bgc_df, _ = load_bgc_preds(args.bgc)
    amp_df = load_amp_table(args.amp)

    # Resolve groups
    if args.groups:
        groups = load_groups(args.groups)
        abund = abund.drop(columns=[c for c in ["group"] if c in abund.columns])
        abund = abund.merge(groups, on="sample", how="left")
    else:
        if "group" not in abund.columns:
            raise ValueError("No groups file provided and no 'Group' column in abundance.")

    # Determine active genome sets
    active_bgc_genomes = set(bgc_df.loc[bgc_df["prob_any_neg"] >= args.prob_cutoff, "genome_norm"])
    active_amp_genomes = set(amp_df.loc[pd.to_numeric(amp_df["AMP_probability"], errors="coerce") >= args.prob_cutoff, "genome_norm"])

    # Flag rows that are active given their class
    a = abund.copy()
    a["is_active"] = False
    a.loc[(a["feature_class"]=="BGC") & (a["genome_norm"].isin(active_bgc_genomes)), "is_active"] = True
    a.loc[(a["feature_class"]=="AMP") & (a["genome_norm"].isin(active_amp_genomes)), "is_active"] = True

    # Build per-sample summaries
    def per_sample(active_df: pd.DataFrame, classes: list[str]):
        rows = []
        for label in classes:
            sub = active_df[(active_df["feature_class"]==label) & (active_df["is_active"])]
            if sub.empty:
                for s in active_df["sample"].unique():
                    rows.append({"sample": s, "what": ("AMP_proxy" if label=="AMP" else label), "n_items": 0, "total_abundance": 0.0, "any_present": 0})
                continue
            sub = sub.copy()
            sub["present"] = (sub["abundance"] > args.presence_threshold).astype(int)
            n_by_sample = sub.groupby("sample").apply(lambda d: d.loc[d["present"]==1, "genome_id"].nunique()).rename("n_items")
            tot_by_sample = sub.groupby("sample")["abundance"].sum().rename("total_abundance")
            any_by_sample = sub.groupby("sample")["present"].max().rename("any_present").astype(int)
            out = pd.concat([n_by_sample, tot_by_sample, any_by_sample], axis=1).reset_index()
            out["what"] = "AMP_proxy" if label=="AMP" else label
            rows.append(out)
        return pd.concat(rows, ignore_index=True)

    ps = per_sample(a, ["BGC","AMP"])

    # BOTH = union across labels
    both = ps.pivot_table(index="sample", columns="what", values=["n_items","total_abundance","any_present"], fill_value=0)
    both_n = (both["n_items"].sum(axis=1)).rename("n_items")
    both_ab = (both["total_abundance"].sum(axis=1)).rename("total_abundance")
    both_any = (both["any_present"].max(axis=1)).rename("any_present").astype(int)
    both_df = pd.concat([both_n, both_ab, both_any], axis=1).reset_index()
    both_df["what"] = "BOTH"

    # Merge back group labels
    sample_groups = abund[["sample"] + (["group"] if "group" in abund.columns else [])].drop_duplicates()
    ps_all = pd.concat([ps, both_df], ignore_index=True).merge(sample_groups, on="sample", how="left")

    # Save per-sample table
    ps_all.to_csv(outdir/"per_sample_counts.csv", index=False)
    pd.DataFrame({"genome_norm": sorted(active_bgc_genomes), "class":"BGC"}).to_csv(outdir/"active_genomes_bgc.csv", index=False)
    pd.DataFrame({"genome_norm": sorted(active_amp_genomes), "class":"AMP_proxy"}).to_csv(outdir/"active_genomes_amp.csv", index=False)

    # Summaries and stats
    summarize(ps_all, outdir)

    # Plots
    prop = ps_all.groupby(["what","group"])["any_present"].mean().reset_index()
    whats = ["BGC","AMP_proxy","BOTH"]
    groups_u = sorted(prop["group"].dropna().unique())
    width = 0.35
    xs = np.arange(len(whats))
    plt.figure(figsize=(7,5))
    for i, grp in enumerate(groups_u):
        vals = [prop[(prop["what"]==w) & (prop["group"]==grp)]["any_present"].mean() if not prop[(prop["what"]==w) & (prop["group"]==grp)].empty else 0 for w in whats]
        plt.bar(xs + (i-0.5)*width, vals, width, label=grp)
        for j,v in enumerate(vals):
            plt.text(xs[j] + (i-0.5)*width, v + 0.02, f"{v:.2f}", ha="center", va="bottom", fontsize=9)
    plt.xticks(xs, whats)
    plt.ylim(0,1)
    plt.ylabel("Proportion with ≥1 active genome")
    plt.xlabel("")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir/"fig_prop_any_active.png", dpi=300)
    plt.close()

    def boxplot_one(df, label, fname, ycol):
        sub = df[df["what"]==label]
        data = [sub[sub["group"]==g][ycol].values for g in groups_u]
        plt.figure(figsize=(7,5))
        plt.boxplot(data, labels=groups_u, showfliers=False)
        plt.ylabel(ycol.replace("_"," "))
        plt.title(f"{label} — {ycol}")
        plt.tight_layout()
        plt.savefig(outdir/fname, dpi=300)
        plt.close()

    for label in whats:
        boxplot_one(ps_all, label, f"fig_counts_box_{label}.png", "n_items")
        boxplot_one(ps_all, label, f"fig_abundance_box_{label}.png", "total_abundance")


if __name__ == "__main__":
    main()

