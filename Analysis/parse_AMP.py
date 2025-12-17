#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 10:50:15 2025

@author: zisan
"""


"""
python plot_amp_density_bubble_violin.py \
  --per_genome_tsv "AMP_density/amp_density_per_genome.tsv" \
  --taxonomy_tsv "/Volumes/Expansion 1/New_transfer/SGI_genomes/Combined_SGI_PiBac/GTDB-Tk_out_on_filtered_SGI_plus_PiBac/gtdbtk.bac120.summary.tsv" \
  --outdir "AMP_density_plots" \
  --taxon_level genus \
  --plots bubble,violin,box \
  --min_genomes_per_group 3 \
  --top_groups_bubble 40

Bubble + Violin plots for AMP density, using amp_density.py outputs.

Input: amp_density_per_genome.tsv (from amp_density.py)
Expected columns: Genome, n_cAMPs, Length_bp, rho_AMP, CI95_lo, CI95_hi, Phylum, Genus, Species
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, re, argparse, numpy as np, pandas as pd, matplotlib.pyplot as plt

try:
    import seaborn as sns
    sns.set_context("talk"); sns.set_style("whitegrid")
except Exception:
    sns = None

# -------------------- Utils --------------------

def makedirs(p): os.makedirs(p, exist_ok=True)

def parse_rank(s, prefix):
    m = re.search(prefix + r"([^;]+)", str(s))
    return m.group(1) if m else np.nan

def load_density(path):
    df = pd.read_csv(path, sep="\t")
    for c in ["rho_AMP", "n_cAMPs", "Length_bp"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    if "Genome" not in df.columns:
        raise ValueError("Density table must have a 'Genome' column.")
    return df

def load_and_parse_gtdb(path):
    tx = pd.read_csv(path, sep="\t")
    ucol = "user_genome" if "user_genome" in tx.columns else (tx.columns[0])
    # If classification present, parse ranks
    if "classification" in tx.columns:
        tx["Phylum"]  = tx["classification"].apply(lambda s: parse_rank(s, "p__"))
        tx["Genus"]   = tx["classification"].apply(lambda s: parse_rank(s, "g__"))
        tx["Species"] = tx["classification"].apply(lambda s: parse_rank(s, "s__"))
    out = tx.rename(columns={ucol: "Genome"})
    keep_cols = [c for c in ["Genome","Phylum","Genus","Species"] if c in out.columns]
    if not keep_cols:
        raise ValueError("Could not find taxonomy columns (Phylum/Genus/Species) or 'classification' in GTDB file.")
    return out[keep_cols].drop_duplicates()

def merge_taxonomy(density_df, tax_df):
    df = density_df.merge(tax_df, on="Genome", how="left")
    for c in ["Phylum","Genus","Species"]:
        if c not in df.columns: df[c] = np.nan
    return df

def mode_or_first(series):
    vc = series.dropna().value_counts()
    return vc.index[0] if len(vc) else np.nan

def nonempty(df, cols):
    return df.dropna(subset=cols) if all(c in df.columns for c in cols) else pd.DataFrame()

# -------------------- Plotters --------------------

def bubble_plot(df, taxon_col, out_png, out_pdf, topN=40, color_by_phylum=True):
    base = nonempty(df, [taxon_col, "rho_AMP"])
    if base.empty:
        print(f"[WARN] No data for bubble plot at {taxon_col}.")
        return False

    stats = (base.groupby(taxon_col, as_index=False)
                  .agg(mean_rho=("rho_AMP","mean"),
                       sd_rho=("rho_AMP","std"),
                       n_genomes=("rho_AMP","count"),
                       prevalence=("n_cAMPs", lambda x: np.mean(pd.to_numeric(x, errors="coerce").fillna(0) >= 1))))
    if stats.empty:
        print(f"[WARN] Empty stats at {taxon_col}.")
        return False

    stats["se_rho"]  = stats["sd_rho"] / np.sqrt(stats["n_genomes"].clip(lower=1))
    stats["ci95_lo"] = stats["mean_rho"] - 1.96*stats["se_rho"]
    stats["ci95_hi"] = stats["mean_rho"] + 1.96*stats["se_rho"]

    # add representative Phylum for coloring
    if color_by_phylum and "Phylum" in df.columns:
        rep_phylum = (df.dropna(subset=[taxon_col])
                        .groupby(taxon_col)["Phylum"]
                        .apply(mode_or_first)
                        .reset_index(name="Phylum"))
        stats = stats.merge(rep_phylum, on=taxon_col, how="left")

    gsub = stats.sort_values("mean_rho", ascending=False).head(topN).reset_index(drop=True)
    if gsub.empty:
        print(f"[WARN] After sorting/topN, no rows to plot at {taxon_col}.")
        return False

    plt.figure(figsize=(10, 7))
    if sns is not None:
        ax = sns.scatterplot(
            data=gsub, x="prevalence", y="mean_rho",
            size="n_genomes",
            hue=("Phylum" if "Phylum" in gsub.columns else None),
            sizes=(60, 900), alpha=0.85, edgecolor="black", linewidth=0.5
        )
    else:
        ax = plt.gca()
        sizes = 60 + 840*(gsub["n_genomes"]-gsub["n_genomes"].min())/max(1,(gsub["n_genomes"].max()-gsub["n_genomes"].min()))
        for i, r in gsub.iterrows():
            ax.scatter(r["prevalence"], r["mean_rho"], s=float(sizes.iloc[i]), alpha=0.85)

    # y-error bars
    ax.errorbar(gsub["prevalence"], gsub["mean_rho"],
                yerr=[gsub["mean_rho"]-gsub["ci95_lo"], gsub["ci95_hi"]-gsub["mean_rho"]],
                fmt="none", ecolor="gray", elinewidth=1, capsize=2, alpha=0.7, zorder=0)

    # annotate top few
    for _, r in gsub.head(10).iterrows():
        ax.text(r["prevalence"]+0.005, r["mean_rho"], str(r[taxon_col]), fontsize=9, va="center")

    ax.set_xlabel("Share of genomes with ≥1 cAMP (prevalence)")
    ax.set_ylabel("Mean AMP density (ρAMP)")
    ax.set_title(f"AMP density by {taxon_col}: bubble plot")
    ax.set_xlim(-0.02, 1.02)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300); plt.savefig(out_pdf); plt.close()
    print(f"[OK] Saved bubble plot ({taxon_col}): {out_png}")
    return True

def violin_plot(df, taxon_col, out_png, out_pdf, min_genomes=3):
    base = nonempty(df, [taxon_col, "rho_AMP"])
    if base.empty:
        print(f"[WARN] No data for violins at {taxon_col}.")
        return False
    counts = base[taxon_col].value_counts()
    keep = counts[counts >= min_genomes].index
    vdf = base[base[taxon_col].isin(keep)].copy()
    if vdf.empty:
        print(f"[WARN] No {taxon_col} meet ≥{min_genomes} genomes; skipping violins.")
        return False

    order = (vdf.groupby(taxon_col)["rho_AMP"].median()
                .sort_values(ascending=False).index.tolist())

    plt.figure(figsize=(max(10, 0.35*len(order)), 6))
    if sns is not None:
        ax = sns.violinplot(data=vdf, x=taxon_col, y="rho_AMP", order=order, inner="quartile", cut=0)
        sns.stripplot(data=vdf, x=taxon_col, y="rho_AMP", order=order, color="k", alpha=0.35, jitter=0.25, size=3)
    else:
        ax = plt.gca()
        data = [vdf.loc[vdf[taxon_col]==g, "rho_AMP"].values for g in order]
        ax.boxplot(data); ax.set_xticks(range(1,len(order)+1)); ax.set_xticklabels(order, rotation=45, ha="right")

    ax.set_xlabel(f"{taxon_col} (≥{min_genomes} genomes)")
    ax.set_ylabel("AMP density (ρAMP)")
    ax.set_title(f"Per-genome AMP density by {taxon_col} (violin)")
    plt.xticks(rotation=45, ha="right"); plt.tight_layout()
    plt.savefig(out_png, dpi=300); plt.savefig(out_pdf); plt.close()
    print(f"[OK] Saved violin plot ({taxon_col}): {out_png}")
    return True

def box_plot(df, taxon_col, out_png, out_pdf, min_genomes=3):
    base = nonempty(df, [taxon_col, "rho_AMP"])
    if base.empty:
        print(f"[WARN] No data for {taxon_col} boxplot.")
        return False
    counts = base[taxon_col].value_counts()
    keep = counts[counts >= min_genomes].index
    bdf = base[base[taxon_col].isin(keep)].copy()
    if bdf.empty:
        print(f"[WARN] No {taxon_col} with ≥{min_genomes} genomes; skipping boxplot.")
        return False

    order = (bdf.groupby(taxon_col)["rho_AMP"].median()
                 .sort_values(ascending=False).index.tolist())

    plt.figure(figsize=(max(10, 0.35*len(order)), 6))
    if sns is not None:
        ax = sns.boxplot(data=bdf, x=taxon_col, y="rho_AMP", order=order, fliersize=5)
        sns.stripplot(data=bdf, x=taxon_col, y="rho_AMP", order=order,
                      color="k", alpha=0.35, jitter=0.25, size=3)
    else:
        ax = plt.gca()
        data = [bdf.loc[bdf[taxon_col]==g, "rho_AMP"].values for g in order]
        ax.boxplot(data); ax.set_xticks(range(1,len(order)+1)); ax.set_xticklabels(order, rotation=45, ha="right")

    ax.set_xlabel(f"{taxon_col} (≥{min_genomes} genomes)")
    ax.set_ylabel("AMP density (ρAMP)")
    ax.set_title(f"Per-genome AMP density by {taxon_col} (box)")
    plt.xticks(rotation=45, ha="right"); plt.tight_layout()
    plt.savefig(out_png, dpi=300); plt.savefig(out_pdf); plt.close()
    print(f"[OK] Saved box plot ({taxon_col}): {out_png}")
    return True

# -------------------- Main --------------------

def main():
    ap = argparse.ArgumentParser(
        description="Generate bubble, violin, and box plots for AMP density with taxonomy merge."
    )
    ap.add_argument("--per_genome_tsv", required=True, help="amp_density_per_genome.tsv from amp_density.py")
    ap.add_argument("--taxonomy_tsv", required=True, help="GTDB summary (e.g., gtdbtk.bac120.summary.tsv)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--taxon_level", choices=["genus","phylum","species"], default="genus",
                    help="Taxonomic level to plot (default: genus)")
    ap.add_argument("--plots", default="bubble,violin,box",
                    help="Comma-separated list of plots to generate: bubble,violin,box (default: all)")
    ap.add_argument("--min_genomes_per_group", type=int, default=3,
                    help="Min genomes required for a group in violin/box (default: 3)")
    ap.add_argument("--top_groups_bubble", type=int, default=40,
                    help="Top N groups to show in bubble plot (by mean rho; default: 40)")
    args = ap.parse_args()

    makedirs(args.outdir); figdir = os.path.join(args.outdir, "figures"); makedirs(figdir)

    dens = load_density(args.per_genome_tsv)
    tax  = load_and_parse_gtdb(args.taxonomy_tsv)
    df   = merge_taxonomy(dens, tax)

    # Select taxon column
    level_map = {"genus": "Genus", "phylum": "Phylum", "species": "Species"}
    taxon_col = level_map[args.taxon_level]
    print(f"[INFO] Using taxon level: {args.taxon_level} → column '{taxon_col}'")

    plots = {p.strip().lower() for p in args.plots.split(",") if p.strip()}
    if not plots:
        print("[WARN] No plots requested via --plots; nothing to do.")
        return

    # Bubble
    if "bubble" in plots:
        ok = bubble_plot(
            df, taxon_col,
            os.path.join(figdir, f"bubble_{args.taxon_level}_density.png"),
            os.path.join(figdir, f"bubble_{args.taxon_level}_density.pdf"),
            topN=args.top_groups_bubble,
            color_by_phylum=True
        )
        if not ok and args.taxon_level != "phylum":
            print("[INFO] Falling back to phylum bubble…")
            bubble_plot(
                df, "Phylum",
                os.path.join(figdir, "bubble_phylum_density.png"),
                os.path.join(figdir, "bubble_phylum_density.pdf"),
                topN=min(args.top_groups_bubble, 30),
                color_by_phylum=False
            )

    # Violin
    if "violin" in plots:
        ok_v = violin_plot(
            df, taxon_col,
            os.path.join(figdir, f"violin_{args.taxon_level}_density.png"),
            os.path.join(figdir, f"violin_{args.taxon_level}_density.pdf"),
            min_genomes=args.min_genomes_per_group
        )
        if not ok_v and args.taxon_level != "phylum":
            print("[INFO] Falling back to phylum violin…")
            violin_plot(
                df, "Phylum",
                os.path.join(figdir, "violin_phylum_density.png"),
                os.path.join(figdir, "violin_phylum_density.pdf"),
                min_genomes=args.min_genomes_per_group
            )

    # Box
    if "box" in plots:
        ok_b = box_plot(
            df, taxon_col,
            os.path.join(figdir, f"box_{args.taxon_level}_density.png"),
            os.path.join(figdir, f"box_{args.taxon_level}_density.pdf"),
            min_genomes=args.min_genomes_per_group
        )
        if not ok_b and args.taxon_level != "phylum":
            print("[INFO] Falling back to phylum box…")
            box_plot(
                df, "Phylum",
                os.path.join(figdir, "box_phylum_density.png"),
                os.path.join(figdir, "box_phylum_density.pdf"),
                min_genomes=args.min_genomes_per_group
            )

if __name__ == "__main__":
    main()
