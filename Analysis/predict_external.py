#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
predict_external.py
===================
Apply saved v6/v7 model bundles to an independent external dataset.

Workflow
--------
1. Load one or more .joblib model bundles (saved by v6 or v7).
2. Read the external dataset in the same long-format TSV as training data.
3. Build the feature matrix using the SAME granularity and feature classes
   as the training run.
4. Align columns to the training feature set (fill missing = 0, drop extras).
5. Predict probabilities and class labels.
6. If true labels are available: compute AUROC, accuracy, confusion matrix.
7. Save all outputs as CSVs + PNG plots for R.

Usage examples
--------------
# Single model, evaluate with true labels:
python predict_external.py \
  --models   results/v6_1_CPM/saved_models/model_Healthy_vs_PWD_RF.joblib \
  --genomics data/external/external_master.tsv \
  --abundance_type CPM \
  --outdir   results/external_validation

# All models from a run, no true labels (deploy mode):
python predict_external.py \
  --models   results/v6_1_CPM/saved_models/*.joblib \
  --genomics data/external/external_master.tsv \
  --abundance_type CPM \
  --outdir   results/external_validation \
  --no_labels

# v7 combo-specific model:
python predict_external.py \
  --models   results/v7_RPKM/all_farms/AMP+BGC/saved_models/model_Healthy_vs_PWD_RF.joblib \
  --genomics data/external/external_master.tsv \
  --abundance_type RPKM \
  --outdir   results/external_validation
"""

import argparse
import sys
import warnings
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.metrics import (
    roc_auc_score,
    roc_curve,
    confusion_matrix,
    accuracy_score,
    classification_report,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")


# =============================================================================
# DATA LOADING  (same logic as v6/v7)
# =============================================================================

def read_long_format(path: str, abundance_type: str = "RPKM") -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df.columns = [c.strip() for c in df.columns]
    col = abundance_type.upper()
    if col not in df.columns:
        raise SystemExit(
            f"[FATAL] Column '{col}' not found. Available: {list(df.columns)}"
        )
    df = df.rename(columns={col: "Abundance"})
    other = "CPM" if col == "RPKM" else "RPKM"
    if other in df.columns:
        df = df.drop(columns=[other])
    if "Group" in df.columns:
        df["Group"] = df["Group"].astype(str).str.strip()
    df["Abundance"] = pd.to_numeric(df["Abundance"], errors="coerce").fillna(0.0)
    df["Sample_ID"] = df["Sample_ID"].astype(str).str.strip()
    return df


def clr_transform(X: pd.DataFrame, pseudocount: float = 1e-6) -> pd.DataFrame:
    """Centered log-ratio transform for compositional data."""
    X_ps = X.values + pseudocount
    log_X = np.log(X_ps)
    geom_mean = log_X.mean(axis=1, keepdims=True)
    return pd.DataFrame(log_X - geom_mean, index=X.index, columns=X.columns)


def make_feature_matrix(df: pd.DataFrame, feature_classes: list,
                         granularity: str = "feature_class") -> pd.DataFrame:
    sub = df[df["Feature_Class"].isin(feature_classes)].copy()
    if sub.empty:
        return pd.DataFrame()
    if granularity == "feature_class":
        mat = (sub.groupby(["Sample_ID", "Feature_Class"])["Abundance"]
                  .sum().unstack(fill_value=0.0))
    elif granularity == "genome_class":
        sub["_feat"] = (sub["Genome_ID"].astype(str) + "|" +
                        sub["Feature_Class"].astype(str))
        mat = (sub.groupby(["Sample_ID", "_feat"])["Abundance"]
                  .sum().unstack(fill_value=0.0))
    else:
        raise ValueError(f"Unknown granularity: {granularity}")
    return mat.reindex(sorted(mat.columns), axis=1)


def get_labels(df: pd.DataFrame) -> pd.Series:
    return (df[["Sample_ID", "Group"]]
            .drop_duplicates()
            .set_index("Sample_ID")["Group"])


# =============================================================================
# PREDICTION HELPERS
# =============================================================================

def _predict_proba(model, X: pd.DataFrame) -> np.ndarray:
    if hasattr(model, "predict_proba"):
        p = model.predict_proba(X)
        return p[:, 1] if p.ndim == 2 else p
    if hasattr(model, "decision_function"):
        d = model.decision_function(X)
        return 1.0 / (1.0 + np.exp(-d))
    clf = model.named_steps.get("clf") if hasattr(model, "named_steps") else None
    if clf is not None and hasattr(clf, "decision_function"):
        d = model.decision_function(X)
        return 1.0 / (1.0 + np.exp(-d))
    return model.predict(X).astype(float)


def align_to_training(X_new: pd.DataFrame, feature_cols: list) -> pd.DataFrame:
    """Reindex new data to exactly the columns seen during training."""
    missing = set(feature_cols) - set(X_new.columns)
    extra   = set(X_new.columns) - set(feature_cols)
    if missing:
        print(f"  [INFO] {len(missing)} training features absent in new data "
              f"=> filled with 0: {sorted(missing)}")
    if extra:
        print(f"  [INFO] {len(extra)} extra features in new data => ignored")
    return X_new.reindex(columns=feature_cols, fill_value=0.0)


# =============================================================================
# PLOTS
# =============================================================================

def plot_roc(fpr, tpr, auc_val, title, outpath):
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.plot([0, 1], [0, 1], "k--", lw=1)
    ax.plot(fpr, tpr, lw=2, label=f"AUROC = {auc_val:.3f}")
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title(title)
    ax.legend(loc="lower right", fontsize=9)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()


def plot_confusion(cm, labels, title, outpath):
    fig, ax = plt.subplots(figsize=(4, 3.5))
    im = ax.imshow(cm, cmap="Blues")
    ax.set_xticks([0, 1]); ax.set_xticklabels(labels, fontsize=10)
    ax.set_yticks([0, 1]); ax.set_yticklabels(labels, fontsize=10)
    for i in range(2):
        for j in range(2):
            ax.text(j, i, str(cm[i, j]), ha="center", va="center",
                    fontsize=14, fontweight="bold",
                    color="white" if cm[i, j] > cm.max() / 2 else "black")
    ax.set_xlabel("Predicted"); ax.set_ylabel("Actual")
    ax.set_title(title, fontsize=10, fontweight="bold")
    plt.colorbar(im, ax=ax)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()


# =============================================================================
# MAIN PREDICTION FUNCTION
# =============================================================================

def predict_with_bundle(
    bundle: dict,
    X_new: pd.DataFrame,
    y_new: pd.Series | None,
    outdir: Path,
) -> dict:
    """
    Apply one model bundle to the external feature matrix.
    Returns a dict with predictions and metrics.
    """
    model       = bundle["model"]
    feat_cols   = bundle["feature_cols"]
    mdl_name    = bundle["model_name"]
    task_tag    = bundle["task_tag"]
    pos_class   = bundle.get("pos_class", "Healthy")
    combo       = bundle.get("combo", "")
    farm_label  = bundle.get("farm_label", "")

    label_parts = [task_tag, mdl_name]
    if combo:
        label_parts.append(combo)
    if farm_label:
        label_parts.append(farm_label)
    label = "_".join(label_parts)

    print(f"\n  Model : {mdl_name}  |  task: {task_tag}"
          + (f"  |  combo: {combo}" if combo else "")
          + (f"  |  farm: {farm_label}" if farm_label else ""))
    print(f"  Training features : {len(feat_cols)}")

    # X_new is already aligned to feat_cols by the caller
    Xf = X_new
    print(f"  Aligned shape     : {Xf.shape}")

    # Predict
    y_prob = _predict_proba(model, Xf)
    y_pred = (y_prob >= 0.5).astype(int)

    result = {
        "model_name":  mdl_name,
        "task_tag":    task_tag,
        "combo":       combo,
        "farm_label":  farm_label,
        "n_samples":   len(Xf),
        "AUROC":       None,
        "Accuracy":    None,
    }

    # Predictions CSV
    pred_df = pd.DataFrame({
        "Sample_ID": Xf.index,
        "y_prob":    y_prob,
        "y_pred":    y_pred,
        "y_pred_label": [pos_class if p == 1 else "PWD" for p in y_pred],
    })

    # Evaluate if labels are available
    if y_new is not None:
        common = pred_df["Sample_ID"][pred_df["Sample_ID"].isin(y_new.index)]
        if len(common) == 0:
            print("  [WARN] No Sample_IDs overlap between predictions and labels.")
        else:
            y_true_bin = (y_new.loc[common] == pos_class).astype(int).values
            y_prob_ev  = y_prob[pred_df["Sample_ID"].isin(common)]
            y_pred_ev  = y_pred[pred_df["Sample_ID"].isin(common)]

            pred_df["y_true"] = np.nan
            mask = pred_df["Sample_ID"].isin(common)
            pred_df.loc[mask, "y_true"] = y_true_bin

            if len(np.unique(y_true_bin)) < 2:
                print("  [WARN] Only one class in labels — AUROC undefined.")
                auc = float("nan")
            else:
                auc = roc_auc_score(y_true_bin, y_prob_ev)
                fpr, tpr, _ = roc_curve(y_true_bin, y_prob_ev)
                plot_roc(fpr, tpr, auc,
                         title=f"External validation — {mdl_name} ({combo})",
                         outpath=outdir / f"roc_{label}.png")

            acc  = accuracy_score(y_true_bin, y_pred_ev)
            cm   = confusion_matrix(y_true_bin, y_pred_ev)
            plot_confusion(cm, ["PWD", pos_class],
                           title=f"CM — {mdl_name} ({combo})",
                           outpath=outdir / f"cm_{label}.png")

            report = classification_report(
                y_true_bin, y_pred_ev,
                target_names=["PWD", pos_class], output_dict=True
            )
            pd.DataFrame(report).transpose().to_csv(
                outdir / f"classification_report_{label}.csv"
            )

            result["AUROC"]    = float(auc)
            result["Accuracy"] = float(acc)
            result["N_pos"]    = int(y_true_bin.sum())
            result["N_neg"]    = int((1 - y_true_bin).sum())
            result["CM_TN"]    = int(cm[0, 0])
            result["CM_FP"]    = int(cm[0, 1])
            result["CM_FN"]    = int(cm[1, 0])
            result["CM_TP"]    = int(cm[1, 1])

            print(f"  AUROC    : {auc:.3f}")
            print(f"  Accuracy : {acc:.3f}")
            print(f"  CM       : TN={cm[0,0]} FP={cm[0,1]} FN={cm[1,0]} TP={cm[1,1]}")

    pred_df.to_csv(outdir / f"predictions_{label}.csv", index=False)
    return result


# =============================================================================
# MAIN
# =============================================================================

def main():
    ap = argparse.ArgumentParser(
        description="Apply saved v6/v7 model bundles to an external dataset"
    )
    ap.add_argument("--models",          required=True, nargs="+",
                    help="Path(s) to .joblib model bundle(s)")
    ap.add_argument("--genomics",        required=True,
                    help="External dataset — same long-format TSV as training")
    ap.add_argument("--abundance_type",  default="RPKM", choices=["RPKM", "CPM"])
    ap.add_argument("--outdir",          required=True,
                    help="Output directory for predictions and metrics")
    ap.add_argument("--pos_class",       default="Healthy",
                    help="Positive class label (default: Healthy)")
    ap.add_argument("--neg_class",       default="PWD",
                    help="Negative class label (default: PWD)")
    ap.add_argument("--no_labels",       action="store_true",
                    help="Set if the external data has no Group/label column")
    ap.add_argument("--feature_classes", default=None,
                    help="Override feature classes (default: use bundle's feature_cols "
                         "to infer). Comma-separated, e.g. AMP,BGC")
    ap.add_argument("--clr", action="store_true",
                    help="Apply CLR transform to feature matrix (required if model was "
                         "trained with --clr)")
    ap.add_argument("--intersect_features", action="store_true",
                    help="Retrain on features present in BOTH training and external data, "
                         "eliminating zero-padding of absent features")
    ap.add_argument("--training_data", default=None,
                    help="Path to original training TSV (required with --intersect_features)")
    ap.add_argument("--exclude_farms", default=None,
                    help="Comma-separated farm IDs to exclude when retraining "
                         "(e.g. Farm_41)")
    ap.add_argument("--farm_col", default="Farm",
                    help="Column name for farm/study in training data (default: Farm)")
    ap.add_argument("--farm_meta", default=None,
                    help="TSV with Sample_ID + farm column (used when training data "
                         "lacks a farm column)")
    ap.add_argument("--topk", type=int, default=0,
                    help="Keep top-k highest-variance features when retraining "
                         "(0 = keep all intersection features, default: 0)")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Resolve glob patterns in --models
    model_paths = []
    for pattern in args.models:
        p = Path(pattern)
        if p.exists():
            model_paths.append(p)
        else:
            matches = sorted(Path(".").glob(pattern))
            model_paths.extend(matches)
    model_paths = [p for p in model_paths if p.exists()]

    if not model_paths:
        raise SystemExit(f"[FATAL] No .joblib files found: {args.models}")
    print(f"[predict] Found {len(model_paths)} model bundle(s)")

    # Load external data
    print(f"[predict] Loading external data: {args.genomics}")
    df = read_long_format(args.genomics, abundance_type=args.abundance_type)

    y_new = None
    if not args.no_labels and "Group" in df.columns:
        y_new = get_labels(df)
        y_new = y_new[y_new.isin([args.pos_class, args.neg_class])]
        print(f"[predict] Labels found: {y_new.value_counts().to_dict()}")
    else:
        print("[predict] Running in deploy mode — no true labels.")

    all_results = []

    for bundle_path in model_paths:
        print(f"\n[predict] Loading bundle: {bundle_path.name}")
        bundle = joblib.load(bundle_path)

        # Determine granularity and feature classes from bundle
        granularity    = bundle.get("granularity", "feature_class")
        feat_cols      = bundle["feature_cols"]

        # Infer feature classes from the training column names
        if args.feature_classes:
            feature_classes = [c.strip() for c in args.feature_classes.split(",")]
        else:
            if granularity == "feature_class":
                feature_classes = feat_cols  # cols ARE the feature classes
            else:
                # genome_class: extract unique Feature_Class from col names (after |)
                feature_classes = list({
                    c.split("|")[1] for c in feat_cols if "|" in c
                })

        # Build external feature matrix (raw counts)
        X_new_raw = make_feature_matrix(df, feature_classes, granularity=granularity)
        if X_new_raw.empty:
            print(f"  [WARN] Empty feature matrix for {bundle_path.name}. Skipping.")
            continue

        if args.intersect_features:
            if not args.training_data:
                raise SystemExit("[FATAL] --intersect_features requires --training_data")

            # ---- Feature intersection: only columns present in BOTH datasets ----
            overlap_cols = [c for c in feat_cols if c in X_new_raw.columns]
            n_drop = len(feat_cols) - len(overlap_cols)
            print(f"  [INFO] Feature intersection: {len(overlap_cols)}/{len(feat_cols)} "
                  f"overlap ({n_drop} training-only features dropped)")
            if len(overlap_cols) == 0:
                print(f"  [WARN] No overlapping features. Skipping.")
                continue

            # ---- Reload training data and rebuild feature matrix ----
            print(f"  [INFO] Loading training data: {args.training_data}")
            df_train = read_long_format(args.training_data,
                                        abundance_type=args.abundance_type)
            X_train = make_feature_matrix(df_train, feature_classes,
                                          granularity=granularity)
            X_train = X_train.reindex(columns=overlap_cols, fill_value=0.0)

            # ---- Training labels ----
            pos_class = bundle.get("pos_class", args.pos_class)
            y_all = get_labels(df_train)
            y_all = y_all[y_all.isin([pos_class, args.neg_class])]
            common_tr = X_train.index.intersection(y_all.index)
            X_train = X_train.loc[common_tr]
            y_train  = (y_all.loc[common_tr] == pos_class).astype(int)

            # ---- Optionally exclude farms ----
            if args.exclude_farms:
                excl = {f.strip() for f in args.exclude_farms.split(",")}
                farm_col = args.farm_col
                # Try training data first, then external farm_meta file
                if farm_col in df_train.columns:
                    farm_src = (df_train[["Sample_ID", farm_col]]
                                .drop_duplicates().set_index("Sample_ID")[farm_col])
                elif args.farm_meta:
                    fm = pd.read_csv(args.farm_meta, sep="\t")
                    farm_src = fm.set_index("Sample_ID")[farm_col]
                else:
                    farm_src = None
                    print(f"  [WARN] Farm column not found and no --farm_meta provided "
                          f"— farm exclusion skipped")
                if farm_src is not None:
                    keep = ~farm_src.reindex(X_train.index).isin(excl)
                    X_train = X_train.loc[keep]
                    y_train  = y_train.loc[keep]
                    print(f"  [INFO] Excluded {excl}: {keep.sum()} training samples remain")

            # ---- Apply CLR to training data (same space) ----
            if args.clr:
                X_train = clr_transform(X_train)

            # ---- Variance-based topk feature selection (on training data) ----
            if args.topk and args.topk > 0:
                variances = X_train.var(axis=0)
                topk_cols = variances.nlargest(min(args.topk, len(variances))).index.tolist()
                X_train = X_train[topk_cols]
                overlap_cols = topk_cols
                print(f"  [INFO] topk={args.topk}: kept {len(topk_cols)} features by variance")

            # ---- Retrain with same hyperparams as saved model ----
            orig = bundle["model"]
            if hasattr(orig, "named_steps"):
                clf_step = orig.named_steps.get(
                    "clf", list(orig.named_steps.values())[-1])
            else:
                clf_step = orig
            new_model = Pipeline([("scaler", StandardScaler()), ("clf", clone(clf_step))])
            print(f"  [INFO] Retraining on {X_train.shape[0]} samples "
                  f"x {X_train.shape[1]} features ...")
            new_model.fit(X_train, y_train)

            # ---- Prepare external data (intersection + topk features only) ----
            X_new = X_new_raw.reindex(columns=overlap_cols, fill_value=0.0)
            if args.clr:
                X_new = clr_transform(X_new)
            print(f"  [INFO] External matrix: {X_new.shape[0]} samples x {X_new.shape[1]} features")

            # Override bundle with retrained model and new feature list
            bundle = dict(bundle)
            bundle["model"]        = new_model
            bundle["feature_cols"] = overlap_cols

        else:
            # ---- Standard path: CLR on external feature space, then align ----
            # CLR is computed over the features actually present in the external
            # dataset before zero-padding absent training features.
            X_new = X_new_raw
            if args.clr:
                print("  [INFO] Applying CLR transform...")
                X_new = clr_transform(X_new)
            X_new = align_to_training(X_new, feat_cols)

        # Predict ALL samples; labels used only for post-hoc evaluation
        y_eval = None
        if y_new is not None:
            y_eval = y_new.reindex(X_new.index).dropna()

        result = predict_with_bundle(bundle, X_new, y_eval, outdir)
        result["bundle_file"] = bundle_path.name
        all_results.append(result)

    # Master results CSV
    if all_results:
        master = pd.DataFrame(all_results)
        master_path = outdir / "EXTERNAL_VALIDATION_SUMMARY.csv"
        master.to_csv(master_path, index=False)
        print(f"\n[predict] Summary saved: {master_path}")

        if "AUROC" in master.columns and master["AUROC"].notna().any():
            print("\n[predict] AUROC summary:")
            print(master[["model_name", "combo", "farm_label", "AUROC", "Accuracy"]]
                  .sort_values("AUROC", ascending=False)
                  .to_string(index=False))

    print(f"\n[predict] All outputs saved to: {outdir}")


if __name__ == "__main__":
    main()
