#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SGI Microbiome Machine-Learning Pipeline
==============================================
  Models: Logistic Regression (LR), Ridge Regression (Ridge),
          Random Forest (RF), XGBoost (XGB)

Usage example:
  python sgi_ml_pipeline_v6_slim.py \
    --genomics  data/merged_rpkm_cpm_with_metadata_long.tsv \
    --outdir    results/v6_CPM \
    --abundance_type CPM \
    --feature_classes AMP,AMR,BGC,CAZyme,VFDB,gutSMASH \
    --study_meta data/sample_farm_mapping.tsv \
    --skip_combo_search
"""

import argparse
import itertools
import sys
import warnings
from datetime import datetime
from pathlib import Path

import joblib

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import mutual_info_classif
from sklearn.inspection import permutation_importance as sk_perm_imp
from sklearn.linear_model import LogisticRegression, RidgeClassifier
from sklearn.metrics import (
    average_precision_score,
    confusion_matrix,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import (
    GridSearchCV,
    StratifiedKFold,
    train_test_split,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

try:
    from xgboost import XGBClassifier
    HAS_XGB = True
except ImportError:
    HAS_XGB = False
    print("[WARN] xgboost not installed — XGB model will be skipped.", file=sys.stderr)

try:
    import shap
    HAS_SHAP = True
except ImportError:
    HAS_SHAP = False

try:
    import seaborn as sns
    HAS_SNS = True
except ImportError:
    HAS_SNS = False

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# =============================================================================
# DATA LOADING
# =============================================================================

def read_long_format(path: str, abundance_type: str = "RPKM") -> pd.DataFrame:
    """Read long-format TSV containing both RPKM and CPM columns."""
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df.columns = [c.strip() for c in df.columns]

    col = abundance_type.upper()
    if col not in df.columns:
        raise SystemExit(
            f"[FATAL] Column '{col}' not found.  Available columns: {list(df.columns)}"
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


def get_labels(df: pd.DataFrame) -> pd.Series:
    return (
        df[["Sample_ID", "Group"]]
        .drop_duplicates()
        .set_index("Sample_ID")["Group"]
    )


def make_feature_matrix(
    df: pd.DataFrame,
    feature_classes=None,
    granularity: str = "feature_class",
) -> pd.DataFrame:
    """Pivot long format → sample × feature matrix."""
    if feature_classes is not None:
        df = df[df["Feature_Class"].isin(feature_classes)].copy()
    if df.empty:
        return pd.DataFrame()

    if granularity == "feature_class":
        mat = (
            df.groupby(["Sample_ID", "Feature_Class"])["Abundance"]
            .sum()
            .unstack(fill_value=0.0)
        )
    elif granularity == "genome_class":
        df = df.copy()
        df["_feat"] = df["Genome_ID"].astype(str) + "|" + df["Feature_Class"].astype(str)
        mat = (
            df.groupby(["Sample_ID", "_feat"])["Abundance"]
            .sum()
            .unstack(fill_value=0.0)
        )
    else:
        raise ValueError(f"Unknown granularity: {granularity}")

    return mat.reindex(sorted(mat.columns), axis=1)


# =============================================================================
# CLR TRANSFORM + FEATURE ENGINEERING
# =============================================================================

def clr_transform(X: pd.DataFrame, pseudocount: float = 1e-6) -> pd.DataFrame:
    """Centered log-ratio transform (Aitchison, compositional data)."""
    X_ps = X.values + pseudocount
    log_X = np.log(X_ps)
    geom_mean = log_X.mean(axis=1, keepdims=True)
    clr_vals = log_X - geom_mean
    return pd.DataFrame(clr_vals, index=X.index, columns=X.columns)


def build_engineered_features(
    df: pd.DataFrame, feature_classes: list
) -> pd.DataFrame:
    """
    Build an enriched feature matrix combining:
      1. Feature-class aggregates (7 cols)
      2. CLR of aggregates (7 cols)
      3. Shannon diversity per class across genomes (7 cols)
      4. Unique-genome count per class (7 cols)
      5. Pairwise ratios between classes (21 cols)
    """
    X_fc = make_feature_matrix(df, feature_classes, granularity="feature_class")
    X_gc = make_feature_matrix(df, feature_classes, granularity="genome_class")

    # 2. CLR of aggregates
    X_clr = clr_transform(X_fc)
    X_clr.columns = [f"clr_{c}" for c in X_clr.columns]

    # 3. Shannon diversity within each class
    div_data = {}
    for fc in feature_classes:
        cols = [c for c in X_gc.columns if c.endswith(f"|{fc}")]
        if not cols:
            continue
        vals = X_gc[cols].values
        totals = vals.sum(axis=1, keepdims=True)
        with np.errstate(divide="ignore", invalid="ignore"):
            p = np.where(totals > 0, vals / totals, 0.0)
            log_p = np.where(p > 0, np.log(p), 0.0)
        div_data[f"shannon_{fc}"] = -(p * log_p).sum(axis=1)
    X_div = pd.DataFrame(div_data, index=X_gc.index).reindex(X_fc.index, fill_value=0.0)

    # 4. Unique genome count per class
    count_data = {}
    for fc in feature_classes:
        cols = [c for c in X_gc.columns if c.endswith(f"|{fc}")]
        if cols:
            count_data[f"n_genomes_{fc}"] = (X_gc[cols] > 0).sum(axis=1)
    X_cnt = pd.DataFrame(count_data, index=X_gc.index).reindex(X_fc.index, fill_value=0.0)

    # 5. Pairwise log-ratios
    ratio_data = {}
    fcs = [c for c in feature_classes if c in X_fc.columns]
    for i, a in enumerate(fcs):
        for b in fcs[i + 1:]:
            ratio_data[f"lratio_{a}_{b}"] = np.log(
                (X_fc[a] + 1e-6) / (X_fc[b] + 1e-6)
            )
    X_rat = pd.DataFrame(ratio_data, index=X_fc.index)

    X_eng = pd.concat([X_fc, X_clr, X_div, X_cnt, X_rat], axis=1).fillna(0.0)
    print(f"  [feature_eng] Shape: {X_eng.shape} "
          f"(base={X_fc.shape[1]}, clr={X_clr.shape[1]}, "
          f"shannon={X_div.shape[1]}, counts={X_cnt.shape[1]}, "
          f"ratios={X_rat.shape[1]})")
    return X_eng


# =============================================================================
# PREPROCESSING / FEATURE FILTERING
# =============================================================================

def drop_constant_and_leaky(
    X: pd.DataFrame, y: pd.Series, auc_thresh: float = 0.95
) -> pd.DataFrame:
    """Remove zero-variance and univariately-leaky features (fitted on X/y)."""
    ybin = _binarize(y)
    keep, flagged = [], []
    for c in X.columns:
        col = X[c].values
        if np.allclose(col, col[0]):
            continue
        try:
            s = roc_auc_score(ybin, col)
            if s >= auc_thresh or (1.0 - s) >= auc_thresh:
                flagged.append(c)
                continue
        except Exception:
            pass
        keep.append(c)

    if not keep and flagged:
        print(
            f"[WARN] All {len(flagged)} non-constant features have univariate "
            f"AUROC ≥ {auc_thresh}. Leakage guard bypassed for this fold.",
            file=sys.stderr,
        )
        return X[flagged]

    return X[keep] if keep else X[[]]


def remove_correlated_features(X: pd.DataFrame, thresh: float = 0.95) -> pd.DataFrame:
    """Remove highly correlated features, keeping higher-variance one."""
    if X.shape[1] <= 1:
        return X
    corr = X.corr().abs()
    upper = corr.where(np.triu(np.ones(corr.shape, dtype=bool), k=1))
    to_drop = set()
    for col in upper.columns:
        partners = upper.index[upper[col] > thresh].tolist()
        for p in partners:
            if p in to_drop or col in to_drop:
                continue
            if X[p].var() <= X[col].var():
                to_drop.add(p)
            else:
                to_drop.add(col)
    return X.drop(columns=list(to_drop))


def _binarize(y: pd.Series) -> pd.Series:
    pos = _detect_pos(y)
    return (y == pos).astype(int)


def _detect_pos(y: pd.Series) -> str:
    uni = list(y.unique())
    for candidate in ["HEALTHY", "Healthy", "HFE", "hfe"]:
        if candidate in uni:
            return candidate
    return sorted(uni)[0]


def _align(X: pd.DataFrame, cols) -> pd.DataFrame:
    return X.reindex(columns=cols, fill_value=0.0)


# =============================================================================
# DATA SPLIT  80 / 10 / 10 (Optional if no hold-out is desired and also for LOSO or Nested CV)
# =============================================================================

def split_80_10_10(X: pd.DataFrame, y: pd.Series, random_state: int = 42):
    """
    Stratified 80 / 10 / 10 split.

    Both train_test_split calls use stratify= so the Healthy:PWD ratio is
    preserved in all three partitions.  The function asserts that every
    partition contains at least one sample from each class before returning.
    """
    X_tv, X_test, y_tv, y_test = train_test_split(
        X, y, test_size=0.10, stratify=y, random_state=random_state
    )
    X_train, X_val, y_train, y_val = train_test_split(
        X_tv, y_tv, test_size=round(1 / 9, 6), stratify=y_tv, random_state=random_state
    )

    # Verify both classes are present in every partition
    for name, yi in [("Train", y_train), ("Val", y_val), ("Test", y_test)]:
        missing = set(y.unique()) - set(yi.unique())
        if missing:
            raise RuntimeError(
                f"[FATAL] Split produced a '{name}' partition missing class(es): "
                f"{missing}.  Increase dataset size or check class labels."
            )

    return X_train, X_val, X_test, y_train, y_val, y_test


# =============================================================================
# MODELS & HYPERPARAMETER GRIDS
# =============================================================================

def build_models(class_ratio: float = 1.0):
    """
    Returns 4 classifiers: LR, Ridge, RF, XGB.
    class_ratio = n_neg / n_pos (used for XGB scale_pos_weight).
    All models use class-weighted training to handle imbalance.

    Ridge vs LR:
      - LR (L2 penalty) uses a probabilistic logistic loss → predict_proba.
      - Ridge minimises squared hinge-like loss → uses decision_function,
        converted to probability via sigmoid in _predict_proba().
      - Ridge is faster and better conditioned in very high-dimensional settings;
        LR is more interpretable for odds-ratio reporting.
    """
    mdls = {
        "LR": Pipeline([
            ("sc", StandardScaler()),
            ("clf", LogisticRegression(
                penalty="l2", C=1.0, solver="liblinear",
                max_iter=2000, class_weight="balanced", random_state=42,
            )),
        ]),
        "Ridge": Pipeline([
            ("sc", StandardScaler()),
            ("clf", RidgeClassifier(
                alpha=1.0, class_weight="balanced",
            )),
        ]),
        "RF": RandomForestClassifier(
            n_estimators=300, max_depth=6, min_samples_leaf=3,
            class_weight="balanced_subsample", random_state=42, n_jobs=-1,
        ),
    }
    if HAS_XGB:
        mdls["XGB"] = XGBClassifier(
            n_estimators=300, max_depth=3, learning_rate=0.05,
            subsample=0.8, colsample_bytree=0.8, reg_lambda=10.0,
            scale_pos_weight=float(class_ratio),
            random_state=42, n_jobs=-1, tree_method="hist",
            eval_metric="auc",
        )
    return mdls


PARAM_GRIDS = {
    "LR":    {"clf__C":        [0.01, 0.1, 1.0, 10.0]},
    "Ridge": {"clf__alpha":    [0.01, 0.1, 1.0, 10.0, 100.0]},
    "RF":    {"n_estimators":  [100, 300], "max_depth": [3, 5, 7]},
    "XGB":   {"n_estimators":  [100, 300], "max_depth": [2, 3], "learning_rate": [0.01, 0.05]},
}

PARAM_GRIDS_EXPANDED = {
    "LR":    {"clf__C":        [0.001, 0.01, 0.1, 1.0, 10.0, 100.0]},
    "Ridge": {"clf__alpha":    [0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0]},
    "RF":    {"n_estimators":  [100, 200, 400], "max_depth": [3, 5, 7, None],
              "min_samples_leaf": [1, 3, 5]},
    "XGB":   {"n_estimators":  [100, 200, 400], "max_depth": [2, 3, 5],
              "learning_rate": [0.005, 0.01, 0.05, 0.1],
              "subsample": [0.7, 0.9], "colsample_bytree": [0.7, 0.9]},
}


# =============================================================================
# ROC / CONFUSION MATRIX HELPERS
# =============================================================================

def plot_roc(fpr, tpr, auc_val, title, outpath):
    plt.figure(figsize=(5, 5))
    plt.plot([0, 1], [0, 1], "k--", lw=1)
    plt.plot(fpr, tpr, lw=2, label=f"AUROC = {auc_val:.3f}")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(title)
    plt.legend(loc="lower right", fontsize=9)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.savefig(str(outpath).replace(".png", ".tiff"), dpi=300)
    plt.close()


def plot_pr_curve(y_true, y_prob, title, outpath):
    prec, rec, _ = precision_recall_curve(y_true, y_prob)
    ap = average_precision_score(y_true, y_prob)
    plt.figure(figsize=(5, 5))
    plt.plot(rec, prec, lw=2, label=f"AP = {ap:.3f}")
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title(title)
    plt.legend(loc="upper right", fontsize=9)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()


def plot_confusion_matrix(y_true, y_pred, title, outpath, labels=None):
    cm = confusion_matrix(y_true, y_pred)
    if labels is None:
        labels = ["Neg (PWD)", "Pos (Healthy)"]
    fig, ax = plt.subplots(figsize=(4, 3))
    if HAS_SNS:
        sns.heatmap(
            cm, annot=True, fmt="d", cmap="Blues",
            xticklabels=labels, yticklabels=labels, ax=ax,
        )
    else:
        im = ax.imshow(cm, cmap="Blues")
        for i in range(cm.shape[0]):
            for j in range(cm.shape[1]):
                ax.text(j, i, str(cm[i, j]), ha="center", va="center")
        ax.set_xticks([0, 1]); ax.set_xticklabels(labels)
        ax.set_yticks([0, 1]); ax.set_yticklabels(labels)
    ax.set_xlabel("Predicted"); ax.set_ylabel("Actual")
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()
    return cm


def mean_roc_from_curves(curves):
    mean_fpr = np.linspace(0, 1, 200)
    tprs, aucs = [], []
    for fpr, tpr, auc in curves:
        if not np.isnan(auc):
            tprs.append(np.interp(mean_fpr, fpr, tpr))
            tprs[-1][0] = 0.0
            aucs.append(auc)
    if not tprs:
        return mean_fpr, np.zeros(200), np.zeros(200), np.nan, np.nan
    m_tpr = np.mean(tprs, axis=0)
    s_tpr = np.std(tprs, axis=0)
    return mean_fpr, m_tpr, s_tpr, float(np.mean(aucs)), float(np.std(aucs))


def plot_mean_roc(mean_fpr, mean_tpr, std_tpr, mean_auc, std_auc, title, outpath):
    plt.figure(figsize=(5, 5))
    plt.plot([0, 1], [0, 1], "k--", lw=1)
    plt.plot(mean_fpr, mean_tpr, lw=2,
             label=f"Mean AUROC = {mean_auc:.3f} ± {std_auc:.3f}")
    plt.fill_between(
        mean_fpr,
        np.maximum(mean_tpr - std_tpr, 0),
        np.minimum(mean_tpr + std_tpr, 1),
        alpha=0.2, label="±1 SD",
    )
    plt.xlabel("FPR"); plt.ylabel("TPR")
    plt.title(title); plt.legend(loc="lower right", fontsize=9)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.savefig(str(outpath).replace(".png", ".tiff"), dpi=300)
    plt.close()


# =============================================================================
# NESTED CROSS-VALIDATION
# =============================================================================

def run_nested_cv(
    X: pd.DataFrame, y: pd.Series, models: dict, outdir: Path, tag: str,
    n_outer: int = 5, n_inner: int = 3,
    leak_auc: float = 0.95, topk: int = 50,
    groups: pd.Series = None,
    param_grids: dict = None,
):
    outdir.mkdir(parents=True, exist_ok=True)
    ybin = _binarize(y)

    if param_grids is None:
        param_grids = PARAM_GRIDS

    if groups is not None:
        from sklearn.model_selection import StratifiedGroupKFold
        grp_arr = groups.reindex(X.index).fillna("unknown").values
        outer_cv = StratifiedGroupKFold(n_splits=n_outer)
        split_kwargs = {"groups": grp_arr}
    else:
        outer_cv = StratifiedKFold(n_splits=n_outer, shuffle=True, random_state=42)
        grp_arr = None
        split_kwargs = {}
    inner_cv = StratifiedKFold(n_splits=n_inner, shuffle=True, random_state=42)

    results = {}
    all_pred_rows, all_fold_rows = [], []

    for mdl_name, mdl in models.items():
        oof_probs = np.zeros(len(ybin))
        curves, fold_aucs = [], []
        fold_counter = 0

        for tr_idx, te_idx in outer_cv.split(X, ybin, **split_kwargs):
            fold_counter += 1
            Xtr, Xte = X.iloc[tr_idx].copy(), X.iloc[te_idx].copy()
            ytr, yte = ybin.iloc[tr_idx], ybin.iloc[te_idx]

            if ytr.nunique() < 2 or yte.nunique() < 2:
                continue

            Xtr_f = drop_constant_and_leaky(Xtr, ytr, auc_thresh=leak_auc)
            Xtr_f = remove_correlated_features(Xtr_f)
            if Xtr_f.shape[1] == 0:
                continue
            if topk > 0 and Xtr_f.shape[1] > topk:
                keep = Xtr_f.var().sort_values(ascending=False).head(topk).index
                Xtr_f = Xtr_f[keep]
            Xte_f = _align(Xte, Xtr_f.columns)

            if mdl_name in param_grids:
                gs = GridSearchCV(
                    clone(mdl), param_grids[mdl_name],
                    cv=inner_cv, scoring="roc_auc", n_jobs=-1, refit=True,
                )
                gs.fit(Xtr_f, ytr)
                best_model = gs.best_estimator_
            else:
                best_model = clone(mdl).fit(Xtr_f, ytr)

            p = _predict_proba(best_model, Xte_f)
            try:
                auc = roc_auc_score(yte, p)
            except Exception:
                auc = np.nan

            fold_aucs.append(auc)
            oof_probs[te_idx] = p
            fpr, tpr, _ = roc_curve(yte, p)
            curves.append((fpr, tpr, auc))

            te_ids = X.index[te_idx]
            all_pred_rows.append(pd.DataFrame({
                "Tag": tag, "Model": mdl_name, "Fold": fold_counter,
                "Sample_ID": te_ids, "y_true": yte.values, "y_prob": p,
            }))
            all_fold_rows.append({
                "Tag": tag, "Model": mdl_name, "Fold": fold_counter, "AUROC": auc,
            })

        try:
            oof_auc = roc_auc_score(ybin, oof_probs)
        except Exception:
            oof_auc = np.nan

        m_fpr, m_tpr, s_tpr, m_auc, s_auc = mean_roc_from_curves(curves)
        plot_mean_roc(
            m_fpr, m_tpr, s_tpr, m_auc, s_auc,
            title=f"{tag} – {mdl_name} Nested CV",
            outpath=outdir / f"roc_nestedCV_{tag}_{mdl_name}.png",
        )

        y_pred_bin = (oof_probs >= 0.5).astype(int)
        cm = plot_confusion_matrix(
            ybin.values, y_pred_bin,
            title=f"CM OOF – {tag} {mdl_name}",
            outpath=outdir / f"cm_nestedCV_{tag}_{mdl_name}.png",
        )

        results[mdl_name] = {
            "fold_aucs": fold_aucs,
            "mean_auroc": float(np.nanmean(fold_aucs)),
            "std_auroc":  float(np.nanstd(fold_aucs)),
            "oof_auroc":  float(oof_auc) if not np.isnan(oof_auc) else None,
            "oof_probs":  oof_probs,
            "cm_oof":     cm.tolist() if cm is not None else None,
        }
        oof_str = f"{oof_auc:.3f}" if not np.isnan(oof_auc) else "nan"
        print(
            f"  [{tag}] {mdl_name}  nested CV: {m_auc:.3f}±{s_auc:.3f} "
            f"(OOF={oof_str})"
        )

    if all_pred_rows:
        pd.concat(all_pred_rows, ignore_index=True).to_csv(
            outdir / f"predictions_nestedCV_{tag}.csv", index=False
        )
    if all_fold_rows:
        pd.DataFrame(all_fold_rows).to_csv(
            outdir / f"fold_aucs_nestedCV_{tag}.csv", index=False
        )

    summary = pd.DataFrame([{
        "Tag": tag, "Model": k, "CV_Method": f"Nested_{n_outer}x{n_inner}fold",
        "Mean_AUROC": v["mean_auroc"], "SD_AUROC": v["std_auroc"],
        "OOF_AUROC": v["oof_auroc"],
    } for k, v in results.items()])
    summary.to_csv(outdir / f"summary_nestedCV_{tag}.csv", index=False)

    return results


# =============================================================================
# HOLDOUT EVALUATION  (train 80% → tune on val 10% → report on test 10%)
# =============================================================================

def evaluate_holdout(
    X_train, y_train, X_val, y_val, X_test, y_test,
    models: dict, outdir: Path, tag: str,
    leak_auc: float = 0.95, topk: int = 50,
):
    outdir.mkdir(parents=True, exist_ok=True)

    pos = _detect_pos(y_train)
    ytr_b = (y_train == pos).astype(int)
    yva_b = (y_val   == pos).astype(int)
    yte_b = (y_test  == pos).astype(int)

    Xtr_f = drop_constant_and_leaky(X_train, ytr_b, auc_thresh=leak_auc)
    Xtr_f = remove_correlated_features(Xtr_f)
    if topk > 0 and Xtr_f.shape[1] > topk:
        keep = Xtr_f.var().sort_values(ascending=False).head(topk).index
        Xtr_f = Xtr_f[keep]

    cols = Xtr_f.columns
    Xva_f = _align(X_val,  cols)
    Xte_f = _align(X_test, cols)

    results = {}
    rows = []

    for mdl_name, mdl in models.items():
        best_model, best_val_auc, best_params = None, -1.0, {}

        if mdl_name in PARAM_GRIDS:
            keys = list(PARAM_GRIDS[mdl_name].keys())
            vals = list(PARAM_GRIDS[mdl_name].values())
            for combo in itertools.product(*vals):
                params = dict(zip(keys, combo))
                m = clone(mdl)
                m.set_params(**params)
                m.fit(Xtr_f, ytr_b)
                p_val = _predict_proba(m, Xva_f)
                try:
                    val_auc = roc_auc_score(yva_b, p_val)
                except Exception:
                    val_auc = 0.0
                if val_auc > best_val_auc:
                    best_val_auc = val_auc
                    best_model = m
                    best_params = params
        else:
            best_model = clone(mdl).fit(Xtr_f, ytr_b)
            p_val = _predict_proba(best_model, Xva_f)
            try:
                best_val_auc = roc_auc_score(yva_b, p_val)
            except Exception:
                best_val_auc = float("nan")

        p_test = _predict_proba(best_model, Xte_f)
        y_pred = (p_test >= 0.5).astype(int)
        try:
            test_auc = roc_auc_score(yte_b, p_test)
        except Exception:
            test_auc = float("nan")

        fpr, tpr, _ = roc_curve(yte_b, p_test)
        plot_roc(fpr, tpr, test_auc,
                 title=f"Holdout Test – {tag} {mdl_name}",
                 outpath=outdir / f"roc_holdout_{tag}_{mdl_name}.png")
        plot_pr_curve(yte_b, p_test,
                      title=f"PR Curve – {tag} {mdl_name}",
                      outpath=outdir / f"pr_holdout_{tag}_{mdl_name}.png")
        cm = plot_confusion_matrix(
            yte_b.values, y_pred,
            title=f"CM Holdout Test – {tag} {mdl_name}",
            outpath=outdir / f"cm_holdout_{tag}_{mdl_name}.png",
        )

        print(
            f"  [{tag}] {mdl_name}  val={best_val_auc:.3f}  "
            f"test={test_auc:.3f}  best_params={best_params}"
        )

        results[mdl_name] = {
            "val_auroc": float(best_val_auc),
            "test_auroc": float(test_auc),
            "best_params": best_params,
            "best_model": best_model,
            "X_cols": cols,
            "cm_test": cm.tolist() if cm is not None else None,
        }
        rows.append({
            "Tag": tag, "Model": mdl_name, "CV_Method": "Holdout_10pct_test",
            "Val_AUROC": float(best_val_auc), "Test_AUROC": float(test_auc),
            "Best_Params": str(best_params),
        })

    pd.DataFrame(rows).to_csv(outdir / f"summary_holdout_{tag}.csv", index=False)

    pred_rows_all = []
    for mdl_name, res in results.items():
        best_model = res["best_model"]
        cols_      = res["X_cols"]
        Xte_f_     = _align(X_test, cols_)
        p_test_    = _predict_proba(best_model, Xte_f_)
        pos_       = _detect_pos(y_train)
        yte_b_     = (y_test == pos_).astype(int)
        pred_rows_all.append(pd.DataFrame({
            "Model":     mdl_name,
            "Sample_ID": X_test.index,
            "y_true":    yte_b_.values,
            "y_prob":    p_test_,
        }))
    if pred_rows_all:
        pd.concat(pred_rows_all, ignore_index=True).to_csv(
            outdir / f"predictions_holdout_{tag}.csv", index=False
        )

    return results


# =============================================================================
# LOSO CROSS-VALIDATION
# =============================================================================

def run_loso(
    X: pd.DataFrame, y: pd.Series, study: pd.Series,
    models: dict, outdir: Path, tag: str,
    leak_auc: float = 0.95, topk: int = 50,
):
    outdir.mkdir(parents=True, exist_ok=True)
    studies = sorted(study.unique())
    if len(studies) < 2:
        print(f"[WARN] LOSO requires ≥2 studies; found {studies}. Skipping.", file=sys.stderr)
        return None

    ybin = _binarize(y)
    results = {}

    for mdl_name, mdl in models.items():
        study_rows, pred_rows = [], []
        oof_probs = np.zeros(len(ybin))

        for held_out in studies:
            te_mask = (study == held_out).values
            tr_mask = ~te_mask

            Xtr, Xte = X[tr_mask].copy(), X[te_mask].copy()
            ytr, yte = ybin[tr_mask], ybin[te_mask]

            if ytr.nunique() < 2 or yte.nunique() < 2:
                print(
                    f"[WARN] LOSO study '{held_out}' missing a class. Skipping fold.",
                    file=sys.stderr,
                )
                continue

            Xtr_f = drop_constant_and_leaky(Xtr, ytr, auc_thresh=leak_auc)
            Xtr_f = remove_correlated_features(Xtr_f)
            if topk > 0 and Xtr_f.shape[1] > topk:
                keep = Xtr_f.var().sort_values(ascending=False).head(topk).index
                Xtr_f = Xtr_f[keep]
            Xte_f = _align(Xte, Xtr_f.columns)

            m = clone(mdl).fit(Xtr_f, ytr)
            p = _predict_proba(m, Xte_f)
            try:
                auc = roc_auc_score(yte, p)
            except Exception:
                auc = float("nan")

            oof_probs[te_mask] = p
            study_rows.append({
                "Study": held_out, "AUROC": auc,
                "N_train": int(tr_mask.sum()), "N_test": int(te_mask.sum()),
                "Pos_in_test": int(yte.sum()),
            })
            pred_rows.append(pd.DataFrame({
                "Model": mdl_name, "Held_Out_Study": held_out,
                "Sample_ID": X.index[te_mask],
                "y_true": yte.values, "y_prob": p,
            }))
            print(f"  LOSO [{tag}/{mdl_name}] held-out={held_out}  AUROC={auc:.3f}")

        df_studies = pd.DataFrame(study_rows)
        df_studies.to_csv(outdir / f"loso_per_study_{tag}_{mdl_name}.csv", index=False)
        if pred_rows:
            pd.concat(pred_rows, ignore_index=True).to_csv(
                outdir / f"loso_preds_{tag}_{mdl_name}.csv", index=False
            )

        try:
            oof_auc = roc_auc_score(ybin, oof_probs)
        except Exception:
            oof_auc = float("nan")

        y_pred_bin = (oof_probs >= 0.5).astype(int)
        plot_confusion_matrix(
            ybin.values, y_pred_bin,
            title=f"LOSO CM – {tag} {mdl_name}",
            outpath=outdir / f"cm_loso_{tag}_{mdl_name}.png",
        )

        results[mdl_name] = {
            "per_study": df_studies,
            "mean_auroc": float(df_studies["AUROC"].mean()),
            "std_auroc":  float(df_studies["AUROC"].std()),
            "oof_auroc":  float(oof_auc),
        }

    summary = pd.DataFrame([{
        "Tag": tag, "Model": k, "CV_Method": "LOSO",
        "Mean_AUROC": v["mean_auroc"], "SD_AUROC": v["std_auroc"],
        "OOF_AUROC": v["oof_auroc"],
    } for k, v in results.items()])
    summary.to_csv(outdir / f"summary_loso_{tag}.csv", index=False)
    return results


# =============================================================================
# FEATURE CLASS COMBINATION SEARCH
# =============================================================================

def run_combo_search(
    df: pd.DataFrame, available_classes: list, y_all: pd.Series,
    models: dict, outdir: Path, tag: str,
    leak_auc: float = 0.95, topk: int = 50,
    max_combo_size: int = None,
    granularity: str = "feature_class",
):
    outdir.mkdir(parents=True, exist_ok=True)
    ybin_all = _binarize(y_all)
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    n = len(available_classes)
    max_r = max_combo_size if max_combo_size else n
    total = sum(
        len(list(itertools.combinations(available_classes, r)))
        for r in range(1, max_r + 1)
    )
    print(
        f"\n[INFO] Combo search: {total} combos × {len(models)} models × 5 folds ...",
        flush=True,
    )

    rows = []
    done = 0

    for r in range(1, max_r + 1):
        for combo in itertools.combinations(available_classes, r):
            combo_name = "+".join(combo)

            X_c = make_feature_matrix(df, feature_classes=list(combo),
                                       granularity=granularity)
            if X_c.empty:
                done += 1
                continue

            common = [i for i in y_all.index if i in X_c.index]
            if len(common) < 10:
                done += 1
                continue

            X_c = X_c.loc[common]
            y_c = ybin_all.loc[common]

            for mdl_name, mdl in models.items():
                fold_aucs, oof_probs = [], np.zeros(len(y_c))

                for tr_idx, te_idx in skf.split(X_c, y_c):
                    Xtr, Xte = X_c.iloc[tr_idx].copy(), X_c.iloc[te_idx].copy()
                    ytr, yte = y_c.iloc[tr_idx], y_c.iloc[te_idx]
                    if ytr.nunique() < 2 or yte.nunique() < 2:
                        continue

                    Xtr_f = drop_constant_and_leaky(Xtr, ytr, auc_thresh=leak_auc)
                    if Xtr_f.shape[1] == 0:
                        continue
                    if topk > 0 and Xtr_f.shape[1] > topk:
                        keep = Xtr_f.var().sort_values(ascending=False).head(topk).index
                        Xtr_f = Xtr_f[keep]
                    Xte_f = _align(Xte, Xtr_f.columns)

                    m = clone(mdl).fit(Xtr_f, ytr)
                    p = _predict_proba(m, Xte_f)
                    try:
                        fold_aucs.append(roc_auc_score(yte, p))
                    except Exception:
                        fold_aucs.append(float("nan"))
                    oof_probs[te_idx] = p

                if not fold_aucs:
                    continue

                try:
                    oof_auc = roc_auc_score(y_c, oof_probs)
                except Exception:
                    oof_auc = float("nan")

                cm_raw = confusion_matrix(y_c, (oof_probs >= 0.5).astype(int))
                tn, fp, fn, tp = (
                    (int(cm_raw[0, 0]), int(cm_raw[0, 1]),
                     int(cm_raw[1, 0]), int(cm_raw[1, 1]))
                    if cm_raw.shape == (2, 2)
                    else (None, None, None, None)
                )
                rows.append({
                    "Combo": combo_name, "N_Classes": r, "Model": mdl_name,
                    "N_Samples": len(y_c),
                    "Mean_AUROC": float(np.nanmean(fold_aucs)),
                    "SD_AUROC":   float(np.nanstd(fold_aucs)),
                    "OOF_AUROC":  float(oof_auc),
                    "TN": tn, "FP": fp, "FN": fn, "TP": tp,
                })

            done += 1
            if done % 20 == 0:
                print(f"  {done}/{total} combos done...", flush=True)

    df_combos = pd.DataFrame(rows)
    if df_combos.empty:
        print("[WARN] No combo results generated.", file=sys.stderr)
        return df_combos

    df_combos = df_combos.sort_values("Mean_AUROC", ascending=False)
    df_combos.to_csv(outdir / f"combo_search_{tag}.csv", index=False)

    top10 = (
        df_combos.groupby("Combo")["Mean_AUROC"]
        .max()
        .sort_values(ascending=False)
        .head(10)
    )
    plt.figure(figsize=(11, 5))
    top10[::-1].plot(kind="barh", color="steelblue")
    plt.xlabel("Mean AUROC (5-fold CV)")
    plt.title(f"Top-10 Feature Class Combinations – {tag}")
    plt.tight_layout()
    plt.savefig(outdir / f"top10_combos_{tag}.png", dpi=150)
    plt.close()

    print(f"\n  Top-5 combinations for {tag}:")
    for combo, auc in top10.head(5).items():
        print(f"    {combo:55s}  AUROC={auc:.3f}")

    return df_combos


# =============================================================================
# UTILITY
# =============================================================================

def _predict_proba(model, X: pd.DataFrame) -> np.ndarray:
    """
    Unified probability extractor.
    - Pipeline/LogisticRegression/XGB: uses predict_proba[:, 1]
    - RidgeClassifier: no predict_proba → sigmoid of decision_function
    """
    if hasattr(model, "predict_proba"):
        return model.predict_proba(X)[:, 1]
    d = model.decision_function(X)
    return 1.0 / (1.0 + np.exp(-d))

# =============================================================================
# MAIN
# =============================================================================

def main():
    ap = argparse.ArgumentParser(
        description="SGI ML Pipeline v6 — LR + Ridge + RF + XGB; verified stratified split."
    )
    ap.add_argument("--genomics", required=True,
                    help="Long-format TSV: Sample_ID, Group, Feature_Class, Genome_ID, RPKM, CPM.")
    ap.add_argument("--outdir", required=True, help="Output directory.")
    ap.add_argument("--abundance_type", default="RPKM", choices=["RPKM", "CPM"],
                    help="Abundance column to use (default: RPKM).")
    ap.add_argument("--feature_classes",
                    default="AMP,AMR,BGC,CAZyme,VFDB,gutSMASH,SGI",
                    help="Comma-separated feature classes to include.")
    ap.add_argument("--topk", type=int, default=50,
                    help="Top-K features by variance per fold (0 = no cap).")
    ap.add_argument("--leak_auc", type=float, default=0.95,
                    help="Univariate AUROC leak-guard threshold.")
    ap.add_argument("--study_col", default=None,
                    help="Column in the input TSV for study/batch (LOSO).")
    ap.add_argument("--study_meta", default=None,
                    help="External TSV with Sample_ID,Study for LOSO.")
    ap.add_argument("--max_combo_size", type=int, default=None,
                    help="Max number of feature classes per combination.")
    ap.add_argument("--granularity", default="feature_class",
                    choices=["feature_class", "genome_class"],
                    help=(
                        "feature_class: sum per Feature_Class → 103×7 (RPKM default). "
                        "genome_class: per Genome_ID|Feature_Class → 103×1774 (required for CPM)."
                    ))
    ap.add_argument("--skip_combo_search", action="store_true",
                    help="Skip the combination search.")
    ap.add_argument("--expand_grid", action="store_true",
                    help="Use expanded hyperparameter grid for nested CV.")
    ap.add_argument("--farm_stratified", action="store_true",
                    help="Use StratifiedGroupKFold (farm-stratified) for nested CV.")
    ap.add_argument("--clr", action="store_true",
                    help="Apply centered log-ratio transform to the feature matrix.")
    ap.add_argument("--feature_engineering", action="store_true",
                    help="Build enriched feature matrix (CLR + Shannon + counts + ratios).")
    ap.add_argument("--random_state", type=int, default=42)
    ap.add_argument("--no_holdout", action="store_true",
                    help="Skip the 80/10/10 holdout split; use all samples for nested CV and feature importance.")
    args = ap.parse_args()

    out_root = Path(args.outdir)
    out_root.mkdir(parents=True, exist_ok=True)

    feature_classes = [x.strip() for x in args.feature_classes.split(",") if x.strip()]

    print(f"\n[INFO] Loading data:  {args.genomics}")
    print(f"[INFO] Abundance:     {args.abundance_type}")
    G = read_long_format(args.genomics, abundance_type=args.abundance_type)

    if "Group" not in G.columns:
        raise SystemExit("[FATAL] 'Group' column not found in input.")

    y_all_raw = get_labels(G)
    groups = sorted(y_all_raw.unique())
    print(f"[INFO] Groups in data: {groups}")

    pos_class = neg_class = None
    for p, n in [("Healthy", "PWD"), ("HEALTHY", "PWD"), ("HFE", "LFE")]:
        if p in groups and n in groups:
            pos_class, neg_class = p, n
            break
    if pos_class is None:
        pos_class, neg_class = groups[0], groups[1]
        print(f"[WARN] Auto-detected: pos={pos_class}, neg={neg_class}", file=sys.stderr)

    y_all = y_all_raw[y_all_raw.isin([pos_class, neg_class])]
    G_f = G[G["Sample_ID"].isin(y_all.index)].copy()

    avail = [c for c in feature_classes if c in G_f["Feature_Class"].unique()]
    print(f"[INFO] Feature classes available: {avail}")

    n_pos   = int((y_all == pos_class).sum())
    n_neg   = int((y_all == neg_class).sum())
    n_total = len(y_all)
    ratio   = round(n_neg / n_pos, 3) if n_pos else float("nan")
    print(f"\n[INFO] ===== SAMPLE SIZES =====")
    print(f"  {pos_class} (positive): {n_pos}")
    print(f"  {neg_class} (negative): {n_neg}")
    print(f"  Total              : {n_total}")
    print(f"  Imbalance ratio    : {ratio:.3f}  (neg/pos)")

    if args.feature_engineering:
        print(f"\n[INFO] Building ENGINEERED feature matrix ...")
        X_all = build_engineered_features(G_f, avail)
        common_ids = [i for i in y_all.index if i in X_all.index]
        X_all = X_all.loc[common_ids]
    else:
        print(f"\n[INFO] Building feature matrix (granularity={args.granularity})...")
        X_all = make_feature_matrix(G_f, feature_classes=avail, granularity=args.granularity)
        common_ids = [i for i in y_all.index if i in X_all.index]
        X_all = X_all.loc[common_ids]
        if args.clr:
            print("[INFO] Applying CLR transform...")
            X_all = clr_transform(X_all)

    y_all = y_all.loc[common_ids]
    feat_desc = "engineered features" if args.feature_engineering else (
        "feature classes" if args.granularity == "feature_class" else "Genome|Class features"
    )
    print(f"  Shape: {X_all.shape}  (samples × {feat_desc})")

    def _split_report(name, y_part, pos_c, neg_c):
        n_p = int((y_part == pos_c).sum())
        n_n = int((y_part == neg_c).sum())
        r   = round(n_n / n_p, 3) if n_p else float("nan")
        print(f"  {name:<8}: n={len(y_part):3d}  "
              f"{pos_c}={n_p}  {neg_c}={n_n}  ratio={r:.3f}")

    if args.no_holdout:
        print(f"\n[INFO] --no_holdout: using all {n_total} samples for nested CV and feature importance.")
        X_train, y_train = X_all, y_all
        X_val,   y_val   = X_all.iloc[:0], y_all.iloc[:0]   # empty
        X_test,  y_test  = X_all.iloc[:0], y_all.iloc[:0]   # empty
        _split_report("All (train)", y_train, pos_class, neg_class)
    else:
        print(f"\n[INFO] Stratified 80/10/10 split ...")
        X_train, X_val, X_test, y_train, y_val, y_test = split_80_10_10(
            X_all, y_all, random_state=args.random_state
        )
        print(f"  {'Partition':<8}  {'n':>4}  class breakdown  (neg/pos ratio)")
        _split_report("Train",  y_train, pos_class, neg_class)
        _split_report("Val",    y_val,   pos_class, neg_class)
        _split_report("Test",   y_test,  pos_class, neg_class)
        print(f"  [OK] Both classes present in all three partitions.")

    models = build_models(class_ratio=ratio)
    print(f"[INFO] Models: {list(models.keys())}")

    study_labels = None
    study_col_used = "not provided"
    n_studies = "N/A"

    if args.study_col and args.study_col in G.columns:
        sm = (
            G[["Sample_ID", args.study_col]]
            .drop_duplicates()
            .set_index("Sample_ID")[args.study_col]
        )
        study_labels = sm.reindex(y_all.index).dropna()
        study_col_used = args.study_col
        n_studies = int(study_labels.nunique())
        print(f"[INFO] Study column: '{args.study_col}'  → {n_studies} studies")

    elif args.study_meta:
        sm = pd.read_csv(args.study_meta, sep=None, engine="python")
        sm.columns = [c.strip() for c in sm.columns]
        sm = sm.rename(columns={sm.columns[0]: "Sample_ID", sm.columns[1]: "Study"})
        sm["Sample_ID"] = sm["Sample_ID"].astype(str).str.strip()
        sm = sm.set_index("Sample_ID")["Study"]
        study_labels = sm.reindex(y_all.index).dropna()
        study_col_used = args.study_meta
        n_studies = int(study_labels.nunique())
        print(f"[INFO] Study metadata: {args.study_meta}  → {n_studies} studies")
    else:
        print("[INFO] No study/batch info — LOSO will be skipped.")

    stats = {
        "dataset": {
            "abundance_type": args.abundance_type,
            "n_total": n_total, "n_pos": n_pos, "n_neg": n_neg,
            "class_ratio": ratio, "scale_pos_weight": ratio,
            "feature_classes": avail,
            "study_col": study_col_used, "n_studies": n_studies,
            "n_train": X_train.shape[0], "n_val": X_val.shape[0],
            "n_test": X_test.shape[0],
            "leak_auc": args.leak_auc, "topk": args.topk,
        },
        "nested_cv": {}, "holdout": {}, "loso": {}, "combo_search": None,
    }
    task_tag = f"{pos_class}_vs_{neg_class}"

    active_grids = PARAM_GRIDS_EXPANDED if args.expand_grid else PARAM_GRIDS
    farm_groups  = study_labels if args.farm_stratified and study_labels is not None else None

    if args.expand_grid:
        print("[INFO] Using EXPANDED hyperparameter grid.")
    if args.farm_stratified:
        if farm_groups is not None:
            print(f"[INFO] Farm-stratified nested CV ({farm_groups.nunique()} farms).")
        else:
            print("[WARN] --farm_stratified requested but no study metadata — using standard CV.")

    print(f"\n[INFO] ===== NESTED 5×3-fold CV =====")
    ncv_dir = out_root / "nested_cv"
    ncv_X = X_all if farm_groups is not None else X_train
    ncv_y = y_all if farm_groups is not None else y_train
    nested_res = run_nested_cv(
        ncv_X, ncv_y, models, ncv_dir, tag=task_tag,
        n_outer=5, n_inner=3, leak_auc=args.leak_auc, topk=args.topk,
        groups=farm_groups, param_grids=active_grids,
    )
    stats["nested_cv"][task_tag] = {
        k: {kk: vv for kk, vv in v.items() if kk != "oof_probs"}
        for k, v in nested_res.items()
    }

    if not args.no_holdout:
        print(f"\n[INFO] ===== HOLDOUT EVALUATION (80/10/10) =====")
        ho_dir = out_root / "holdout"
        holdout_res = evaluate_holdout(
            X_train, y_train, X_val, y_val, X_test, y_test,
            models, ho_dir, tag=task_tag,
            leak_auc=args.leak_auc, topk=args.topk,
        )
        stats["holdout"][task_tag] = {
            k: {kk: vv for kk, vv in v.items() if kk not in ("best_model", "X_cols")}
            for k, v in holdout_res.items()
        }
    else:
        print(f"\n[INFO] ===== HOLDOUT EVALUATION skipped (--no_holdout) =====")

    print(f"\n[INFO] ===== FEATURE IMPORTANCE =====")
    imp_dir = out_root / "feature_importance"
    imp_dir.mkdir(parents=True, exist_ok=True)

    ytr_bin = (y_train == pos_class).astype(int)
    Xtr_clean = drop_constant_and_leaky(X_train, ytr_bin, auc_thresh=args.leak_auc)
    Xtr_clean = remove_correlated_features(Xtr_clean)
    if args.topk > 0 and Xtr_clean.shape[1] > args.topk:
        keep = Xtr_clean.var().sort_values(ascending=False).head(args.topk).index
        Xtr_clean = Xtr_clean[keep]


    # Save directory for serialised models
    models_dir = out_root / "saved_models"
    models_dir.mkdir(parents=True, exist_ok=True)

    for mdl_name, mdl in models.items():
        print(f"  {mdl_name} importances...")
        final_model = clone(mdl).fit(Xtr_clean, ytr_bin)

        # ── Save model + feature columns for external prediction ──
        model_bundle = {
            "model":         final_model,
            "feature_cols":  list(Xtr_clean.columns),
            "task_tag":      task_tag,
            "model_name":    mdl_name,
            "abundance_type": args.abundance_type,
            "granularity":   args.granularity,
            "pos_class":     _detect_pos(y_train),
        }
        joblib.dump(
            model_bundle,
            models_dir / f"model_{task_tag}_{mdl_name}.joblib",
        )
        print(f"    Saved: saved_models/model_{task_tag}_{mdl_name}.joblib")


    if study_labels is not None and study_labels.nunique() >= 2:
        print(f"\n[INFO] ===== LOSO CV =====")
        loso_dir = out_root / "loso"
        loso_ids = [i for i in y_all.index if i in study_labels.index]
        loso_res = run_loso(
            X_all.loc[loso_ids], y_all.loc[loso_ids],
            study_labels.loc[loso_ids],
            models, loso_dir, tag=task_tag,
            leak_auc=args.leak_auc, topk=args.topk,
        )
        if loso_res:
            stats["loso"][task_tag] = {
                k: {kk: vv for kk, vv in v.items() if kk != "per_study"}
                for k, v in loso_res.items()
            }
    else:
        print(f"\n[INFO] Skipping LOSO (no study metadata).")

    if not args.skip_combo_search:
        cs_dir = out_root / "combo_search"
        df_combos = run_combo_search(
            G_f, avail, y_all, models, cs_dir, tag=task_tag,
            leak_auc=args.leak_auc, topk=args.topk,
            max_combo_size=args.max_combo_size,
            granularity=args.granularity,
        )
        stats["combo_search"] = df_combos
    else:
        print("[INFO] Skipping combo search (--skip_combo_search).")

    print(f"\n[INFO] Writing master summary...")
    summary_rows = []
    for task_k, res in stats["nested_cv"].items():
        for mdl, r in res.items():
            summary_rows.append({
                "Task": task_k, "Model": mdl, "CV_Method": "Nested_5x3fold",
                "Mean_AUROC": r.get("mean_auroc"), "SD_AUROC": r.get("std_auroc"),
                "OOF_AUROC": r.get("oof_auroc"),
            })
    for task_k, res in stats["holdout"].items():
        for mdl, r in res.items():
            summary_rows.append({
                "Task": task_k, "Model": mdl, "CV_Method": "Holdout_10pct_test",
                "Mean_AUROC": r.get("test_auroc"), "SD_AUROC": None,
                "OOF_AUROC": None,
            })
    for task_k, res in stats["loso"].items():
        for mdl, r in res.items():
            summary_rows.append({
                "Task": task_k, "Model": mdl, "CV_Method": "LOSO",
                "Mean_AUROC": r.get("mean_auroc"), "SD_AUROC": r.get("std_auroc"),
                "OOF_AUROC": r.get("oof_auroc"),
            })
    pd.DataFrame(summary_rows).to_csv(out_root / "MASTER_AUROC_SUMMARY.csv", index=False)


    manifest = []
    for f in sorted(out_root.rglob("*")):
        if f.is_file() and f.suffix in (".csv", ".txt", ".png", ".tiff"):
            manifest.append({"File": str(f.relative_to(out_root)), "Type": f.suffix.lstrip(".")})
    pd.DataFrame(manifest).to_csv(out_root / "OUTPUT_MANIFEST.csv", index=False)

    print(f"\n[INFO] ===== PIPELINE COMPLETE =====")
    print(f"  Output root:      {out_root}")
    print(f"  Master summary:   {out_root / 'MASTER_AUROC_SUMMARY.csv'}")
    print(f"  Rebuttal report:  {out_root / 'REBUTTAL_REPORT.txt'}")
    print(f"  Output manifest:  {out_root / 'OUTPUT_MANIFEST.csv'}")
    print(f"  Total files:      {len(manifest)}")


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    main()
