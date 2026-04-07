#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Train ensemble models for general, TF-only, TF-tuned, and TF+transformer.

This script is training-only. It does not run evaluation, testing, or plotting.
"""

from __future__ import annotations

import argparse
import gc
import json
import logging
import os
import pickle
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd
import xgboost as xgb

from config import ENHANCERS_BED, EMBEDDINGS_CSV, TRAINING_SETS_PK, MODEL_DIR


LOGGER = logging.getLogger("ensemble_training")

MODEL_ROOT = f"{MODEL_DIR}/ensemble"
GENERAL_DIR = f"{MODEL_ROOT}/general"
TFONLY_DIR = f"{MODEL_ROOT}/tf_only"
TFTUNED_DIR = f"{MODEL_ROOT}/tf_tuned"
TFTRANS_DIR = f"{MODEL_ROOT}/tf_transformer"
GENERAL_STD2_TSV = f"{MODEL_ROOT}/general_std2_stats.tsv"
TF_STD2_TSV = f"{MODEL_ROOT}/tf_families_std2_stats.tsv"

EMB_COLS = [f"emb_{i}" for i in range(1024)]
STD2_FEATURES = ["tf_exp", "tf_activity"]

FEATURES30 = [
    "maxpwm", "cons", "crup", "crup_mean", "crup_delta", "remap", "tf_exp", "tf_activity",
    "tfact_crupcor_coef", "tfact_crupcor_pval", "trap", "atac_min", "atac_max", "atac_mean", "atac_mean_mean",
    "atac_delta_min", "atac_delta_max", "atac_delta_mean", "tobias_avg", "delta_tobias_avg", "tobias_mean_mean", "tobias_count",
    "cot_hits_0", "cot_hits_1", "cot_hits_2", "cot_hits_3",
    "cot_maxpwm_0", "cot_maxpwm_1", "cot_maxpwm_2", "cot_maxpwm_3",
]

TFONLY_FEATS = FEATURES30
TFTUNED_FEATS = FEATURES30 + ["xgb_general_model"]
TFTRANS_FEATS = EMB_COLS + FEATURES30


def setup_logging(level: str = "INFO") -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%H:%M:%S",
    )


def ensure_dirs() -> None:
    for path in [MODEL_ROOT, GENERAL_DIR, TFONLY_DIR, TFTUNED_DIR, TFTRANS_DIR]:
        os.makedirs(path, exist_ok=True)


def ensure_columns(df: pd.DataFrame, cols: List[str], fill: float = 0.0) -> pd.DataFrame:
    miss = [c for c in cols if c not in df.columns]
    if not miss:
        return df
    return pd.concat([df, pd.DataFrame(fill, index=df.index, columns=miss)], axis=1)


def dmatrix_from_df(df: pd.DataFrame, label: np.ndarray | None = None) -> xgb.DMatrix:
    arr = df.to_numpy(dtype=np.float32, copy=False)
    if label is None:
        return xgb.DMatrix(arr, feature_names=list(df.columns))
    return xgb.DMatrix(arr, label=label, feature_names=list(df.columns))


def get_row_enh_ids(df: pd.DataFrame) -> pd.Series:
    if "enh_id" in df.columns:
        return df["enh_id"].astype(str)
    return pd.Series(df.index.astype(str), index=df.index)


def split_two(values: List[str]) -> Tuple[List[str], List[str]]:
    vals = sorted(list(dict.fromkeys(values)))
    if len(vals) <= 1:
        return vals, vals
    mid = len(vals) // 2
    a, b = vals[:mid], vals[mid:]
    if len(a) == 0:
        a = b
    if len(b) == 0:
        b = a
    return a, b


def split_random_two(values: List[str], rng: np.random.Generator) -> Tuple[List[str], List[str]]:
    vals = list(dict.fromkeys(values))
    if len(vals) <= 1:
        return vals, vals
    perm = list(rng.permutation(vals))
    mid = len(perm) // 2
    a, b = perm[:mid], perm[mid:]
    if len(a) == 0:
        a = b
    if len(b) == 0:
        b = a
    return a, b


def split_random_k(values: List[str], k: int, rng: np.random.Generator) -> List[List[str]]:
    vals = list(dict.fromkeys(values))
    if len(vals) == 0:
        return [[] for _ in range(k)]
    perm = list(rng.permutation(vals))
    out = [[] for _ in range(k)]
    for i, v in enumerate(perm):
        out[i % k].append(v)
    return [g for g in out if len(g) > 0]


def load_enhancer_chr_map() -> pd.Series:
    enh = pd.read_csv(
        ENHANCERS_BED,
        sep="\t",
        header=None,
        usecols=[0, 3],
        names=["chr", "enh_id"],
        dtype=str,
    )
    return enh.drop_duplicates(subset=["enh_id"], keep="first").set_index("enh_id")["chr"]


def load_training_sets() -> Dict[str, pd.DataFrame]:
    with open(TRAINING_SETS_PK, "rb") as f:
        return pickle.load(f)


def load_embeddings() -> pd.DataFrame:
    nt = pd.read_csv(EMBEDDINGS_CSV)
    nt.index = nt["enh_id"].astype(str)
    drop_cols = [c for c in ["Unnamed: 0", "chr", "start", "end", "enh_id", "coord"] if c in nt.columns]
    nt = nt.drop(columns=drop_cols)
    nt.columns = EMB_COLS
    nt = nt.apply(pd.to_numeric, errors="coerce").fillna(0.0).astype(np.float32)
    return nt


def compute_std2_stats(df: pd.DataFrame) -> Dict[str, float]:
    stats = {}
    for feat in STD2_FEATURES:
        arr = pd.to_numeric(df[feat], errors="coerce").fillna(0.0).to_numpy(dtype=np.float32)
        mean = float(np.mean(arr))
        std = float(np.std(arr))
        if std == 0.0:
            std = 1.0
        stats[f"{feat}_mean"] = mean
        stats[f"{feat}_std"] = std
    return stats


def apply_std2(df: pd.DataFrame, stats: Dict[str, float]) -> pd.DataFrame:
    out = df.copy()
    for feat in STD2_FEATURES:
        mean = float(stats.get(f"{feat}_mean", 0.0))
        std = float(stats.get(f"{feat}_std", 1.0))
        if std == 0.0:
            std = 1.0
        out[feat] = (pd.to_numeric(out[feat], errors="coerce").fillna(0.0) - mean) / std
    return out


def train_one_fold(
    df: pd.DataFrame,
    features: List[str],
    train_mask: np.ndarray,
    valid_mask: np.ndarray,
    out_json: str,
    seed: int,
    device: str,
    max_rounds: int,
    early_stopping_rounds: int,
) -> Dict[str, Any]:
    if int(train_mask.sum()) == 0 or int(valid_mask.sum()) == 0:
        raise RuntimeError("Empty fold split")

    work = ensure_columns(df.copy(), features, 0.0)
    work[features] = work[features].apply(pd.to_numeric, errors="coerce").fillna(0.0)

    X_tr = work.loc[train_mask, features].copy()
    X_va = work.loc[valid_mask, features].copy()
    y_tr = work.loc[train_mask, "label_bin"].to_numpy(dtype=np.float32)
    y_va = work.loc[valid_mask, "label_bin"].to_numpy(dtype=np.float32)

    std_stats = compute_std2_stats(X_tr)
    X_tr = apply_std2(X_tr, std_stats)
    X_va = apply_std2(X_va, std_stats)

    dtrain = dmatrix_from_df(X_tr, y_tr)
    dvalid = dmatrix_from_df(X_va, y_va)

    params = {
        "booster": "gbtree",
        "objective": "reg:logistic",
        "learning_rate": 0.05,
        "subsample": 0.8,
        "colsample_bytree": 0.8,
        "max_depth": 6,
        "min_child_weight": 1,
        "eval_metric": "logloss",
        "random_state": seed,
        "device": device,
    }

    booster = xgb.train(
        params,
        dtrain,
        num_boost_round=max_rounds,
        evals=[(dvalid, "valid")],
        early_stopping_rounds=early_stopping_rounds,
        verbose_eval=False,
    )

    os.makedirs(os.path.dirname(out_json), exist_ok=True)
    booster.save_model(out_json)

    return {
        "model_path": out_json,
        "train_rows": int(train_mask.sum()),
        "valid_rows": int(valid_mask.sum()),
        "train_pos_rate": float(np.mean(y_tr)),
        "valid_pos_rate": float(np.mean(y_va)),
        "best_iteration": int(booster.best_iteration) if booster.best_iteration is not None else int(max_rounds),
        "best_score": float(booster.best_score) if booster.best_score is not None else np.nan,
        **std_stats,
    }


def build_general_matrix(training_sets: Dict[str, pd.DataFrame], enh_chr: pd.Series) -> pd.DataFrame:
    rows = []
    for tf in sorted(training_sets.keys()):
        raw = training_sets[tf].copy()
        if "tissue" not in raw.columns or "label" not in raw.columns:
            continue
        raw = ensure_columns(raw, FEATURES30, 0.0)
        raw[FEATURES30] = raw[FEATURES30].apply(pd.to_numeric, errors="coerce").fillna(0.0).astype(np.float32)

        enh_ids = get_row_enh_ids(raw)
        row_chr = enh_ids.map(enh_chr).fillna("chrUNK").astype(str)
        label_bin = (pd.to_numeric(raw["label"], errors="coerce").fillna(0.0) > 2).astype(np.int8)

        part = pd.DataFrame(
            {
                "tf": tf,
                "tissue": raw["tissue"].astype(str).to_numpy(),
                "chr": row_chr.to_numpy(),
                "label_bin": label_bin.to_numpy(),
            }
        )
        part = pd.concat([part.reset_index(drop=True), raw[FEATURES30].reset_index(drop=True)], axis=1)
        rows.append(part)

    if len(rows) == 0:
        raise RuntimeError("No rows available for general matrix")
    out = pd.concat(rows, axis=0, ignore_index=True)
    LOGGER.info(
        "General matrix | rows=%d TFs=%d tissues=%d chrs=%d pos=%.5f",
        len(out),
        int(out["tf"].nunique()),
        int(out["tissue"].nunique()),
        int(out["chr"].nunique()),
        float(out["label_bin"].to_numpy(dtype=np.float32).mean()),
    )
    return out


def train_general_model(
    general_df: pd.DataFrame,
    seed: int,
    max_rounds: int,
    early_stopping_rounds: int,
    device: str,
) -> None:
    rng = np.random.default_rng(seed)
    tfs = sorted(general_df["tf"].astype(str).unique().tolist())
    tf_groups = split_random_k(tfs, 4, rng)
    all_tfs = set(tfs)

    chrs = sorted(general_df["chr"].astype(str).unique().tolist())
    chr_a, chr_b = split_random_two(chrs, rng)
    chr_halves = [set(chr_a), set(chr_b)]

    rows = []
    fold_i = 0
    for gi, valid_tf_group in enumerate(tf_groups):
        valid_tf = set(valid_tf_group)
        train_tf = all_tfs.difference(valid_tf)
        for ci in [0, 1]:
            valid_chr = chr_halves[ci]
            train_chr = chr_halves[1 - ci]
            tr = (general_df["tf"].isin(train_tf) & general_df["chr"].isin(train_chr)).to_numpy(dtype=bool)
            va = (general_df["tf"].isin(valid_tf) & general_df["chr"].isin(valid_chr)).to_numpy(dtype=bool)
            fold_name = f"g{gi+1}_c{ci+1}"

            info = train_one_fold(
                df=general_df,
                features=FEATURES30,
                train_mask=tr,
                valid_mask=va,
                out_json=f"{GENERAL_DIR}/{fold_name}.json",
                seed=seed + fold_i,
                device=device,
                max_rounds=max_rounds,
                early_stopping_rounds=early_stopping_rounds,
            )
            info.update({"fold": fold_name, "valid_tf_group": gi + 1, "valid_chr_half": ci + 1})
            rows.append(info)
            LOGGER.info("General fold=%s train=%d valid=%d", fold_name, info["train_rows"], info["valid_rows"])
            fold_i += 1

    summary_df = pd.DataFrame(rows)
    summary_df.to_csv(f"{GENERAL_DIR}/fold_summary.tsv", sep="\t", index=False)
    summary_df[["fold", "model_path", "tf_exp_mean", "tf_exp_std", "tf_activity_mean", "tf_activity_std"]].to_csv(
        GENERAL_STD2_TSV,
        sep="\t",
        index=False,
    )
    with open(f"{GENERAL_DIR}/metadata.json", "w") as f:
        json.dump({"tf_groups": tf_groups, "chr_split": {"A": chr_a, "B": chr_b}}, f, indent=2)


def build_tf_folds(row_tissue: pd.Series, row_chr: pd.Series) -> List[Tuple[str, List[str], List[str], List[str], List[str]]]:
    tissues = sorted(row_tissue.dropna().astype(str).unique().tolist())
    chrs = sorted(row_chr.dropna().astype(str).unique().tolist())
    ta, tb = split_two(tissues)
    c1, c2 = split_two(chrs)

    if len(tissues) <= 1:
        return [
            ("c1", tissues, c1, tissues, c2),
            ("c2", tissues, c2, tissues, c1),
        ]

    return [
        ("a1", ta, c1, tb, c2),
        ("a2", ta, c2, tb, c1),
        ("b1", tb, c1, ta, c2),
        ("b2", tb, c2, ta, c1),
    ]


def get_feature_names_for_booster(booster: xgb.Booster, fallback: List[str]) -> List[str]:
    try:
        names = booster.feature_names
        if names is not None:
            return list(names)
    except Exception:
        pass
    return list(fallback)


def predict_booster(booster: xgb.Booster, X: pd.DataFrame, fallback_features: List[str]) -> np.ndarray:
    feats = get_feature_names_for_booster(booster, fallback_features)
    X2 = ensure_columns(X, feats, 0.0)[feats]
    return booster.predict(dmatrix_from_df(X2))


def load_general_models_with_stats() -> List[Dict[str, Any]]:
    summary = pd.read_csv(f"{GENERAL_DIR}/fold_summary.tsv", sep="\t")
    out = []
    for _, row in summary.iterrows():
        p = str(row["model_path"])
        if not os.path.exists(p):
            continue
        m = xgb.Booster()
        m.load_model(p)
        out.append(
            {
                "model": m,
                "tf_exp_mean": float(row["tf_exp_mean"]),
                "tf_exp_std": float(row["tf_exp_std"]),
                "tf_activity_mean": float(row["tf_activity_mean"]),
                "tf_activity_std": float(row["tf_activity_std"]),
            }
        )
    return out


def predict_with_std2_models(models: List[Dict[str, Any]], X: pd.DataFrame, fallback_features: List[str]) -> np.ndarray:
    if len(models) == 0:
        return np.full(len(X), np.nan, dtype=float)
    preds = []
    for entry in models:
        stats = {
            "tf_exp_mean": entry["tf_exp_mean"],
            "tf_exp_std": entry["tf_exp_std"],
            "tf_activity_mean": entry["tf_activity_mean"],
            "tf_activity_std": entry["tf_activity_std"],
        }
        X2 = apply_std2(ensure_columns(X.copy(), fallback_features, 0.0), stats)
        preds.append(predict_booster(entry["model"], X2, fallback_features))
    return np.mean(np.vstack(preds), axis=0)


def train_tf_family(
    tf: str,
    tf_df: pd.DataFrame,
    folds: List[Tuple[str, List[str], List[str], List[str], List[str]]],
    feature_map: Dict[str, Tuple[List[str], str]],
    seed_base: int,
    max_rounds: int,
    early_stopping_rounds: int,
    device: str,
) -> List[Dict[str, Any]]:
    rows = []
    row_tissue = tf_df["tissue"].astype(str).reset_index(drop=True)
    row_chr = tf_df["chr"].astype(str).reset_index(drop=True)

    for fam_i, (family_name, (features, out_dir)) in enumerate(feature_map.items()):
        for fold_i, (fold_name, tr_tis, tr_chr, va_tis, va_chr) in enumerate(folds):
            tr = (row_tissue.isin(tr_tis) & row_chr.isin(tr_chr)).to_numpy(dtype=bool)
            va = (row_tissue.isin(va_tis) & row_chr.isin(va_chr)).to_numpy(dtype=bool)

            info = train_one_fold(
                df=tf_df,
                features=features,
                train_mask=tr,
                valid_mask=va,
                out_json=f"{out_dir}/{tf}_{fold_name}.json",
                seed=seed_base + 1000 * fam_i + fold_i,
                device=device,
                max_rounds=max_rounds,
                early_stopping_rounds=early_stopping_rounds,
            )
            info.update({"tf": tf, "family": family_name, "fold": fold_name})
            rows.append(info)
            LOGGER.info("TF=%s family=%s fold=%s train=%d valid=%d", tf, family_name, fold_name, info["train_rows"], info["valid_rows"])

    return rows


def train_tf_models(
    training_sets: Dict[str, pd.DataFrame],
    nt_embeddings: pd.DataFrame,
    enh_chr: pd.Series,
    general_models: List[Dict[str, Any]],
    seed: int,
    max_rounds: int,
    early_stopping_rounds: int,
    device: str,
) -> None:
    all_rows = []

    feature_map = {
        "tf_only": (TFONLY_FEATS, TFONLY_DIR),
        "tf_tuned": (TFTUNED_FEATS, TFTUNED_DIR),
        "tf_transformer": (TFTRANS_FEATS, TFTRANS_DIR),
    }

    for tf_i, tf in enumerate(sorted(training_sets.keys())):
        raw = training_sets[tf].copy()
        if "tissue" not in raw.columns or "label" not in raw.columns:
            continue

        enh_ids = get_row_enh_ids(raw)
        row_chr = enh_ids.map(enh_chr).fillna("chrUNK").astype(str)

        base = ensure_columns(raw, FEATURES30, 0.0).copy()
        base[FEATURES30] = base[FEATURES30].apply(pd.to_numeric, errors="coerce").fillna(0.0).astype(np.float32)
        base["tissue"] = raw["tissue"].astype(str)
        base["label_bin"] = (pd.to_numeric(raw["label"], errors="coerce").fillna(0.0) > 2).astype(np.int8)
        base["chr"] = row_chr.values

        emb_aligned = nt_embeddings.reindex(enh_ids.values)
        emb_aligned.index = base.index
        emb_aligned = emb_aligned.astype(np.float32)
        base = pd.concat([base, emb_aligned], axis=1)
        base = ensure_columns(base, EMB_COLS, 0.0)
        base[EMB_COLS] = base[EMB_COLS].apply(pd.to_numeric, errors="coerce").fillna(0.0).astype(np.float32)

        base["xgb_general_model"] = predict_with_std2_models(general_models, base, FEATURES30).astype(np.float32)

        folds = build_tf_folds(base["tissue"], base["chr"])
        if len(folds) < 4:
            LOGGER.warning("TF=%s has %d folds (single-tissue fallback)", tf, len(folds))

        tf_rows = train_tf_family(
            tf=tf,
            tf_df=base,
            folds=folds,
            feature_map=feature_map,
            seed_base=seed + 10000 + tf_i,
            max_rounds=max_rounds,
            early_stopping_rounds=early_stopping_rounds,
            device=device,
        )
        all_rows.extend(tf_rows)

        del raw, base, emb_aligned
        gc.collect()

    tf_summary_df = pd.DataFrame(all_rows)
    tf_summary_df.to_csv(f"{MODEL_ROOT}/tf_families_fold_summary.tsv", sep="\t", index=False)
    tf_summary_df[
        [
            "tf",
            "family",
            "fold",
            "model_path",
            "tf_exp_mean",
            "tf_exp_std",
            "tf_activity_mean",
            "tf_activity_std",
        ]
    ].to_csv(
        TF_STD2_TSV,
        sep="\t",
        index=False,
    )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Train ensemble models (general, tf_only, tf_tuned, tf_transformer)")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--max-rounds", type=int, default=1000)
    p.add_argument("--early-stopping-rounds", type=int, default=50)
    p.add_argument("--device", type=str, default="cuda", choices=["cuda", "cpu"])
    p.add_argument("--log-level", type=str, default="INFO")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    setup_logging(args.log_level)
    ensure_dirs()

    LOGGER.info("Loading training sets")
    training_sets = load_training_sets()
    enh_chr = load_enhancer_chr_map()
    nt_embeddings = load_embeddings()

    LOGGER.info("Step 1: train general ensemble")
    general_df = build_general_matrix(training_sets, enh_chr)
    train_general_model(
        general_df=general_df,
        seed=args.seed,
        max_rounds=args.max_rounds,
        early_stopping_rounds=args.early_stopping_rounds,
        device=args.device,
    )

    LOGGER.info("Loading trained general ensemble for downstream TF families")
    general_models = load_general_models_with_stats()

    LOGGER.info("Step 2: train TF family ensembles")
    train_tf_models(
        training_sets=training_sets,
        nt_embeddings=nt_embeddings,
        enh_chr=enh_chr,
        general_models=general_models,
        seed=args.seed,
        max_rounds=args.max_rounds,
        early_stopping_rounds=args.early_stopping_rounds,
        device=args.device,
    )


if __name__ == "__main__":
    main()
