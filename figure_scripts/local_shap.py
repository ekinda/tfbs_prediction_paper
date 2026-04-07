import os
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import shap
import xgboost as xgb

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
import config as cfg

TR_DIR = cfg.TF_TRANSFORMER_MODEL_DIR
OUTDIR = cfg.LOCAL_SHAP_OUTPUT_DIR

selected_tf = "HNF4A"
selected_tissue = "liver"
enhancer_ids = ["EH38E1355849", "EH38E2817083", "EH38E2817084", "EH38E1355852"]
top_k = 8

emb_cols = [f"emb_{i}" for i in range(1024)]
gen_features = [
    "maxpwm", "cons", "crup", "crup_mean", "crup_delta", "remap", "tf_exp", "tf_activity",
    "tfact_crupcor_coef", "tfact_crupcor_pval", "trap", "atac_min", "atac_max", "atac_mean", "atac_mean_mean",
    "atac_delta_min", "atac_delta_max", "atac_delta_mean", "tobias_avg", "delta_tobias_avg", "tobias_mean_mean", "tobias_count",
    "cot_hits_0", "cot_hits_1", "cot_hits_2", "cot_hits_3",
    "cot_maxpwm_0", "cot_maxpwm_1", "cot_maxpwm_2", "cot_maxpwm_3",
]
nt_features = emb_cols + gen_features


def ensure_columns(df: pd.DataFrame, cols, fill=0.0) -> pd.DataFrame:
    miss = [c for c in cols if c not in df.columns]
    if miss:
        df = df.copy()
        for c in miss:
            df[c] = fill
    return df


def dmat_from_df(df: pd.DataFrame):
    arr = df.to_numpy(dtype=np.float32, copy=False)
    try:
        return xgb.DMatrix(arr, feature_names=list(df.columns))
    except Exception:
        return xgb.DMatrix(arr)


def load_tr_models(tf: str):
    first = [f"{TR_DIR}/{tf}_{s}.json" for s in ["a1", "a2", "b1", "b2"]]
    fallback = [f"{TR_DIR}/{tf}_{s}.json" for s in ["a", "b"]]
    paths = first if all(os.path.exists(p) for p in first) else fallback
    if not all(os.path.exists(p) for p in paths):
        raise FileNotFoundError(f"No complete transformer ensemble for {tf}")
    out = []
    for p in paths:
        m = xgb.Booster()
        m.load_model(p)
        out.append(m)
    return out


def ensemble_local_contribs(models, X: pd.DataFrame):
    d = dmat_from_df(X)
    stack = []
    for m in models:
        c = m.predict(d, pred_contribs=True)
        if c.ndim == 3:
            c = c[:, 0, :] if c.shape[1] == 1 else c[:, 1, :]
        stack.append(c)
    cmean = np.mean(np.stack(stack, axis=0), axis=0)
    shap_vals = cmean[:, :-1]
    base_vals = cmean[:, -1]
    scores = base_vals + shap_vals.sum(axis=1)
    return shap_vals, base_vals, scores


def save_shap_waterfall(shap_s: pd.Series, feat_values: pd.Series, base: float, top_k: int, title: str, out_png: str):
    explanation = shap.Explanation(
        values=shap_s.to_numpy(dtype=float),
        base_values=float(base),
        data=feat_values.loc[shap_s.index].to_numpy(dtype=float),
        feature_names=list(shap_s.index),
    )

    plt.figure(figsize=(9, 6))
    shap.plots.waterfall(explanation, max_display=top_k, show=False)
    ax = plt.gca()
    ax.grid(False)
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    print("Saved:", out_png)
    plt.close()


def run() -> pd.DataFrame:
    plt.style.use("default")
    plt.rcParams["axes.grid"] = False
    os.makedirs(OUTDIR, exist_ok=True)

    with open(cfg.TEST_SETS_PK, "rb") as f:
        test_sets = pickle.load(f)
    for key in test_sets:
        test_sets[key]["label"] = test_sets[key]["label"] > 2

    nt = pd.read_csv(cfg.EMBEDDINGS_CSV)
    nt.index = nt["enh_id"].astype(str)
    drop_cols = [c for c in ["Unnamed: 0", "chr", "start", "end", "enh_id", "coord"] if c in nt.columns]
    nt = nt.drop(columns=drop_cols)
    nt.columns = emb_cols

    cur = test_sets[selected_tf][test_sets[selected_tf]["tissue"].astype(str) == selected_tissue].copy()
    cur.index = cur.index.astype(str)
    cur = cur.merge(nt, left_index=True, right_index=True, how="left")
    cur = ensure_columns(cur, nt_features, 0.0)

    region_df = cur.loc[[eid for eid in enhancer_ids if eid in cur.index]].copy()
    if len(region_df) == 0:
        raise ValueError("None of the requested enhancer IDs found for selected TF/tissue")

    models = load_tr_models(selected_tf)
    S, B, _ = ensemble_local_contribs(models, region_df[nt_features])

    print("Enhancers found:", len(region_df))
    print(region_df[["chr", "start", "end", "label"]])

    rows = []

    for i, (eid, row) in enumerate(region_df.iterrows()):
        shap_s = pd.Series(S[i], index=nt_features, dtype=float)
        final = float(B[i] + shap_s.sum())

        out_png = f"{OUTDIR}/local_shap_waterfall_{eid}_tf_transformer.png"
        save_shap_waterfall(
            shap_s=shap_s,
            feat_values=row[nt_features],
            base=float(B[i]),
            top_k=top_k,
            title=f"{eid} ({selected_tf}, {selected_tissue})",
            out_png=out_png,
        )

        top = shap_s.abs().sort_values(ascending=False).head(top_k).index
        for f in top:
            rows.append(
                {
                    "enh_id": eid,
                    "feature": f,
                    "feature_value": row[f],
                    "shap_contribution": float(shap_s[f]),
                    "baseline": float(B[i]),
                    "score": float(final),
                }
            )

    top_df = pd.DataFrame(rows)
    out_tsv = f"{OUTDIR}/local_shap_top_features_4enhancers.tsv"
    top_df.to_csv(out_tsv, sep="\t", index=False)
    print("Saved:", out_tsv)
    print(top_df)
    return top_df


if __name__ == "__main__":
    run()
