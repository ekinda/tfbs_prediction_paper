import os
import glob
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import xgboost as xgb

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
import config as cfg

TR_DIR = cfg.TF_TRANSFORMER_MODEL_DIR
OUTDIR = cfg.GLOBAL_SHAP_OUTPUT_DIR

MODEL_KEYS = ["general", "tf_tuned", "tf_only", "tf+transformer"]
DEFAULT_SELECTED_TFS = ["CTCF", "HNF4A", "YY1", "GABPA", "REST"]

MAX_ROWS_PER_TF_TISSUE = 1000
PLOT_TOP_N = 20
RANDOM_STATE = 7
DEPENDENCE_MAX_POINTS = 50000
DEPENDENCE_TOP_GENERAL_TFONLY = 5
DEPENDENCE_TOP_EMBEDDINGS = 3

emb_cols = [f"emb_{i}" for i in range(1024)]
gen_features = [
    "maxpwm", "cons", "crup", "crup_mean", "crup_delta", "remap", "tf_exp", "tf_activity",
    "tfact_crupcor_coef", "tfact_crupcor_pval", "trap", "atac_min", "atac_max", "atac_mean", "atac_mean_mean",
    "atac_delta_min", "atac_delta_max", "atac_delta_mean", "tobias_avg", "delta_tobias_avg", "tobias_mean_mean", "tobias_count",
    "cot_hits_0", "cot_hits_1", "cot_hits_2", "cot_hits_3",
    "cot_maxpwm_0", "cot_maxpwm_1", "cot_maxpwm_2", "cot_maxpwm_3",
]
tfonly_features = [
    "maxpwm", "cons", "crup", "crup_mean", "crup_delta", "remap", "tf_exp", "tf_activity",
    "tfact_crupcor_coef", "tfact_crupcor_pval", "trap", "atac_min", "atac_max", "atac_mean", "atac_mean_mean",
    "atac_delta_min", "atac_delta_max", "atac_delta_mean", "tobias_avg", "delta_tobias_avg", "tobias_mean_mean", "tobias_count",
    "cot_hits_0", "cot_hits_1", "cot_hits_2", "cot_hits_3", "cot_hits_4", "cot_hits_5",
    "cot_maxpwm_0", "cot_maxpwm_1", "cot_maxpwm_2", "cot_maxpwm_3", "cot_maxpwm_4", "cot_maxpwm_5",
]
tftuned_default = gen_features + ["xgb_general"]
nt_features = emb_cols + gen_features


def ensure_columns(df: pd.DataFrame, cols, fill=0.0) -> pd.DataFrame:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        df = df.copy()
        for col in missing:
            df[col] = fill
    return df


def to_booster(model):
    if isinstance(model, xgb.Booster):
        return model
    if hasattr(model, "get_booster"):
        return model.get_booster()
    raise TypeError(type(model))


def get_feats(model, fallback):
    try:
        feature_names = to_booster(model).feature_names
        if feature_names is not None:
            return list(feature_names)
    except Exception:
        pass

    if hasattr(model, "feature_names_in_"):
        return list(model.feature_names_in_)

    return list(fallback)


def dmat_from_df(df: pd.DataFrame):
    arr = df.to_numpy(dtype=np.float32, copy=False)
    try:
        return xgb.DMatrix(arr, feature_names=list(df.columns))
    except Exception:
        return xgb.DMatrix(arr)


def pred_contribs(model, X: pd.DataFrame):
    contribs = to_booster(model).predict(dmat_from_df(X), pred_contribs=True)
    if contribs.ndim == 3:
        contribs = contribs[:, 0, :] if contribs.shape[1] == 1 else contribs[:, 1, :]
    return contribs[:, :-1], contribs[:, -1]


def predict_scores(model, X: pd.DataFrame) -> np.ndarray:
    return to_booster(model).predict(dmat_from_df(X))


def load_tr_models(tf_name: str):
    first = [f"{TR_DIR}/{tf_name}_{s}.json" for s in ["a1", "a2", "b1", "b2"]]
    fallback = [f"{TR_DIR}/{tf_name}_{s}.json" for s in ["a", "b"]]
    paths = first if all(os.path.exists(p) for p in first) else fallback
    if not all(os.path.exists(p) for p in paths):
        return []

    models = []
    for path in paths:
        model = xgb.Booster()
        model.load_model(path)
        models.append(model)
    return models


def ensemble_contribs(models, X: pd.DataFrame):
    shap_stack = []
    for model in models:
        shap_values, _ = pred_contribs(model, X)
        shap_stack.append(shap_values)
    return np.mean(np.stack(shap_stack, axis=0), axis=0)


def empty_agg_dict():
    return {k: {"sum_abs": pd.Series(dtype=float), "sum_signed": pd.Series(dtype=float), "n": 0} for k in MODEL_KEYS}


def update_agg(agg, key: str, feats, shap: np.ndarray):
    sum_abs = pd.Series(np.abs(shap).sum(axis=0), index=feats, dtype=float)
    sum_signed = pd.Series(shap.sum(axis=0), index=feats, dtype=float)
    agg[key]["sum_abs"] = agg[key]["sum_abs"].add(sum_abs, fill_value=0.0)
    agg[key]["sum_signed"] = agg[key]["sum_signed"].add(sum_signed, fill_value=0.0)
    agg[key]["n"] += shap.shape[0]


def agg_to_df(agg) -> pd.DataFrame:
    rows = []
    for model_key, stats in agg.items():
        n = max(int(stats["n"]), 1)
        feats = sorted(set(stats["sum_abs"].index) | set(stats["sum_signed"].index))
        for feat in feats:
            rows.append(
                {
                    "model": model_key,
                    "feature": feat,
                    "mean_abs_shap": float(stats["sum_abs"].get(feat, 0.0) / n),
                    "mean_shap": float(stats["sum_signed"].get(feat, 0.0) / n),
                    "n_rows": n,
                }
            )

    out = pd.DataFrame(rows)
    if len(out) == 0:
        return out
    out["rank"] = out.groupby("model")["mean_abs_shap"].rank(ascending=False, method="first")
    return out.sort_values(["model", "mean_abs_shap"], ascending=[True, False]).reset_index(drop=True)


def top_n_by_group(df: pd.DataFrame, group_cols, n: int) -> pd.DataFrame:
    if len(df) == 0:
        return df.copy()
    asc = [True] * len(group_cols) + [False]
    return (
        df.sort_values(group_cols + ["mean_abs_shap"], ascending=asc)
        .groupby(group_cols, group_keys=False)
        .head(n)
        .reset_index(drop=True)
    )


def top_features_for_model(global_df: pd.DataFrame, model_key: str, k: int):
    sub = global_df[global_df["model"] == model_key].nlargest(k, "mean_abs_shap")
    return sub["feature"].tolist()


def is_embedding_feature(feature_name: str) -> bool:
    return str(feature_name).startswith("emb_")


def non_embedding_features_for_model(global_df: pd.DataFrame, model_key: str):
    sub = global_df[global_df["model"] == model_key].sort_values("mean_abs_shap", ascending=False)
    return [f for f in sub["feature"].tolist() if not is_embedding_feature(f)]


def top_embedding_features_for_model(global_df: pd.DataFrame, model_key: str, k: int):
    sub = global_df[global_df["model"] == model_key]
    sub = sub[sub["feature"].apply(is_embedding_feature)].nlargest(k, "mean_abs_shap")
    return sub["feature"].tolist()


def append_dependence_parts(parts, model_key: str, feature_list, feature_names, X_df: pd.DataFrame, shap_values: np.ndarray, split_name: str, tf_name: str, tissue: str):
    if not feature_list:
        return

    feat_to_idx = {f: i for i, f in enumerate(feature_names)}
    for feat in feature_list:
        if feat not in feat_to_idx or feat not in X_df.columns:
            continue
        idx = feat_to_idx[feat]
        part = pd.DataFrame(
            {
                "model": model_key,
                "feature": feat,
                "feature_value": X_df[feat].to_numpy(dtype=float, copy=False),
                "shap_value": shap_values[:, idx].astype(float),
                "split": split_name,
                "tf": tf_name,
                "tissue": tissue,
            }
        )
        parts.append(part)


def collect_dependence_points(split_map, nt_df, general_model, tf_tuned_models, tf_only_models, tr_cache, dependence_features):
    parts = []
    pair_i = 0

    for split_name, split_sets in split_map.items():
        for tf_name in sorted(split_sets):
            if tf_name not in tf_tuned_models:
                continue

            cur_all = split_sets[tf_name].copy()
            tissues = sorted(cur_all["tissue"].astype(str).unique())

            if tf_name not in tr_cache:
                tr_cache[tf_name] = load_tr_models(tf_name)

            for tissue in tissues:
                pair_i += 1
                cur = cur_all[cur_all["tissue"].astype(str) == tissue].copy()
                cur.index = cur.index.astype(str)
                cur = cur.merge(nt_df, left_index=True, right_index=True, how="left")
                cur = cur.dropna(subset=gen_features)
                if MAX_ROWS_PER_TF_TISSUE is not None and len(cur) > MAX_ROWS_PER_TF_TISSUE:
                    cur = cur.sample(MAX_ROWS_PER_TF_TISSUE, random_state=RANDOM_STATE)
                if len(cur) == 0:
                    continue

                try:
                    gf = get_feats(general_model, gen_features)
                    cur = ensure_columns(cur, gf, 0.0)
                    Xg = cur[gf]
                    shap_general, _ = pred_contribs(general_model, Xg)
                    append_dependence_parts(
                        parts,
                        "general",
                        dependence_features.get("general", []),
                        gf,
                        Xg,
                        shap_general,
                        split_name,
                        tf_name,
                        tissue,
                    )

                    cur["xgb_general"] = predict_scores(general_model, Xg)

                    tf_tuned_model = tf_tuned_models[tf_name]
                    tf_tuned_feats = get_feats(tf_tuned_model, tftuned_default)
                    cur = ensure_columns(cur, tf_tuned_feats, 0.0)
                    Xt = cur[tf_tuned_feats]
                    shap_tftuned, _ = pred_contribs(tf_tuned_model, Xt)
                    append_dependence_parts(
                        parts,
                        "tf_tuned",
                        dependence_features.get("tf_tuned", []),
                        tf_tuned_feats,
                        Xt,
                        shap_tftuned,
                        split_name,
                        tf_name,
                        tissue,
                    )

                    if tf_name in tf_only_models:
                        try:
                            tf_only_model = tf_only_models[tf_name]
                            tf_only_feats = get_feats(tf_only_model, tfonly_features)
                            cur = ensure_columns(cur, tf_only_feats, 0.0)
                            Xo = cur[tf_only_feats]
                            shap_tfonly, _ = pred_contribs(tf_only_model, Xo)
                            append_dependence_parts(
                                parts,
                                "tf_only",
                                dependence_features.get("tf_only", []),
                                tf_only_feats,
                                Xo,
                                shap_tfonly,
                                split_name,
                                tf_name,
                                tissue,
                            )
                        except Exception as e:
                            print(f"[WARN] dependence tf_only skipped {split_name}:{tf_name}-{tissue}: {e}")

                    tr_models = tr_cache.get(tf_name, [])
                    if len(tr_models) > 0:
                        cur = ensure_columns(cur, nt_features, 0.0)
                        Xtr = cur[nt_features]
                        shap_transformer = ensemble_contribs(tr_models, Xtr)
                        append_dependence_parts(
                            parts,
                            "tf+transformer",
                            dependence_features.get("tf+transformer", []),
                            nt_features,
                            Xtr,
                            shap_transformer,
                            split_name,
                            tf_name,
                            tissue,
                        )

                except Exception as e:
                    print(f"[WARN] dependence pair failed {split_name}:{tf_name}-{tissue}: {e}")

                if pair_i % 20 == 0:
                    print("Dependence pass processed TF-tissue pairs:", pair_i)

    if parts:
        return pd.concat(parts, axis=0, ignore_index=True)

    return pd.DataFrame(columns=["model", "feature", "feature_value", "shap_value", "split", "tf", "tissue"])


def plot_topn_abs_grid(shap_df: pd.DataFrame, out_png: str, title_prefix: str, top_n: int):
    fig, axes = plt.subplots(2, 2, figsize=(18, 16))

    for ax, model_key in zip(axes.ravel(), MODEL_KEYS):
        sub = (
            shap_df[shap_df["model"] == model_key]
            .nlargest(top_n, "mean_abs_shap")
            .sort_values("mean_abs_shap", ascending=True)
            .reset_index(drop=True)
        )

        if len(sub) == 0:
            ax.text(0.5, 0.5, "No data", ha="center", va="center")
            ax.set_axis_off()
            continue

        y = np.arange(len(sub))
        ax.barh(y, sub["mean_abs_shap"], color="#2563eb", alpha=0.85, label="mean(|SHAP|)")
        ax.set_yticks(y)
        ax.set_yticklabels(sub["feature"])
        ax.set_title(f"{title_prefix} | {model_key} | top {top_n}")
        ax.set_xlabel("mean(|SHAP contribution|)")
        ax.set_ylabel("")
        ax.legend(loc="lower right", fontsize=9)

    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    print("Saved:", out_png)
    plt.show()


def plot_dependence_panels(dep_df: pd.DataFrame, model_key: str, features, out_png: str, title_prefix: str):
    features = [f for f in features if f in set(dep_df[dep_df["model"] == model_key]["feature"])]
    if not features:
        print(f"No dependence features to plot for {model_key}")
        return

    n = len(features)
    ncols = min(5, n)
    nrows = int(np.ceil(n / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 3.8 * nrows), squeeze=False)

    for i, feat in enumerate(features):
        ax = axes[i // ncols][i % ncols]
        sub = dep_df[(dep_df["model"] == model_key) & (dep_df["feature"] == feat)]

        if len(sub) == 0:
            ax.text(0.5, 0.5, "No data", ha="center", va="center")
            ax.set_axis_off()
            continue

        if len(sub) > DEPENDENCE_MAX_POINTS:
            sub = sub.sample(DEPENDENCE_MAX_POINTS, random_state=RANDOM_STATE)

        ax.scatter(
            sub["feature_value"],
            sub["shap_value"],
            s=7,
            alpha=0.25,
            color="#1d4ed8",
            edgecolors="none",
        )
        ax.axhline(0.0, color="black", linewidth=1.0)
        ax.set_xlabel(feat)
        ax.set_ylabel("SHAP value")
        ax.set_title(f"{model_key} | {feat}")

    for j in range(n, nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")

    fig.suptitle(title_prefix, fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    print("Saved:", out_png)
    plt.show()


def run(selected_tfs=None, top_n: int = PLOT_TOP_N) -> pd.DataFrame:
    sns.set(style="whitegrid", context="talk")
    np.random.seed(RANDOM_STATE)
    os.makedirs(OUTDIR, exist_ok=True)

    if selected_tfs is None:
        selected_tfs = DEFAULT_SELECTED_TFS
    requested_tfs = [str(x).upper() for x in selected_tfs]

    with open(cfg.TEST_SETS_PK, "rb") as f:
        test_sets = pickle.load(f)
    with open(cfg.LEADERBOARD_SETS_PK, "rb") as f:
        leaderboard_sets = pickle.load(f)

    for split_sets in [test_sets, leaderboard_sets]:
        for key in split_sets:
            split_sets[key]["label"] = split_sets[key]["label"] > 2

    if "NANOG" in test_sets:
        del test_sets["NANOG"]

    split_map = {
        "test": test_sets,
        "leaderboard": leaderboard_sets,
    }

    nt = pd.read_csv(cfg.EMBEDDINGS_CSV)
    nt.index = nt["enh_id"].astype(str)
    drop_cols = [c for c in ["Unnamed: 0", "chr", "start", "end", "enh_id", "coord"] if c in nt.columns]
    nt = nt.drop(columns=drop_cols)
    nt.columns = emb_cols

    general = xgb.Booster()
    general.load_model(cfg.GENERAL_BASE_MODEL_JSON)

    tf_tuned = {}
    for path in glob.glob(cfg.TF_TUNED_MODEL_GLOB):
        tf_name = os.path.basename(path).split("_")[0]
        model = xgb.Booster()
        model.load_model(path)
        tf_tuned[tf_name] = model
    tf_tuned.pop("general", None)

    with open(cfg.TF_ONLY_MODELS_PK, "rb") as f:
        tf_only = pickle.load(f)

    available_tfs = set()
    for split_sets in split_map.values():
        available_tfs.update(split_sets.keys())

    selected_present = [tf for tf in requested_tfs if tf in available_tfs]
    selected_tfs = [tf for tf in selected_present if tf in tf_tuned]
    missing_tfs = sorted(set(requested_tfs) - set(selected_present))
    missing_model_tfs = sorted(set(selected_present) - set(selected_tfs))

    print("Loaded test TFs:", len(test_sets), "| leaderboard TFs:", len(leaderboard_sets))
    print("Loaded tf-tuned models:", len(tf_tuned), "| tf-only models:", len(tf_only))
    print("Selected TFs for breakdown:", selected_tfs)
    if missing_tfs:
        print("Selected TFs absent from both splits:", missing_tfs)
    if missing_model_tfs:
        print("Selected TFs without tf-tuned models:", missing_model_tfs)

    global_agg = empty_agg_dict()
    per_tf_agg = {tf: empty_agg_dict() for tf in selected_tfs}

    tr_cache = {}
    pair_i = 0

    for split_name, split_sets in split_map.items():
        for tf_name in sorted(split_sets):
            if tf_name not in tf_tuned:
                continue

            cur_all = split_sets[tf_name].copy()
            tissues = sorted(cur_all["tissue"].astype(str).unique())

            if tf_name not in tr_cache:
                tr_cache[tf_name] = load_tr_models(tf_name)

            for tissue in tissues:
                pair_i += 1
                cur = cur_all[cur_all["tissue"].astype(str) == tissue].copy()
                cur.index = cur.index.astype(str)
                cur = cur.merge(nt, left_index=True, right_index=True, how="left")
                cur = cur.dropna(subset=gen_features)
                if MAX_ROWS_PER_TF_TISSUE is not None and len(cur) > MAX_ROWS_PER_TF_TISSUE:
                    cur = cur.sample(MAX_ROWS_PER_TF_TISSUE, random_state=RANDOM_STATE)
                if len(cur) == 0:
                    continue

                try:
                    gf = get_feats(general, gen_features)
                    cur = ensure_columns(cur, gf, 0.0)
                    Xg = cur[gf]
                    shap_general, _ = pred_contribs(general, Xg)
                    update_agg(global_agg, "general", gf, shap_general)
                    if tf_name in per_tf_agg:
                        update_agg(per_tf_agg[tf_name], "general", gf, shap_general)

                    cur["xgb_general"] = predict_scores(general, Xg)

                    tf_tuned_model = tf_tuned[tf_name]
                    tf_tuned_feats = get_feats(tf_tuned_model, tftuned_default)
                    cur = ensure_columns(cur, tf_tuned_feats, 0.0)
                    shap_tftuned, _ = pred_contribs(tf_tuned_model, cur[tf_tuned_feats])
                    update_agg(global_agg, "tf_tuned", tf_tuned_feats, shap_tftuned)
                    if tf_name in per_tf_agg:
                        update_agg(per_tf_agg[tf_name], "tf_tuned", tf_tuned_feats, shap_tftuned)

                    if tf_name in tf_only:
                        try:
                            tf_only_model = tf_only[tf_name]
                            tf_only_feats = get_feats(tf_only_model, tfonly_features)
                            cur = ensure_columns(cur, tf_only_feats, 0.0)
                            shap_tfonly, _ = pred_contribs(tf_only_model, cur[tf_only_feats])
                            update_agg(global_agg, "tf_only", tf_only_feats, shap_tfonly)
                            if tf_name in per_tf_agg:
                                update_agg(per_tf_agg[tf_name], "tf_only", tf_only_feats, shap_tfonly)
                        except Exception as e:
                            print(f"[WARN] tf_only skipped {split_name}:{tf_name}-{tissue}: {e}")

                    tr_models = tr_cache.get(tf_name, [])
                    if len(tr_models) > 0:
                        cur = ensure_columns(cur, nt_features, 0.0)
                        shap_transformer = ensemble_contribs(tr_models, cur[nt_features])
                        update_agg(global_agg, "tf+transformer", nt_features, shap_transformer)
                        if tf_name in per_tf_agg:
                            update_agg(per_tf_agg[tf_name], "tf+transformer", nt_features, shap_transformer)

                except Exception as e:
                    print(f"[WARN] pair failed {split_name}:{tf_name}-{tissue}: {e}")

                if pair_i % 20 == 0:
                    print("Processed TF-tissue pairs:", pair_i)

    global_shap_df = agg_to_df(global_agg)
    global_shap_df["tf_scope"] = "ALL"
    global_shap_df["split_scope"] = "test+leaderboard"

    per_tf_frames = []
    for tf_name in selected_tfs:
        tf_df = agg_to_df(per_tf_agg[tf_name])
        if len(tf_df) == 0:
            continue
        tf_df["tf_scope"] = tf_name
        tf_df["split_scope"] = "test+leaderboard"
        per_tf_frames.append(tf_df)

    tf_breakdown_df = pd.concat(per_tf_frames, axis=0, ignore_index=True) if per_tf_frames else pd.DataFrame(columns=global_shap_df.columns)
    combined_df = pd.concat([global_shap_df, tf_breakdown_df], axis=0, ignore_index=True)

    all_out_tsv = f"{OUTDIR}/global_shap_test_leaderboard_all4models.tsv"
    combined_out_tsv = f"{OUTDIR}/global_shap_test_leaderboard_all4models_with_selected_tf_breakdown.tsv"
    top_global_tsv = f"{OUTDIR}/global_shap_top{top_n}_test_leaderboard_all4models.tsv"
    top_tf_tsv = f"{OUTDIR}/global_shap_top{top_n}_selected_tfs_test_leaderboard_all4models.tsv"

    global_shap_df.to_csv(all_out_tsv, sep="\t", index=False)
    combined_df.to_csv(combined_out_tsv, sep="\t", index=False)

    top_global_df = top_n_by_group(global_shap_df, ["model"], top_n)
    top_global_df.to_csv(top_global_tsv, sep="\t", index=False)

    top_tf_df = top_n_by_group(tf_breakdown_df, ["tf_scope", "model"], top_n)
    top_tf_df.to_csv(top_tf_tsv, sep="\t", index=False)

    print("Saved:", all_out_tsv)
    print("Saved:", combined_out_tsv)
    print("Saved:", top_global_tsv)
    print("Saved:", top_tf_tsv)

    overall_png = f"{OUTDIR}/global_shap_top{top_n}_abs_test_leaderboard_all4models.png"
    plot_topn_abs_grid(global_shap_df, overall_png, "ALL TFs (test+leaderboard)", top_n)

    for tf_name in selected_tfs:
        tf_subset = tf_breakdown_df[tf_breakdown_df["tf_scope"] == tf_name]
        if len(tf_subset) == 0:
            continue
        tf_png = f"{OUTDIR}/global_shap_top{top_n}_abs_{tf_name.lower()}_test_leaderboard_all4models.png"
        plot_topn_abs_grid(tf_subset, tf_png, f"{tf_name} (test+leaderboard)", top_n)

    dependence_features = {
        "general": top_features_for_model(global_shap_df, "general", DEPENDENCE_TOP_GENERAL_TFONLY),
        "tf_only": top_features_for_model(global_shap_df, "tf_only", DEPENDENCE_TOP_GENERAL_TFONLY),
        "tf_tuned": list(
            dict.fromkeys(
                non_embedding_features_for_model(global_shap_df, "tf_tuned")
                + top_embedding_features_for_model(global_shap_df, "tf_tuned", DEPENDENCE_TOP_EMBEDDINGS)
            )
        ),
        "tf+transformer": list(
            dict.fromkeys(
                non_embedding_features_for_model(global_shap_df, "tf+transformer")
                + top_embedding_features_for_model(global_shap_df, "tf+transformer", DEPENDENCE_TOP_EMBEDDINGS)
            )
        ),
    }

    for model_key in MODEL_KEYS:
        present = set(global_shap_df[global_shap_df["model"] == model_key]["feature"])
        dependence_features[model_key] = [f for f in dependence_features.get(model_key, []) if f in present]

    print("Dependence feature sets:", dependence_features)

    dependence_df = collect_dependence_points(
        split_map=split_map,
        nt_df=nt,
        general_model=general,
        tf_tuned_models=tf_tuned,
        tf_only_models=tf_only,
        tr_cache=tr_cache,
        dependence_features=dependence_features,
    )

    dependence_out_tsv = f"{OUTDIR}/global_shap_dependence_points_test_leaderboard.tsv"
    dependence_df.to_csv(dependence_out_tsv, sep="\t", index=False)
    print("Saved:", dependence_out_tsv)

    dep_plot_specs = [
        (
            "general",
            dependence_features["general"],
            f"{OUTDIR}/dependence_general_top5_test_leaderboard.png",
            "SHAP dependence | general | top 5",
        ),
        (
            "tf_only",
            dependence_features["tf_only"],
            f"{OUTDIR}/dependence_tf_only_top5_test_leaderboard.png",
            "SHAP dependence | tf_only | top 5",
        ),
        (
            "tf_tuned",
            dependence_features["tf_tuned"],
            f"{OUTDIR}/dependence_tf_tuned_nonnt_plus_top3nt_test_leaderboard.png",
            "SHAP dependence | tf_tuned | non-NT + top 3 NT",
        ),
        (
            "tf+transformer",
            dependence_features["tf+transformer"],
            f"{OUTDIR}/dependence_tf_transformer_nonnt_plus_top3nt_test_leaderboard.png",
            "SHAP dependence | tf+transformer | non-NT + top 3 NT",
        ),
    ]

    for model_key, feat_list, out_png, title in dep_plot_specs:
        plot_dependence_panels(dependence_df, model_key, feat_list, out_png, title)

    return combined_df


if __name__ == "__main__":
    run()
