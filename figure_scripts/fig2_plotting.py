#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import math
import os
import pickle
import sys
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from scipy.stats import ttest_rel
from sklearn.metrics import average_precision_score, precision_recall_curve

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
import config as cfg

DEFAULT_INPUT = cfg.FIG2_INFERENCE_EXPORT_PK
DEFAULT_PNG = str(Path(cfg.FIGURES_DIR) / "fig2.png")
DEFAULT_PDF = str(Path(cfg.FIGURES_DIR) / "fig2.pdf")

MODEL_ORDER = [
    "general",
    "tf_only",
    "tf_tuned",
    "tf_only_cc",
    "tf_tuned_cc",
    "tf_transformer",
]

DISPLAY_NAMES = {
    "general": "General",
    "tf_only": "TF-only",
    "tf_tuned": "TF-tuned",
    "tf_only_cc": "TF-only\n(ensemble)",
    "tf_tuned_cc": "TF-tuned\n(ensemble)",
    "tf_transformer": "TF+transformer",
}

PRED_COLS = {
    "general": "pred_general",
    "tf_only": "pred_tf_only",
    "tf_tuned": "pred_tf_tuned",
    "tf_only_cc": "pred_tf_only_cc",
    "tf_tuned_cc": "pred_tf_tuned_cc",
    "tf_transformer": "pred_tf_transformer",
}

PRED_COL_ALIASES = {
    "pred_tf_only_cc": ["pred_tf_only_cc", "pred_tf_only_ensemble"],
    "pred_tf_tuned_cc": ["pred_tf_tuned_cc", "pred_tf_tuned_ensemble"],
}

COLORS = {
    "general": "#f4a259",
    "tf_only": "#b2182b",
    "tf_tuned": "#1b9e77",
    "tf_only_cc": "#d66c7a",
    "tf_tuned_cc": "#63c29c",
    "tf_transformer": "#9467bd",
}

BASELINE_COLORS = {
    "no_skill": "#9e9e9e",
    "motif_scan": "#7f7f7f",
    "atac_seq": "#5f5f5f",
}

TTEST_COMPARISONS = [
    ("general", "tf_only"),
    ("tf_only", "tf_only_cc"),
    ("tf_tuned", "tf_tuned_cc"),
    ("tf_only_cc", "tf_transformer"),
    ("tf_transformer", "tf_tuned_cc"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Figure 2 from paper revision evaluation output")
    parser.add_argument("--input", default=DEFAULT_INPUT, help="Input pickle from paper_revision_evaluation.py")
    parser.add_argument("--out-png", default=DEFAULT_PNG)
    parser.add_argument("--out-pdf", default=DEFAULT_PDF)
    return parser.parse_args()


def resolve_prediction_columns(df: pd.DataFrame) -> Dict[str, str]:
    resolved = dict(PRED_COLS)
    for canonical, options in PRED_COL_ALIASES.items():
        for opt in options:
            if opt in df.columns:
                for model_name, col_name in list(resolved.items()):
                    if col_name == canonical:
                        resolved[model_name] = opt
                break
    return resolved


def aupr_safe(y_true: np.ndarray, y_score: np.ndarray) -> float:
    mask = np.isfinite(y_score)
    if int(mask.sum()) == 0:
        return np.nan
    y_sub = y_true[mask]
    s_sub = y_score[mask]
    if np.unique(y_sub).size < 2:
        return np.nan
    return float(average_precision_score(y_sub, s_sub))


def compute_pair_metrics(combined_df: pd.DataFrame, pred_cols: Dict[str, str]) -> pd.DataFrame:
    rows = []
    for pair_id, g in combined_df.groupby("pair_id", sort=False):
        y = g["label_bin"].to_numpy(dtype=np.int8)
        row = {
            "pair_id": pair_id,
            "tf": str(g["tf"].iloc[0]),
            "tissue": str(g["tissue"].iloc[0]),
            "split": str(g["split"].iloc[0]),
            "no_skill": float(np.mean(y)),
            "motif_scan": aupr_safe(y, g["motif_scan"].to_numpy(dtype=np.float32)),
            "atac_seq": aupr_safe(y, g["atac_seq"].to_numpy(dtype=np.float32)),
            "n_rows": int(len(y)),
        }
        for m in MODEL_ORDER:
            row[m] = aupr_safe(y, g[pred_cols[m]].to_numpy(dtype=np.float32))
        rows.append(row)
    return pd.DataFrame(rows)


def compute_global_pr(combined_df: pd.DataFrame, pred_cols: Dict[str, str]) -> Dict[str, Dict[str, np.ndarray | float]]:
    out: Dict[str, Dict[str, np.ndarray | float]] = {}
    y = combined_df["label_bin"].to_numpy(dtype=np.int8)
    for model in MODEL_ORDER:
        s = combined_df[pred_cols[model]].to_numpy(dtype=np.float32)
        mask = np.isfinite(s)
        if int(mask.sum()) == 0:
            continue
        y_sub = y[mask]
        s_sub = s[mask]
        if np.unique(y_sub).size < 2:
            continue
        p, r, _ = precision_recall_curve(y_sub, s_sub)
        out[model] = {
            "precision": p,
            "recall": r,
            "aupr": float(average_precision_score(y_sub, s_sub)),
        }
    return out


def add_baseline_lines(
    ax,
    no_skill: float,
    motif_mean: float,
    atac_mean: float,
    x_min: float,
    x_max: float,
    no_skill_draw_y: float | None = None,
) -> None:
    lines = [
        ("No-skill", no_skill_draw_y if no_skill_draw_y is not None else no_skill, BASELINE_COLORS["no_skill"], no_skill, 1.0, 0.85),
        ("Motif", motif_mean, BASELINE_COLORS["motif_scan"], motif_mean, 0.8, 0.45),
        ("ATAC", atac_mean, BASELINE_COLORS["atac_seq"], atac_mean, 0.8, 0.45),
    ]
    for i, (name, y_draw, color, y_label, lw, alpha) in enumerate(lines):
        if pd.isna(y_draw):
            continue
        ax.axhline(y_draw, linestyle="--", linewidth=lw, color=color, alpha=alpha, zorder=1)
        ax.text(x_min + (x_max - x_min) * 0.01, y_draw + 0.006 + i * 0.0015, f"{name} ({y_label:.3f})", color=color, fontsize=9)


def format_p_number(p: float) -> str:
    if pd.isna(p):
        return ""
    if p < 0.05:
        return "*"
    return ""


def add_sig_bracket(ax, x1: float, x2: float, y: float, text: str, h: float = 0.010) -> None:
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color="black", lw=1.0)
    if text != "":
        ax.text((x1 + x2) / 2, y + h + 0.003, text, ha="center", va="bottom", fontsize=10, fontweight="bold")


def run_selected_paired_ttests(pair_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for a, b in TTEST_COMPARISONS:
        if a not in pair_df.columns or b not in pair_df.columns:
            rows.append({"model_a": a, "model_b": b, "n_pairs": 0, "p_value": np.nan})
            continue
        sub = pair_df[[a, b]].dropna()
        n = len(sub)
        if n < 2:
            rows.append({"model_a": a, "model_b": b, "n_pairs": n, "p_value": np.nan})
            continue
        _, p_val = ttest_rel(sub[a], sub[b], alternative="two-sided", nan_policy="omit")
        rows.append({"model_a": a, "model_b": b, "n_pairs": n, "p_value": float(p_val)})
    return pd.DataFrame(rows)


def plot_tf_panel(ax, pair_df: pd.DataFrame, tf_subset: List[str]) -> None:
    if len(tf_subset) == 0:
        ax.axis("off")
        return

    sub = pair_df[pair_df["tf"].isin(tf_subset)].copy()
    if sub.empty:
        ax.axis("off")
        return

    tf_to_x = {tf: i for i, tf in enumerate(tf_subset)}
    offsets = np.linspace(-0.30, 0.30, len(MODEL_ORDER))

    y_values = []

    for tf in tf_subset:
        cur = sub[sub["tf"] == tf].copy()
        for _, row in cur.iterrows():
            xs = []
            ys = []
            bx = tf_to_x[tf]
            tissue = str(row["tissue"])
            for i, m in enumerate(MODEL_ORDER):
                y = row[m]
                if pd.isna(y):
                    continue
                x = bx + offsets[i]
                xs.append(x)
                ys.append(float(y))
                y_values.append(float(y))
                ax.scatter(x, y, s=20, color=COLORS[m], edgecolor="white", linewidth=0.4, zorder=3)
            if len(xs) >= 2:
                ax.plot(xs, ys, color="#8a8a8a", alpha=0.65, linewidth=0.9, zorder=2)
                ax.text(xs[-1] + 0.03, ys[-1], tissue, fontsize=6, color="#666666", va="center")

        bx = tf_to_x[tf]
        ns = float(cur["no_skill"].mean()) if len(cur) else np.nan
        motif = float(cur["motif_scan"].mean()) if len(cur) else np.nan
        atac = float(cur["atac_seq"].mean()) if len(cur) else np.nan
        for yy in [ns, motif, atac]:
            if np.isfinite(yy):
                y_values.append(float(yy))

        x1 = bx + offsets[0]
        x2 = bx + offsets[-1]
        label_x = x2 + 0.015
        baselines = [
            ("No-skill", ns, BASELINE_COLORS["no_skill"]),
            ("Motif", motif, BASELINE_COLORS["motif_scan"]),
            ("ATAC", atac, BASELINE_COLORS["atac_seq"]),
        ]
        for name, yy, color in baselines:
            if not np.isfinite(yy):
                continue
            ax.plot([x1, x2], [yy, yy], linestyle="--", color=color, linewidth=0.75, alpha=0.6, zorder=2)
            ax.text(label_x, yy, name, fontsize=6, color=color, va="center", ha="left")

    ax.set_xticks(np.arange(len(tf_subset)))
    ax.set_xticklabels(tf_subset, rotation=45, ha="right", fontsize=8)
    ax.set_xlim(-0.8, len(tf_subset) + 0.40)
    if len(y_values) > 0:
        ymin = max(0.0, float(np.nanmin(y_values)) - 0.02)
        ymax = min(1.0, float(np.nanmax(y_values)) + 0.02)
        if ymax - ymin < 0.08:
            mid = 0.5 * (ymin + ymax)
            ymin = max(0.0, mid - 0.04)
            ymax = min(1.0, mid + 0.04)
        ax.set_ylim(ymin, ymax)
    ax.set_ylabel("AUPR")
    ax.grid(False)
    sns.despine(ax=ax)


def main() -> None:
    args = parse_args()

    if not os.path.exists(args.input):
        raise RuntimeError(f"Input file not found: {args.input}")

    with open(args.input, "rb") as f:
        payload = pickle.load(f)

    combined_df = payload.get("combined_table")
    if combined_df is None:
        combined_df = payload.get("combined_all")
    if combined_df is None or len(combined_df) == 0:
        raise RuntimeError("Missing or empty combined dataframe in input payload (expected `combined_table` or `combined_all`)")

    pred_cols = resolve_prediction_columns(combined_df)

    for c in ["label_bin", "motif_scan", "atac_seq"] + list(pred_cols.values()):
        if c in combined_df.columns:
            combined_df[c] = pd.to_numeric(combined_df[c], errors="coerce")

    pair_df = compute_pair_metrics(combined_df, pred_cols)
    pair_complete = pair_df.dropna(subset=MODEL_ORDER).copy()
    if len(pair_complete) == 0:
        raise RuntimeError("No complete TF-tissue pairs after filtering requested model predictions")

    long_violin = pair_complete.melt(
        id_vars=["pair_id", "tf", "tissue", "split", "no_skill", "motif_scan", "atac_seq", "n_rows"],
        value_vars=MODEL_ORDER,
        var_name="model",
        value_name="aupr",
    ).dropna(subset=["aupr"])

    global_pr = compute_global_pr(combined_df, pred_cols)

    baseline_pr = {}
    y_global = combined_df["label_bin"].to_numpy(dtype=np.int8)
    for baseline_name, baseline_col in [("Motif", "motif_scan"), ("ATAC", "atac_seq")]:
        s = combined_df[baseline_col].to_numpy(dtype=np.float32)
        mask = np.isfinite(s)
        if int(mask.sum()) == 0:
            continue
        y_sub = y_global[mask]
        s_sub = s[mask]
        if np.unique(y_sub).size < 2:
            continue
        p, r, _ = precision_recall_curve(y_sub, s_sub)
        baseline_pr[baseline_name] = {
            "precision": p,
            "recall": r,
            "aupr": float(average_precision_score(y_sub, s_sub)),
        }

    no_skill_mean = float(pair_complete["no_skill"].mean())
    motif_mean = float(pair_complete["motif_scan"].mean())
    atac_mean = float(pair_complete["atac_seq"].mean())

    tf_rank = pair_complete.groupby("tf")["tf_transformer"].mean().dropna().sort_values()
    tf_order = tf_rank.index.tolist()
    half = int(math.ceil(len(tf_order) / 2))

    sns.set_theme(style="ticks")
    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.labelsize": 11,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 9,
        }
    )

    fig = plt.figure(figsize=(18, 13))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.05, 1.20, 1.20], hspace=0.34, wspace=0.18)
    ax2a = fig.add_subplot(gs[0, 0])
    ax2b = fig.add_subplot(gs[0, 1])
    ax2c_top = fig.add_subplot(gs[1, :])
    ax2c_bottom = fig.add_subplot(gs[2, :])

    violin_palette = [COLORS[m] for m in MODEL_ORDER]
    sns.violinplot(
        data=long_violin,
        x="model",
        y="aupr",
        order=MODEL_ORDER,
        palette=violin_palette,
        inner=None,
        cut=0,
        linewidth=1.0,
        ax=ax2a,
    )
    sns.stripplot(
        data=long_violin,
        x="model",
        y="aupr",
        order=MODEL_ORDER,
        color="black",
        size=1.8,
        alpha=0.18,
        jitter=0.09,
        ax=ax2a,
    )

    means = long_violin.groupby("model")["aupr"].mean().reindex(MODEL_ORDER)
    for i, m in enumerate(MODEL_ORDER):
        mu = means[m]
        if pd.isna(mu):
            continue
        ax2a.hlines(mu, i - 0.24, i + 0.24, colors="black", linestyles="-", linewidth=2.1, zorder=5)
        ax2a.text(i, mu + 0.013, f"{mu:.3f}", ha="center", va="bottom", fontsize=11, fontweight="bold", color="black")

    add_baseline_lines(ax2a, no_skill_mean, motif_mean, atac_mean, -0.5, len(MODEL_ORDER) - 0.5)

    ttest_df = run_selected_paired_ttests(pair_complete)
    x_lookup = {m: i for i, m in enumerate(MODEL_ORDER)}
    y_start = float(np.nanmax(long_violin["aupr"])) + 0.03
    y_step = 0.035
    y_max_needed = y_start
    for i, row in ttest_df.reset_index(drop=True).iterrows():
        a = str(row["model_a"])
        b = str(row["model_b"])
        if a not in x_lookup or b not in x_lookup:
            continue
        y_here = y_start + i * y_step
        y_max_needed = max(y_max_needed, y_here + 0.02)
        add_sig_bracket(ax2a, x_lookup[a], x_lookup[b], y_here, format_p_number(float(row["p_value"])), h=0.010)

    ax2a.set_xticklabels([DISPLAY_NAMES[m] for m in MODEL_ORDER], rotation=24, ha="right")
    ax2a.set_xlabel("")
    ax2a.set_ylabel("AUPR")
    ax2a.set_ylim(0.0, max(0.8, y_max_needed))
    ax2a.text(0.02, 0.97, f"n={len(pair_complete)}", transform=ax2a.transAxes, ha="left", va="top", fontsize=11, fontweight="bold")
    ax2a.grid(False)
    sns.despine(ax=ax2a)
    ax2a.text(-0.13, 1.02, "a", transform=ax2a.transAxes, fontsize=14, fontweight="bold")

    for m in MODEL_ORDER:
        if m not in global_pr:
            continue
        cur = global_pr[m]
        display_name_2b = DISPLAY_NAMES[m].replace("\n", " ")
        ax2b.plot(
            cur["recall"],
            cur["precision"],
            color=COLORS[m],
            linewidth=2.0,
            label=f"{display_name_2b} (AUPR={cur['aupr']:.3f})",
        )
    for baseline_name in ["Motif", "ATAC"]:
        if baseline_name not in baseline_pr:
            continue
        cur = baseline_pr[baseline_name]
        color = BASELINE_COLORS["motif_scan"] if baseline_name == "Motif" else BASELINE_COLORS["atac_seq"]
        ax2b.plot(
            cur["recall"],
            cur["precision"],
            color=color,
            linewidth=1.2,
            linestyle="--",
            alpha=0.7,
            label=f"{baseline_name} (AUPR={cur['aupr']:.3f})",
        )
    ax2b.plot([0, 1], [no_skill_mean, no_skill_mean], "--", color=BASELINE_COLORS["no_skill"], linewidth=1.4)
    ax2b.set_xlim(0, 1)
    ax2b.set_ylim(0, 1)
    ax2b.set_xlabel("Recall")
    ax2b.set_ylabel("Precision")
    ax2b.legend(loc="upper right", frameon=False)
    ax2b.grid(False)
    sns.despine(ax=ax2b)
    ax2b.text(-0.13, 1.02, "b", transform=ax2b.transAxes, fontsize=14, fontweight="bold")

    plot_tf_panel(ax2c_top, pair_complete, tf_order[:half])
    plot_tf_panel(ax2c_bottom, pair_complete, tf_order[half:])
    ax2c_top.text(-0.01, 1.04, "c", transform=ax2c_top.transAxes, fontsize=14, fontweight="bold")

    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=COLORS[m],
            markeredgecolor="white",
            markeredgewidth=0.5,
            markersize=7,
            label=DISPLAY_NAMES[m],
        )
        for m in MODEL_ORDER
    ]
    legend_handles.append(Line2D([0], [0], color="#8a8a8a", linewidth=1.0, label="Same TF-tissue pair"))
    legend_handles.append(Line2D([0], [0], color=BASELINE_COLORS["atac_seq"], linestyle="--", linewidth=1.2, label="Baselines"))
    ax2c_top.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, 1.26), ncol=4, fontsize=9, frameon=False)

    plt.tight_layout()

    os.makedirs(os.path.dirname(args.out_png), exist_ok=True)
    fig.savefig(args.out_png, dpi=300, bbox_inches="tight")
    fig.savefig(args.out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved figure: {args.out_png}")
    print(f"Saved figure: {args.out_pdf}")


if __name__ == "__main__":
    main()
