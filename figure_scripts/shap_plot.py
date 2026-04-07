#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
import config as cfg

DEFAULT_TSV = cfg.GLOBAL_SHAP_TSV
DEFAULT_OUTDIR = cfg.GLOBAL_SHAP_PLOT_DIR
DEFAULT_TOP_N = 15
DEFAULT_WIDTH = 6.2
DEFAULT_HEIGHT = 5.0


def safe_name(model_name: str) -> str:
    return (
        str(model_name)
        .replace("+", "plus")
        .replace(" ", "_")
        .replace("/", "_")
        .replace("(", "")
        .replace(")", "")
        .lower()
    )


def plot_model(sub: pd.DataFrame, model_name: str, outdir: str, top_n: int, width: float, height: float) -> None:
    if len(sub) == 0:
        return

    top = (
        sub.nlargest(top_n, "mean_abs_shap")
        .sort_values("mean_abs_shap", ascending=True)
        .reset_index(drop=True)
    )

    fig, ax = plt.subplots(figsize=(width, height))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.grid(False)
    y = np.arange(len(top))
    ax.barh(y, top["mean_abs_shap"], color="#2563eb", alpha=0.92)
    ax.set_yticks(y)
    ax.set_yticklabels(top["feature"], fontsize=12)
    ax.set_xlabel("Mean contribution (SHAP)", fontsize=12)
    ax.set_title(f"Top {top_n} absolute SHAP | {model_name}", fontsize=13)
    ax.tick_params(axis="x", labelsize=11)

    plt.tight_layout()

    stem = f"fig4_shap_plotter_top{top_n}_{safe_name(model_name)}"
    pdf_path = os.path.join(outdir, f"{stem}.pdf")
    plt.savefig(pdf_path, dpi=300, bbox_inches="tight")
    plt.close()

    print("Saved:", pdf_path)


def run(tsv_path: str, outdir: str, top_n: int, width: float, height: float) -> None:
    os.makedirs(outdir, exist_ok=True)
    sns.set(style="white", context="talk")

    df = pd.read_csv(tsv_path, sep="\t")
    needed = {"model", "feature", "mean_abs_shap"}
    if not needed.issubset(set(df.columns)):
        missing = sorted(list(needed.difference(set(df.columns))))
        raise ValueError(f"Missing required columns in TSV: {missing}")

    models = [m for m in sorted(df["model"].dropna().astype(str).unique().tolist())]
    if len(models) == 0:
        raise ValueError("No models found in TSV")

    for model_name in models:
        sub = df[df["model"].astype(str) == model_name].copy()
        plot_model(sub, model_name, outdir, top_n, width, height)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot compact per-model top-N absolute SHAP barplots from fig4 TSV")
    p.add_argument("--tsv", type=str, default=DEFAULT_TSV, help="Input TSV with columns model, feature, mean_abs_shap")
    p.add_argument("--outdir", type=str, default=DEFAULT_OUTDIR, help="Output directory for PDF figures")
    p.add_argument("--top-n", type=int, default=DEFAULT_TOP_N, help="Top N features per model")
    p.add_argument("--width", type=float, default=DEFAULT_WIDTH, help="Figure width in inches")
    p.add_argument("--height", type=float, default=DEFAULT_HEIGHT, help="Figure height in inches")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(args.tsv, args.outdir, args.top_n, args.width, args.height)
