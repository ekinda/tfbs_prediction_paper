#!/usr/bin/env python3
"""Fetch ENCODE QC metrics for ENCODE accessions and merge into all test/leaderboard TF-tissue pairs.

Key features:
- Reads accession input from the revision ENCODE ID table.
- Programmatic ENCODE API traversal from file accession -> BAM ancestors -> derived IDR peak files.
- Extracts FRIP and peak_count from peak-file quality metrics when available.
- Uses multiprocessing for faster API collection.
- Produces one dataframe covering all test + leaderboard TF-tissue pairs with quality metrics.

Example:
    python 13_encode_quality_metrics_parallel.py --workers 8
"""

from __future__ import annotations

import argparse
import os
import pickle
import re
import time
from multiprocessing import Pool
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import requests
from config import REPO_ROOT, FINAL_ROUND_SETS_PK, LEADERBOARD_SETS_PK


ENCODE_BASE = "https://www.encodeproject.org"
OUTDIR = str(REPO_ROOT / "analysis_outputs" / "reviewer_encode_quality_fixed")
ENCODE_IDS_CSV = str(REPO_ROOT.parent / "revision" / "encode_chipseq_ids_fixed.csv")
os.makedirs(OUTDIR, exist_ok=True)


def norm_tissue(x: Any) -> str:
    s = str(x).strip().lower().replace("_", "-").replace(" ", "-")
    s = re.sub(r"[^a-z0-9-]", "", s)
    aliases = {
        "ipsc": "induced-pluripotent-stem-cell",
        "induced-pluripotent-stem-cell": "induced-pluripotent-stem-cell",
        "h1": "h1",
        "h1-hesc": "h1",
        "sk-n-sh": "sk-n-sh",
        "wtc11": "wtc11",
        "pc-3": "pc-3",
    }
    return aliases.get(s, s)


def norm_tf(x: Any) -> str:
    return str(x).strip().upper()


def encode_get_json(url: str, retries: int = 4, sleep_s: float = 1.2, timeout: int = 35) -> Dict[str, Any]:
    headers = {"accept": "application/json"}
    for i in range(retries):
        try:
            r = requests.get(url, headers=headers, timeout=timeout)
            if r.status_code == 200:
                return r.json()
            if r.status_code in (429, 500, 502, 503, 504):
                time.sleep(sleep_s * (i + 1))
                continue
            return {}
        except Exception:
            time.sleep(sleep_s * (i + 1))
    return {}


def flatten_numeric(obj: Any, prefix: str = "") -> Dict[str, float]:
    out: Dict[str, float] = {}
    if isinstance(obj, dict):
        for k, v in obj.items():
            key = f"{prefix}.{k}" if prefix else str(k)
            out.update(flatten_numeric(v, key))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            key = f"{prefix}[{i}]" if prefix else f"[{i}]"
            out.update(flatten_numeric(v, key))
    else:
        if isinstance(obj, (int, float)) and np.isfinite(obj):
            out[prefix] = float(obj)
    return out


def pick_metric(flat_map: Dict[str, float], patterns: Iterable[str]) -> float:
    pats = [p.lower() for p in patterns]
    candidates: List[Tuple[str, float]] = []
    for k, v in flat_map.items():
        lk = k.lower()
        if any(p in lk for p in pats):
            candidates.append((k, v))
    if not candidates:
        return np.nan
    candidates.sort(key=lambda kv: (len(kv[0]), kv[0].count(".")))
    return float(candidates[0][1])


def pick_metric_all_terms(flat_map: Dict[str, float], terms: Iterable[str]) -> float:
    terms_l = [t.lower() for t in terms]
    candidates: List[Tuple[str, float]] = []
    for k, v in flat_map.items():
        lk = k.lower()
        if all(t in lk for t in terms_l):
            candidates.append((k, v))
    if not candidates:
        return np.nan
    candidates.sort(key=lambda kv: (len(kv[0]), kv[0].count(".")))
    return float(candidates[0][1])


def _norm_metric_token(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def pick_metric_keys(flat_map: Dict[str, float], keys: Iterable[str]) -> float:
    """Pick by explicit metric-key names (e.g., auc, syn_auc, idr_cutoff).

    This is more reliable than substring search for ENCODE quality metric JSON,
    where many numeric fields coexist in the same object.
    """
    wanted = {_norm_metric_token(k) for k in keys}
    candidates: List[Tuple[str, float]] = []
    for full_key, value in flat_map.items():
        parts = str(full_key).split(".")
        for part in reversed(parts):
            part_clean = re.sub(r"\[[^\]]*\]", "", part)
            if _norm_metric_token(part_clean) in wanted:
                candidates.append((full_key, value))
                break
    if not candidates:
        return np.nan
    candidates.sort(key=lambda kv: (len(kv[0]), kv[0].count(".")))
    return float(candidates[0][1])


def _extract_accession(ref: Any) -> Optional[str]:
    if ref is None:
        return None
    if isinstance(ref, dict):
        ref = ref.get("accession") or ref.get("@id") or ref.get("uuid")
    s = str(ref).strip()
    m = re.search(r"(ENCFF[0-9A-Z]{6})", s)
    return m.group(1) if m else (s if s.startswith("ENCFF") else None)


def _get_file_json(acc: str) -> Dict[str, Any]:
    return encode_get_json(f"{ENCODE_BASE}/files/{acc}/?format=json") if acc else {}


def _extract_experiment_accession(ref: Any) -> Optional[str]:
    if ref is None:
        return None
    s = str(ref)
    m = re.search(r"(ENCSR[0-9A-Z]{6})", s)
    return m.group(1) if m else None


def _get_experiment_json(exp_acc: str) -> Dict[str, Any]:
    return encode_get_json(f"{ENCODE_BASE}/experiments/{exp_acc}/?format=json") if exp_acc else {}


def _is_control_experiment(exp_json: Dict[str, Any]) -> bool:
    assay = str(exp_json.get("assay_title", "")).lower()
    control_type = str(exp_json.get("control_type", "")).lower()
    desc = str(exp_json.get("description", "")).lower()
    target = exp_json.get("target", {}) if isinstance(exp_json.get("target", {}), dict) else {}
    target_label = str(target.get("label", "")).lower()
    text = " | ".join([assay, control_type, desc, target_label])
    return ("control chip-seq" in text) or ("igg" in text)


def _replicate_str(file_json: Dict[str, Any]) -> str:
    b = file_json.get("biological_replicates", [])
    t = file_json.get("technical_replicates", [])
    btxt = ",".join(map(str, sorted(b))) if isinstance(b, list) else str(b)
    ttxt = ",".join(map(str, sorted(t))) if isinstance(t, list) else str(t)
    return f"bio:{btxt}|tech:{ttxt}"


def _iter_derived_from_accessions(file_json: Dict[str, Any]) -> List[str]:
    out: List[str] = []
    for ref in file_json.get("derived_from", []):
        acc = _extract_accession(ref)
        if acc:
            out.append(acc)
    return out


def _collect_flat_from_qm(file_json: Dict[str, Any], prefix: str) -> Dict[str, float]:
    flat = flatten_numeric(file_json, prefix=prefix)
    qms = file_json.get("quality_metrics", [])
    for i, qm in enumerate(qms):
        if isinstance(qm, dict) and "@id" in qm:
            qurl = qm["@id"]
            if qurl.startswith("/"):
                qurl = ENCODE_BASE + qurl
            if "?format=json" not in qurl:
                qurl += "&format=json" if "?" in qurl else "?format=json"
            qmj = encode_get_json(qurl)
            if qmj:
                flat.update(flatten_numeric(qmj, prefix=f"{prefix}.qm[{i}]"))
        elif isinstance(qm, dict):
            flat.update(flatten_numeric(qm, prefix=f"{prefix}.qm_inline[{i}]"))
    return flat


def _search_files_derived_from(parent_acc: str) -> List[Dict[str, Any]]:
    # reverse lookup: children that derive from this file
    url = (
        f"{ENCODE_BASE}/search/?type=File"
        f"&derived_from=/files/{parent_acc}/"
        f"&limit=all"
        f"&field=accession&field=file_format&field=file_format_type&field=output_type&field=status&field=assembly&field=dataset&field=quality_metrics"
        f"&field=biological_replicates&field=technical_replicates"
        f"&format=json"
    )
    sj = encode_get_json(url)
    graph = sj.get("@graph", []) if isinstance(sj, dict) else []
    return [g for g in graph if isinstance(g, dict) and g.get("accession")]


def _is_idr_peak_file(file_json: Dict[str, Any]) -> bool:
    ff = str(file_json.get("file_format", "")).lower()
    ot = str(file_json.get("output_type", "")).lower()
    if ff != "bed":
        return False
    if "conservative" in ot:
        return False
    return "idr" in ot and "peak" in ot and "threshold" in ot


def _is_any_peak_file(file_json: Dict[str, Any]) -> bool:
    ff = str(file_json.get("file_format", "")).lower()
    ot = str(file_json.get("output_type", "")).lower()
    return ff in {"bed", "bigbed", "bed.gz"} and "peak" in ot


def _fetch_encode_metrics_one(task: Tuple[str, str]) -> Dict[str, Any]:
    file_acc, file_type = task
    file_type_l = str(file_type).lower()
    is_peak_root = "peak" in file_type_l
    is_bigwig_root = "bigwig" in file_type_l

    root = _get_file_json(file_acc)
    if not root:
        return {"file_id": file_acc}

    row: Dict[str, Any] = {
        "file_id": file_acc,
        "file_accession": root.get("accession", file_acc),
        "file_format": root.get("file_format"),
        "output_type": root.get("output_type"),
        "status": root.get("status"),
        "assembly": root.get("assembly"),
        "dataset": root.get("dataset"),
    }

    # backward traversal: root -> ancestors (find BAMs)
    queue: List[Tuple[str, int]] = [(file_acc, 0)]
    seen: set[str] = set()
    traversed: List[Dict[str, Any]] = []
    max_depth = 4

    while queue:
        acc, depth = queue.pop(0)
        if acc in seen:
            continue
        seen.add(acc)

        fj = _get_file_json(acc)
        if not fj:
            continue
        traversed.append(fj)

        if depth < max_depth:
            for parent in _iter_derived_from_accessions(fj):
                if parent not in seen:
                    queue.append((parent, depth + 1))

    bam_files = [f for f in traversed if str(f.get("file_format", "")).lower() == "bam"]

    # Exclude IgG/control BAMs using experiment metadata
    control_bams: List[Dict[str, Any]] = []
    non_control_bams: List[Dict[str, Any]] = []
    for b in bam_files:
        exp_acc = _extract_experiment_accession(b.get("dataset"))
        expj = _get_experiment_json(exp_acc) if exp_acc else {}
        if _is_control_experiment(expj):
            control_bams.append(b)
        else:
            non_control_bams.append(b)

    # forward traversal from non-control BAMs to derived peak files (primary path for bigWig roots)
    derived_peak_files: List[Dict[str, Any]] = []
    peak_seen: set[str] = set()
    bam_to_peak_candidates: Dict[str, List[Dict[str, Any]]] = {}

    if is_bigwig_root or not is_peak_root:
        for b in non_control_bams:
            bacc = b.get("accession")
            if not bacc:
                continue
            children = _search_files_derived_from(bacc)
            bam_to_peak_candidates[bacc] = []
            for child in children:
                cacc = child.get("accession")
                if not cacc:
                    continue
                if _is_idr_peak_file(child):
                    bam_to_peak_candidates[bacc].append(child)
                    if cacc not in peak_seen:
                        derived_peak_files.append(child)
                        peak_seen.add(cacc)

    # fallback to any peaks only for peak-root rows (never for bigWig roots)
    if (not is_bigwig_root) and is_peak_root and not derived_peak_files:
        for b in non_control_bams:
            bacc = b.get("accession")
            if not bacc:
                continue
            children = _search_files_derived_from(bacc)
            bam_to_peak_candidates.setdefault(bacc, [])
            for child in children:
                cacc = child.get("accession")
                if not cacc:
                    continue
                if _is_any_peak_file(child):
                    bam_to_peak_candidates[bacc].append(child)
                    if cacc not in peak_seen:
                        derived_peak_files.append(child)
                        peak_seen.add(cacc)

    # Choose peak files: prefer intersection across BAM replicates (joint file), else one per BAM then mean
    chosen_peak_files: List[Dict[str, Any]] = []

    # First preference: any IDR-thresholded file explicitly representing multiple biological replicates.
    # This is usually the recommended joint/all-replicates peak file.
    preferred_joint = [
        f for f in derived_peak_files
        if isinstance(f.get("biological_replicates", []), list) and len(f.get("biological_replicates", [])) >= 2
    ]
    if preferred_joint:
        preferred_joint.sort(
            key=lambda f: (
                -len(f.get("biological_replicates", [])),
                str(f.get("accession", "")),
            )
        )
        chosen_peak_files = [preferred_joint[0]]

    if bam_to_peak_candidates:
        per_bam_sets = []
        for bacc, lst in bam_to_peak_candidates.items():
            if lst:
                per_bam_sets.append(set([str(x.get("accession")) for x in lst if x.get("accession")]))

        if (not chosen_peak_files) and len(per_bam_sets) >= 2:
            inter = set.intersection(*per_bam_sets)
            if inter:
                inter_files = [f for f in derived_peak_files if str(f.get("accession")) in inter]
                # prefer files explicitly marked with multi-replicate support
                inter_files.sort(
                    key=lambda f: (
                        -len(f.get("biological_replicates", []) if isinstance(f.get("biological_replicates", []), list) else []),
                        str(f.get("accession", "")),
                    )
                )
                chosen_peak_files = [inter_files[0]]

    if not chosen_peak_files:
        # fallback: choose one best candidate per BAM replicate
        for bacc, lst in bam_to_peak_candidates.items():
            if not lst:
                continue
            lst_sorted = sorted(
                lst,
                key=lambda f: (
                    -len(f.get("biological_replicates", []) if isinstance(f.get("biological_replicates", []), list) else []),
                    str(f.get("accession", "")),
                ),
            )
            chosen_peak_files.append(lst_sorted[0])

    # deduplicate chosen peaks
    chosen_seen: set[str] = set()
    chosen_peak_files = [
        f for f in chosen_peak_files
        if f.get("accession") and not (f.get("accession") in chosen_seen or chosen_seen.add(f.get("accession")))
    ]

    # for rows explicitly marked as peaks, include the root peak file itself for FRIP/peak_count
    root_as_peak_files: List[Dict[str, Any]] = []
    if is_peak_root and _is_any_peak_file(root):
        root_as_peak_files.append(root)

    bam_metric_files = non_control_bams if len(non_control_bams) > 0 else bam_files
    peak_metric_files = chosen_peak_files if len(chosen_peak_files) > 0 else derived_peak_files
    source_files = bam_metric_files + peak_metric_files + root_as_peak_files
    if not source_files:
        source_files = traversed

    bam_flat: Dict[str, float] = {}
    for i, fj in enumerate(bam_metric_files):
        acc = fj.get("accession", f"bam_{i}")
        bam_flat.update(_collect_flat_from_qm(fj, prefix=f"bam[{i}].{acc}"))

    peak_flat: Dict[str, float] = {}
    peak_sources = peak_metric_files + root_as_peak_files
    for i, fj in enumerate(peak_sources):
        acc = fj.get("accession", f"peak_{i}")
        peak_flat.update(_collect_flat_from_qm(fj, prefix=f"peak[{i}].{acc}"))

    flat: Dict[str, float] = {}
    for i, fj in enumerate(source_files):
        acc = fj.get("accession", f"file_{i}")
        flat.update(_collect_flat_from_qm(fj, prefix=f"source[{i}].{acc}"))

    row["n_files_traversed"] = len(traversed)
    row["n_bam_found"] = len(bam_files)
    row["n_bam_control_ignored"] = len(control_bams)
    row["n_bam_used"] = len(bam_metric_files)
    row["n_candidate_peak_files"] = len(derived_peak_files)
    row["n_derived_peak_files"] = len(derived_peak_files)
    row["n_root_peak_files"] = len(root_as_peak_files)
    row["derived_bam_accessions"] = ";".join(sorted([str(f.get("accession")) for f in bam_files if f.get("accession")]))
    row["bam_used_accessions"] = ";".join(sorted([str(f.get("accession")) for f in bam_metric_files if f.get("accession")]))
    row["bam_control_accessions"] = ";".join(sorted([str(f.get("accession")) for f in control_bams if f.get("accession")]))
    row["derived_peak_accessions"] = ";".join(sorted([str(f.get("accession")) for f in derived_peak_files if f.get("accession")]))
    row["chosen_peak_accessions"] = ";".join(sorted([str(f.get("accession")) for f in chosen_peak_files if f.get("accession")]))
    row["chosen_peak_replicates"] = ";".join([
        f"{f.get('accession')}|{_replicate_str(f)}" for f in chosen_peak_files if f.get("accession")
    ])

    # canonical metrics
    peak_src = peak_flat if peak_flat else flat
    row["frip"] = pick_metric_keys(peak_src, ["frip", "fraction_of_reads_in_peaks"])
    row["reproducible_peaks"] = pick_metric_keys(peak_src, ["reproducible_peaks", "reproducible_peak"])
    row["peak_count"] = pick_metric_keys(peak_src, ["peak_count", "number_of_peaks", "npeaks"])
    if np.isnan(row["peak_count"]):
        row["peak_count"] = row["reproducible_peaks"]
    row["average_peak_size"] = pick_metric_keys(peak_src, ["average_peak_size", "mean"])
    row["peak_min_size"] = pick_metric_keys(peak_src, ["min_size"])
    row["peak_25_pct"] = pick_metric_keys(peak_src, ["25_pct"])
    row["peak_50_pct"] = pick_metric_keys(peak_src, ["50_pct"])
    row["peak_75_pct"] = pick_metric_keys(peak_src, ["75_pct"])
    row["peak_max_size"] = pick_metric_keys(peak_src, ["max_size"])
    row["idr_cutoff"] = pick_metric_keys(peak_src, ["idr_cutoff"])
    row["idr_rescue_ratio"] = pick_metric_keys(peak_src, ["rescue_ratio"])
    row["idr_self_consistency_ratio"] = pick_metric_keys(peak_src, ["self_consistency_ratio"])

    row["nrf"] = pick_metric(bam_flat if bam_flat else flat, ["nrf", "non_redundant_fraction"])
    row["pbc1"] = pick_metric(bam_flat if bam_flat else flat, ["pbc1", "pbc_1"])
    row["pbc2"] = pick_metric(bam_flat if bam_flat else flat, ["pbc2", "pbc_2"])
    row["nsc"] = pick_metric(bam_flat if bam_flat else flat, ["nsc", "normalized_strand_cross_correlation"])
    row["rsc"] = pick_metric(bam_flat if bam_flat else flat, ["rsc", "relative_strand_cross_correlation"])
    row["mapped_reads"] = pick_metric(bam_flat if bam_flat else flat, ["mapped_read", "mapped_reads"])
    row["usable_fragments"] = pick_metric(bam_flat if bam_flat else flat, ["usable_fragment", "usable_fragments"])
    deep_src = bam_flat if bam_flat else flat
    row["deeptools_auc"] = pick_metric_keys(deep_src, ["auc"])
    row["deeptools_synthetic_auc"] = pick_metric_keys(deep_src, ["syn_auc", "synthetic_auc"])
    row["enrichment_x_intercept"] = pick_metric_keys(deep_src, ["x_intercept"])
    row["enrichment_syn_x_intercept"] = pick_metric_keys(deep_src, ["syn_x_intercept", "synthetic_x_intercept"])
    row["enrichment_elbow_pt"] = pick_metric_keys(deep_src, ["elbow_pt"])
    row["enrichment_syn_elbow_pt"] = pick_metric_keys(deep_src, ["syn_elbow_pt", "synthetic_elbow_pt"])
    row["enrichment_jsd"] = pick_metric_keys(deep_src, ["jsd"])
    row["enrichment_syn_jsd"] = pick_metric_keys(deep_src, ["syn_jsd", "synthetic_jsd"])
    row["pct_genome_enrich"] = pick_metric_keys(deep_src, ["pct_genome_enrich"])
    row["diff_enrich"] = pick_metric_keys(deep_src, ["diff_enrich"])
    row["ch_div"] = pick_metric_keys(deep_src, ["ch_div"])
    row["estimated_fragment_len"] = pick_metric_keys(deep_src, ["estimated_fragment_len"])
    row["corr_estimated_fragment_len"] = pick_metric_keys(deep_src, ["corr_estimated_fragment_len"])
    row["phantom_peak"] = pick_metric_keys(deep_src, ["phantom_peak"])
    row["corr_phantom_peak"] = pick_metric_keys(deep_src, ["corr_phantom_peak"])

    return row


def load_sup2() -> pd.DataFrame:
    sup2 = pd.read_csv(ENCODE_IDS_CSV)
    rename = {
        "TF": "tf",
        "Cell Type": "cell_type",
        "File ID": "file_id",
        "File Type": "file_type",
    }
    sup2 = sup2.rename(columns=rename)
    required = ["tf", "cell_type", "file_id", "file_type"]
    missing = [c for c in required if c not in sup2.columns]
    if missing:
        raise ValueError(f"Missing required columns in {ENCODE_IDS_CSV}: {missing}")
    sup2["tf_norm"] = sup2["tf"].map(norm_tf)
    sup2["tissue_norm"] = sup2["cell_type"].map(norm_tissue)
    sup2["pair_key"] = sup2["tf_norm"] + "||" + sup2["tissue_norm"]
    return sup2


def load_pair_table() -> pd.DataFrame:
    with open(FINAL_ROUND_SETS_PK, "rb") as f:
        test_sets = pickle.load(f)
    with open(LEADERBOARD_SETS_PK, "rb") as f:
        leaderboard_sets = pickle.load(f)

    rows: List[Dict[str, str]] = []

    for tf, df in test_sets.items():
        tissues = df["tissue"].astype(str).dropna().unique()
        for tissue in tissues:
            rows.append({"split": "test", "tf": tf, "tissue": tissue})

    for tf, df in leaderboard_sets.items():
        tissues = df["tissue"].astype(str).dropna().unique()
        for tissue in tissues:
            rows.append({"split": "leaderboard", "tf": tf, "tissue": tissue})

    pairs = pd.DataFrame(rows).drop_duplicates()
    pairs["tf_norm"] = pairs["tf"].map(norm_tf)
    pairs["tissue_norm"] = pairs["tissue"].map(norm_tissue)
    pairs["pair_key"] = pairs["tf_norm"] + "||" + pairs["tissue_norm"]
    return pairs


def build_pair_level_metrics(encode_metrics: pd.DataFrame, sup2: pd.DataFrame) -> pd.DataFrame:
    # attach pair mapping
    enc = sup2[["pair_key", "file_id", "file_type"]].merge(encode_metrics, on="file_id", how="left")

    metric_cols = [
        "frip", "peak_count", "reproducible_peaks", "average_peak_size",
        "peak_min_size", "peak_25_pct", "peak_50_pct", "peak_75_pct", "peak_max_size",
        "idr_cutoff", "idr_rescue_ratio", "idr_self_consistency_ratio",
        "nrf", "pbc1", "pbc2", "nsc", "rsc", "mapped_reads", "usable_fragments",
        "deeptools_auc", "deeptools_synthetic_auc",
        "enrichment_x_intercept", "enrichment_syn_x_intercept",
        "enrichment_elbow_pt", "enrichment_syn_elbow_pt",
        "enrichment_jsd", "enrichment_syn_jsd",
        "pct_genome_enrich", "diff_enrich", "ch_div",
        "estimated_fragment_len", "corr_estimated_fragment_len",
        "phantom_peak", "corr_phantom_peak",
        "n_files_traversed", "n_bam_found", "n_bam_control_ignored", "n_bam_used",
        "n_candidate_peak_files", "n_derived_peak_files", "n_root_peak_files",
    ]

    def _first_non_null(series: pd.Series):
        s = series.dropna()
        return s.iloc[0] if len(s) else np.nan

    agg_map = {c: "mean" for c in metric_cols}
    agg_map.update({
        "file_id": lambda s: ";".join(sorted(set(map(str, s.dropna())))),
        "file_type": lambda s: ";".join(sorted(set(map(str, s.dropna())))),
        "bam_used_accessions": _first_non_null,
        "bam_control_accessions": _first_non_null,
        "derived_bam_accessions": _first_non_null,
        "derived_peak_accessions": _first_non_null,
        "chosen_peak_accessions": _first_non_null,
        "chosen_peak_replicates": _first_non_null,
    })

    for c in [
        "bam_used_accessions", "bam_control_accessions",
        "derived_bam_accessions", "derived_peak_accessions",
        "chosen_peak_accessions", "chosen_peak_replicates",
    ]:
        if c not in enc.columns:
            enc[c] = np.nan

    pair_metrics = enc.groupby("pair_key", as_index=False).agg(agg_map)
    pair_metrics = pair_metrics.rename(columns={
        "file_id": "sup2_file_ids",
        "file_type": "sup2_file_types",
    })
    return pair_metrics


def run(workers: int, force_refresh: bool, sample_n: Optional[int]) -> None:
    sup2 = load_sup2()
    pairs = load_pair_table()

    tasks = list(
        sup2[["file_id", "file_type"]]
        .dropna(subset=["file_id"])
        .astype(str)
        .itertuples(index=False, name=None)
    )
    if sample_n is not None and sample_n > 0:
        tasks = tasks[:sample_n]

    cache_raw = os.path.join(OUTDIR, "encode_quality_metrics_raw.tsv")
    if os.path.exists(cache_raw) and not force_refresh and sample_n is None:
        encode_metrics = pd.read_csv(cache_raw, sep="\t")
    else:
        if workers <= 1:
            rows = [_fetch_encode_metrics_one(task) for task in tasks]
        else:
            with Pool(processes=workers) as pool:
                rows = pool.map(_fetch_encode_metrics_one, tasks)
        encode_metrics = pd.DataFrame(rows)
        if sample_n is None:
            encode_metrics.to_csv(cache_raw, sep="\t", index=False)

    pair_metrics = build_pair_level_metrics(encode_metrics, sup2)
    final_df = pairs.merge(pair_metrics, on="pair_key", how="left")

    # Save outputs
    out_pairs = os.path.join(OUTDIR, "encode_quality_metrics_all_test_leaderboard_pairs.tsv")
    out_raw = os.path.join(OUTDIR, "encode_quality_metrics_sup2_files_current_run.tsv")

    final_df.to_csv(out_pairs, sep="\t", index=False)
    encode_metrics.to_csv(out_raw, sep="\t", index=False)

    metric_cols = [
        "frip", "peak_count", "reproducible_peaks", "average_peak_size",
        "peak_min_size", "peak_25_pct", "peak_50_pct", "peak_75_pct", "peak_max_size",
        "idr_cutoff", "idr_rescue_ratio", "idr_self_consistency_ratio",
        "nrf", "pbc1", "pbc2", "nsc", "rsc", "mapped_reads", "usable_fragments",
        "deeptools_auc", "deeptools_synthetic_auc",
        "enrichment_x_intercept", "enrichment_syn_x_intercept",
        "enrichment_elbow_pt", "enrichment_syn_elbow_pt",
        "enrichment_jsd", "enrichment_syn_jsd",
        "pct_genome_enrich", "diff_enrich", "ch_div",
        "estimated_fragment_len", "corr_estimated_fragment_len",
        "phantom_peak", "corr_phantom_peak",
    ]
    print("Saved:")
    print(" -", out_raw)
    print(" -", out_pairs)
    print("\nRaw ENCODE metric non-NaN counts:")
    existing = [c for c in metric_cols if c in encode_metrics.columns]
    print(encode_metrics[existing].notna().sum())
    print("\nFinal pair table shape:", final_df.shape)
    print("Pairs with any quality metric:", int(final_df[existing].notna().any(axis=1).sum()) if existing else 0)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Fetch ENCODE QC metrics in parallel and merge into all test/leaderboard TF-tissue pairs")
    p.add_argument("--workers", type=int, default=8, help="Number of processes for ENCODE API calls (e.g. 4 or 8)")
    p.add_argument("--force-refresh", action="store_true", help="Ignore cached raw metrics and refetch")
    p.add_argument("--sample-n", type=int, default=None, help="Optional debug mode: only process first N accessions")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(workers=args.workers, force_refresh=args.force_refresh, sample_n=args.sample_n)
