#!/usr/bin/env python3
"""Central configuration for TFBS prediction scripts."""

from pathlib import Path
import glob

REPO_ROOT = Path(__file__).resolve().parent
DATA_DIR = REPO_ROOT / "data"
TRAINING_DATA_DIR = DATA_DIR / "training_data"
TEST_DATA_DIR = DATA_DIR / "test_data"
FIGURES_DIR = REPO_ROOT / "figures"
SCRIPTS_OUTPUT_DIR = REPO_ROOT / "scripts" / "outputs"

# -----------------------------
# Binaries
# -----------------------------
TRAP_BIN = "/project/ngsvin/bin/TRAP/TRAPv1.04"
MOODS_BIN = "moods-dna.py"
TOBIAS_BIN = "TOBIAS"
BWTOOL_BIN = "bwtool"
BEDTOOLS_BIN = "bedtools"

# -----------------------------
# Reference resources
# -----------------------------
HG38 = str(DATA_DIR / "reference" / "hg38.fa")
BLACKLIST = str(DATA_DIR / "reference" / "GRCh38_unified_blacklist.bed")

ENHANCERS_BED = str(DATA_DIR / "enhancers" / "all_regions.bed")
ENHANCERS_FASTA = str(DATA_DIR / "enhancers" / "all_regions.fa")

TRAP_BG = str(DATA_DIR / "motifs" / "trap_background" / "e24_random30k_bgregions2.GEVparams")
MOTIF_FILE = str(DATA_DIR / "motifs" / "TFP2022_recommended_ekin.fa")
MOTIF_DIR = str(DATA_DIR / "motifs" / "tfp2022_ekinrecommended_pfms")

GENE_TO_TRANSFAC_FILE = str(DATA_DIR / "metadata" / "gene_symbol_to_transfac.txt")
ENCODE_TISSUES_FILE = str(DATA_DIR / "metadata" / "chosen_encode_tissues.txt")
BULK_MRNA_FILE = str(DATA_DIR / "expression" / "all_tissues_meantpm.tsv")
TF_ACTIVITY_FILE = str(DATA_DIR / "tf_links" / "tf_activities_fixed.tsv")
CONS_FILE = str(DATA_DIR / "enhancers" / "all_regions_phastcons.tsv")
REMAP_FILE = str(DATA_DIR / "enhancers" / "all_regions_remap.tsv")
CRUP_FILE = str(DATA_DIR / "enhancers" / "all_regions_crupscores.bed")

# NT embeddings used by feature extraction
EMBEDDINGS_CSV = str(TRAINING_DATA_DIR / "embeddings.csv")
TRAINING_SETS_PK = str(TRAINING_DATA_DIR / "training_sets.pk")
FINAL_ROUND_SETS_PK = str(TEST_DATA_DIR / "final_round_sets.pk")
LEADERBOARD_SETS_PK = str(TEST_DATA_DIR / "leaderboard_sets.pk")
TEST_SETS_PK = FINAL_ROUND_SETS_PK

# Manuscript training/validation sets
MANUSCRIPT_TRAINING_SETS_PK = TRAINING_SETS_PK
MANUSCRIPT_VALIDATION_SETS_PK = str(DATA_DIR / "manuscript" / "validation_sets.pk")
MANUSCRIPT_EMBEDDINGS_CSV = EMBEDDINGS_CSV
PCA_N_COMPONENTS = 50
MANUSCRIPT_EMBEDDINGS_PCA_MODEL_PK = str(DATA_DIR / "models" / "embeddings_pca100_model.pk")
MANUSCRIPT_EMBEDDINGS_PCA100_PK = str(DATA_DIR / "models" / "embeddings_pca100_float32.pk")

# -----------------------------
# Shared artifacts
# -----------------------------
GLOBAL_ATAC_MEAN = str(DATA_DIR / "atac" / "atac_mean_mean.csv")
GLOBAL_ATAC_QN_MIN = str(DATA_DIR / "atac" / "atac_min_quantilenormalized_values.csv")
GLOBAL_ATAC_QN_MAX = str(DATA_DIR / "atac" / "atac_max_quantilenormalized_values.csv")
GLOBAL_ATAC_QN_MEAN = str(DATA_DIR / "atac" / "atac_mean_quantilenormalized_values.csv")
GLOBAL_ATAC_QN_MEAN_TABLE = str(DATA_DIR / "atac" / "atac_qnorm_mean_dream.csv")
GLOBAL_TOBIAS_MEAN = str(DATA_DIR / "tobias" / "tobias_mean_mean.csv")

# -----------------------------
# Model paths
# -----------------------------
MODEL_DIR = str(DATA_DIR / "models")
GENERAL_BASE_MODEL_JSON = f"{MODEL_DIR}/general_models/general_base_model.json"
TF_TUNED_MODEL_GLOB = f"{MODEL_DIR}/tf_models/*_base_model.json"
TF_ONLY_MODELS_PK = f"{MODEL_DIR}/tf_models_nogeneral.pk"
TF_TRANSFORMER_MODEL_DIR = f"{MODEL_DIR}/tf_transformer_models"

ENSEMBLE_MODEL_ROOT = f"{MODEL_DIR}/ensemble"
ENSEMBLE_GENERAL_DIR = f"{ENSEMBLE_MODEL_ROOT}/general_model_tf4_chr2_8"
ENSEMBLE_TF_SUMMARY = f"{ENSEMBLE_MODEL_ROOT}/tf_families_fold_summary.tsv"
ENSEMBLE_TF_ONLY_DIR = f"{ENSEMBLE_MODEL_ROOT}/tf_only"
ENSEMBLE_TF_TUNED_DIR = f"{ENSEMBLE_MODEL_ROOT}/tf_tuned"
ENSEMBLE_TF_TRANSFORMER_DIR = f"{ENSEMBLE_MODEL_ROOT}/tf_transformer"

PCA100_TUNED_MODEL_ROOT = f"{MODEL_DIR}/pca100_tuned"
PCA100_TUNED_GENERAL_DIR = f"{PCA100_TUNED_MODEL_ROOT}/general_pc100_tf4_chr2_8"
PCA100_TUNED_TF_SUMMARY = f"{PCA100_TUNED_MODEL_ROOT}/tf_families_fold_summary.tsv"
PCA100_TUNED_PCA_MODEL_PK = f"{PCA100_TUNED_MODEL_ROOT}/pca_nt_100.pkl"

FIG2_INFERENCE_EXPORT_PK = str(DATA_DIR / "paper_figure_2_prediction_export.pkl")
FIG2_PREDICTION_EXPORT_PKL = FIG2_INFERENCE_EXPORT_PK
GLOBAL_SHAP_TSV = str(SCRIPTS_OUTPUT_DIR / "fig4" / "global_shap_test_leaderboard_all4models.tsv")
GLOBAL_SHAP_PLOT_DIR = str(SCRIPTS_OUTPUT_DIR / "fig4")

GLOBAL_SHAP_OUTPUT_DIR = str(SCRIPTS_OUTPUT_DIR / "fig_sz")
LOCAL_SHAP_OUTPUT_DIR = str(SCRIPTS_OUTPUT_DIR / "fig_st")

# Notebook/script compatibility aliases
GENERAL_MODELS_DIR = f"{MODEL_DIR}/general_models"
TF_MODELS_DIR = f"{MODEL_DIR}/tf_models"
TF_TRANSFORMER_MODELS_DIR = TF_TRANSFORMER_MODEL_DIR
ENSEMBLE_TF_TRANSFORMER_TUNED_DIR = ENSEMBLE_TF_TRANSFORMER_DIR

DREAM_CONTESTANTS_DIR = str(DATA_DIR / "dream_contestants")
PREDICT_CHIPSEQ_WITH_CHIPSEQ_DIR = str(DATA_DIR / "predict_chipseq_with_chipseq")
NEW_TFS_DIR = str(DATA_DIR / "new_tfs")
GENERAL_MODEL_NT = f"{MODEL_DIR}/general_model_nt.json"
GENERAL_MODEL_TRANSFORMER_ALLPAIRS = f"{MODEL_DIR}/general_model_transformer_allpairs.json"
TF_TUNED_TRANSFORMER_MODEL_DIR = f"{MODEL_DIR}/tf_tuned_transformer"
GENERAL_MODEL_LEGACY = f"{MODEL_DIR}/legacy/general_base_model.json"
TF_TUNED_TRANSFORMER_LEGACY_GENERAL_RAWEMB_MODEL_DIR = f"{MODEL_DIR}/tf_tuned_transformer_legacy_general_rawemb"
TF_TUNED_TRANSFORMER_GENERAL_NT_RAWEMB_MODEL_DIR = f"{MODEL_DIR}/tf_tuned_transformer_general_nt_rawemb"

# Feature extraction outputs used for training
FEATURE_SET_GLOB = str(DATA_DIR / "new_tfs" / "*" / "training_set.pk")

# Leaderboard sets
LEADERBOARD_SETS_PATH = LEADERBOARD_SETS_PK
LEADERBOARD_MODEL_DIR = f"{MODEL_DIR}/leaderboard_tf_only"

# -----------------------------
# Features
# -----------------------------
EMB_COLS = [f"emb_{i}" for i in range(1024)]
GEN_FEATURES = [
    "maxpwm", "cons", "crup", "crup_mean", "crup_delta", "remap", "tf_exp", "tf_activity",
    "tfact_crupcor_coef", "tfact_crupcor_pval", "trap", "atac_min", "atac_max", "atac_mean", "atac_mean_mean",
    "atac_delta_min", "atac_delta_max", "atac_delta_mean", "tobias_avg", "delta_tobias_avg", "tobias_mean_mean", "tobias_count",
    "cot_hits_0", "cot_hits_1", "cot_hits_2", "cot_hits_3", "cot_hits_4", "cot_hits_5",
    "cot_maxpwm_0", "cot_maxpwm_1", "cot_maxpwm_2", "cot_maxpwm_3", "cot_maxpwm_4", "cot_maxpwm_5",
]


def get_feature_training_files():
    return sorted(glob.glob(FEATURE_SET_GLOB))


def ensure_model_dirs():
    Path(MODEL_DIR).mkdir(parents=True, exist_ok=True)
    Path(LEADERBOARD_MODEL_DIR).mkdir(parents=True, exist_ok=True)
    Path(TF_TUNED_TRANSFORMER_MODEL_DIR).mkdir(parents=True, exist_ok=True)
    Path(TF_TUNED_TRANSFORMER_LEGACY_GENERAL_RAWEMB_MODEL_DIR).mkdir(parents=True, exist_ok=True)
    Path(TF_TUNED_TRANSFORMER_GENERAL_NT_RAWEMB_MODEL_DIR).mkdir(parents=True, exist_ok=True)
