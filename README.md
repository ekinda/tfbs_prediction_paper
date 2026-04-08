# TFBS prediction paper repository

Code and data layout for reproducing model outputs and figures from:
"Transfer learning and DNA language models enhance transcription factor binding predictions".

## Data layout (required)

Download the [Zenodo tarball](https://zenodo.org/records/19457863), then unpack into [data](data) so that these files exist:

- [data/training_data/training_sets.pk](data/training_data/training_sets.pk)
- [data/training_data/embeddings.csv](data/training_data/embeddings.csv)
- [data/test_data/final_round_sets.pk](data/test_data/final_round_sets.pk)
- [data/test_data/leaderboard_sets.pk](data/test_data/leaderboard_sets.pk)
- model files under [data/models](data/models)

## Environment setup

Run from repository root.

1) Create and activate virtual environment

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip setuptools wheel
```

2) Install Python dependencies

```bash
pip install -r requirements.txt
```

Optional variants:

```bash
# Equivalent explicit form
pip install --requirement requirements.txt

# Editable local install only if you add package metadata (pyproject.toml/setup.py)
pip install -e . -r requirements.txt
```

3) External command-line tools (install yourself)

The scripts expect these tools and versions (or compatible versions):

- TRAP v1.04 (`TRAP_BIN`)
- MOODS `moods-dna.py` 1.9.4
- TOBIAS 0.17.1
- bwtool 1.0
- bedtools 2.31.1

Tool paths and file locations are centralized in [config.py](config.py).

## Training tutorial (scripts 1-7)

Run scripts from repository root with `python scripts/<script>.py ...`.

### Important optional/required split

Scripts 1-4 are **optional** when using the downloadable data pack.

The data pack already includes full feature matrices and embeddings required for model training, so most users should **start directly from script 5**.

- Optional preprocessing: scripts 1-4
- Core model training: scripts 5-7

### Script 1: NT embeddings (optional)

File: [scripts/1_nucleotide_transformer_embeddings.py](scripts/1_nucleotide_transformer_embeddings.py)

Purpose:
- Runs the InstaDeep Nucleotide Transformer model over enhancer FASTA.
- Produces per-enhancer embedding vectors used by downstream feature extraction.

Inputs from [config.py](config.py):
- `ENHANCERS_FASTA`

Output:
- `EMBEDDINGS_CSV` (typically [data/training_data/embeddings.csv](data/training_data/embeddings.csv))

Command:

```bash
python scripts/1_nucleotide_transformer_embeddings.py
```

### Script 2: ATAC + TOBIAS processing (optional)

File: [scripts/2_atac_tobias_processing.py](scripts/2_atac_tobias_processing.py)

Purpose:
- Runs TOBIAS `ATACorrect` and `ScoreBigwig` on ATAC BAM + peak set.
- Summarizes footprint and ATAC tracks over enhancer regions.

Outputs:
- Per-sample TOBIAS footprint bigWig and summaries in [data/atac](data/atac)
- Tissue-level ATAC summary table in [data/atac](data/atac)

Command template:

```bash
python scripts/2_atac_tobias_processing.py \
	<tissue> <placeholder_tf_arg> <atac_bigwig> <atac_bam> <atac_peak_bed>
```

Note: argument 2 is currently unused in this script and kept for CLI compatibility.

### Script 3: global ATAC/TOBIAS means (optional)

File: [scripts/3_atac_tobias_means.py](scripts/3_atac_tobias_means.py)

Purpose:
- Aggregates global TOBIAS footprint mean profiles.
- Builds quantile-normalized ATAC reference values used as normalization targets.

Outputs (via [config.py](config.py)):
- `GLOBAL_TOBIAS_MEAN`
- `GLOBAL_ATAC_MEAN`
- `GLOBAL_ATAC_QN_MIN`, `GLOBAL_ATAC_QN_MAX`, `GLOBAL_ATAC_QN_MEAN`
- `GLOBAL_ATAC_QN_MEAN_TABLE`

Command:

```bash
python scripts/3_atac_tobias_means.py
```

### Script 4: TF/tissue feature extraction (optional)

File: [scripts/4_feature_extraction.py](scripts/4_feature_extraction.py)

Purpose:
- Creates per-pair training and held-out sets with motif, ATAC, TOBIAS, TF-activity, CRUP, conservation, and co-binding features.
- Writes pair-specific `training_set.pk` and `test_set.pk`.

Required arguments:
1. tissue
2. tf
3. atac_bigwig
4. atac_bam
5. atac_peak
6. chip_peak
7. encode_tissue
8. outdir

Command template:

```bash
python scripts/4_feature_extraction.py \
	<tissue> <tf> <atac_bigwig> <atac_bam> <atac_peak_bed> <chip_peak_bed> <encode_tissue> <outdir>
```

### Script 5: general model training (recommended start)

File: [scripts/5_general_model_training.py](scripts/5_general_model_training.py)

Purpose:
- Loads all pair training matrices listed by `FEATURE_SET_GLOB`.
- Trains and saves the base general XGBoost model.

Output:
- `GENERAL_MODEL_NT` (typically [data/models/general_model_nt.json](data/models/general_model_nt.json))

Command:

```bash
python scripts/5_general_model_training.py
```

### Script 6: TF-specific model training

File: [scripts/6_tf_model_training.py](scripts/6_tf_model_training.py)

Purpose:
- Loads a pair-specific `training_set.pk`.
- Uses general-model predictions as one feature for TF-tuned model.
- Trains TF-tuned and TF-only models for a specific tissue-TF pair.

Required arguments:
1. tissue
2. tf
3. encode_tissue

Command template:

```bash
python scripts/6_tf_model_training.py <tissue> <tf> <encode_tissue>
```

### Script 7: ensemble model training

File: [scripts/7_ensemble_models.py](scripts/7_ensemble_models.py)

Purpose:
- Trains ensemble families for general, TF-only, TF-tuned, and TF+transformer variants.
- Writes fold-level summaries and model metadata under [data/models/ensemble](data/models/ensemble).

Main options:
- `--seed` (default 42)
- `--max-rounds` (default 1000)
- `--early-stopping-rounds` (default 50)
- `--device` (`cuda` or `cpu`)

Command:

```bash
python scripts/7_ensemble_models.py --device cuda
```

## Recommended run path with downloaded data pack

Most users should run:

```bash
python scripts/5_general_model_training.py
python scripts/6_tf_model_training.py <tissue> <tf> <encode_tissue>
python scripts/7_ensemble_models.py --device cuda
```

Only run scripts 1-4 if you need to regenerate features from raw FASTA/ATAC/ChIP resources.

## Prediction table artifact for figure reproduction

Figure reproduction is driven by:

- [data/all_features_predictions_leaderboard_test.parquet](data/all_features_predictions_leaderboard_test.parquet)

This table is expected to contain final-round + leaderboard rows, features, and model prediction columns including:

- general model
- tf only
- tf tuned
- tf only (ensemble)
- tf tuned (ensemble)
- tf transformer

This artifact is then used by scripts/notebooks under [figure_scripts](figure_scripts).
