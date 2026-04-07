#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 13:48:40 2024

@author: aksu
"""
import pandas as pd
import numpy as np
import qnorm
import glob
from pathlib import Path
from config import (
    ENHANCERS_BED,
    GLOBAL_ATAC_MEAN,
    GLOBAL_ATAC_QN_MIN,
    GLOBAL_ATAC_QN_MAX,
    GLOBAL_ATAC_QN_MEAN,
    GLOBAL_ATAC_QN_MEAN_TABLE,
    GLOBAL_TOBIAS_MEAN,
)

outdir = str(Path(ENHANCERS_BED).resolve().parents[1] / 'atac')
tissues = [] # List of tissues

enhancers = pd.read_csv(ENHANCERS_BED, sep='\t', names=['chr', 'start', 'end', 'enh_id'], header=None, dtype=str)
enhancers['coord'] = enhancers['chr'] + ':' + enhancers['start'] + '-' + enhancers['end']
enhancers.index = enhancers.enh_id

# TOBIAS Footprints
tobias_summaries = []
for tobias_file in glob.glob(f'{outdir}/*_footprints_summary.tsv'):
    tobias_summaries.append(pd.read_csv(tobias_file, sep='\t', header='infer'))

avg = np.nanmean([df['mean'] for df in tobias_summaries], axis=0)
tobias_mean_mean = pd.DataFrame()
tobias_mean_mean['coord'] = tobias_summaries[0]['#chrom'] + ':' + tobias_summaries[0]['start'].astype(str) + '-' + tobias_summaries[0]['end'].astype(str)
tobias_mean_mean['score'] = avg
tobias_mean_mean = tobias_mean_mean.merge(enhancers, on='coord')
tobias_mean_mean.index = tobias_mean_mean.enh_id
tobias_mean_mean = tobias_mean_mean[['coord', 'score']]
Path(GLOBAL_TOBIAS_MEAN).parent.mkdir(parents=True, exist_ok=True)
tobias_mean_mean.to_csv(GLOBAL_TOBIAS_MEAN, sep='\t', header=False)

# ATAC-seq quantile normalization
atac = {}

for tissue in tissues:
    tmp = pd.read_csv(f'{outdir}/{tissue}_atacfeatures.txt', sep='\t', header='infer')
    tmp['coord'] = tmp['#chrom'] + ':' + tmp['start'].astype(str) + '-' + tmp['end'].astype(str)
    tmp = tmp.merge(enhancers, on='coord', how='right')
    tmp.index = tmp.enh_id
    atac[tissue] = tmp[['min','max','mean']].fillna(0)

atac_qnorm = {}

for typ in ['min', 'max', 'mean']:
    df = pd.concat([b[typ] for b in atac.values()], axis=1)
    df.columns = atac.keys()
    atac_qnorm[typ] = qnorm.quantile_normalize(df)
    target = {
        'min': GLOBAL_ATAC_QN_MIN,
        'max': GLOBAL_ATAC_QN_MAX,
        'mean': GLOBAL_ATAC_QN_MEAN,
    }[typ]
    Path(target).parent.mkdir(parents=True, exist_ok=True)
    atac_qnorm[typ].iloc[:,0].to_csv(target, sep='\t', index=False, header=False)

atac_qnorm['mean'].to_csv(GLOBAL_ATAC_QN_MEAN_TABLE, sep='\t', index=True, header=True)

# Mean ATAC signal across all cells
atac_mean_mean = atac_qnorm['mean'].mean(axis=1)
atac_mean_mean.to_csv(GLOBAL_ATAC_MEAN, sep='\t', header=False)




    
