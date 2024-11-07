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

outdir = ''
tissues = [] # List of tissues
ENHANCERS_BED = ''

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
tobias_mean_mean.to_csv(f'{outdir}/tobias_mean_mean.tsv', sep='\t')

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
    atac_qnorm[typ].iloc[:,0].to_csv(f'{outdir}/atac_{typ}_qn_values.tsv', sep='\t', index=False, header=False)

atac_qnorm['mean'].to_csv(f'{outdir}/atac_qnorm_mean.tsv', sep='\t', index=True, header=True)

# Mean ATAC signal across all cells
atac_mean_mean = atac_qnorm['mean'].mean(axis=1)
atac_mean_mean.to_csv(f'{outdir}/atac_mean_mean.tsv', sep='\t')




    
