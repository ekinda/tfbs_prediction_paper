#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 13:48:40 2024

@author: aksu
"""
import pickle
import pandas as pd
import numpy as np
import xgboost as xgb

ENHANCERS_BED = 'data/all_regions.bed'
outdirs = [] # List of directories of feature extraction
model_dir = 'models'

enhancers = pd.read_csv(ENHANCERS_BED, sep='\t', names=['chr', 'start', 'end', 'enh_id'], header=None, dtype=str)
enhancers['coord'] = enhancers['chr'] + ':' + enhancers['start'] + '-' + enhancers['end']
enhancers.index = enhancers.enh_id

# Feature set
gen_features = ['maxpwm', 'cons', 'crup', 'crup_mean', 'crup_delta', 'remap', 'tf_exp', 'tf_activity',
            'tfact_crupcor_coef', 'tfact_crupcor_pval', 'trap', 'atac_min', 'atac_max', 'atac_mean', 'atac_mean_mean',
            'atac_delta_min', 'atac_delta_max', 'atac_delta_mean', 'tobias_avg', 'delta_tobias_avg', 'tobias_mean_mean', 'tobias_count',
           'cot_hits_0', 'cot_hits_1', 'cot_hits_2', 'cot_hits_3',
           'cot_maxpwm_0', 'cot_maxpwm_1', 'cot_maxpwm_2', 'cot_maxpwm_3']

training_sets = []
for outdir in outdirs:
    with open(f'{outdir}/training_set.pk', 'rb') as f:
        tmp = pickle.load(f)
        training_sets.append(tmp[gen_features])
        
training_set = pd.concat(training_sets)
del training_sets

params = {
    'booster':'gbtree',
    'objective':'reg:logistic',
    'random_state':1,
    'subsample':0.8,
    'colsample_bytree':0.8,
    'learning_rate': 0.05,
    'device':'cuda',
    'max_depth': 6,
    'min_child_weight':1,
    'n_estimators':500,
}

general_model = xgb.XGBRegressor(**params)
general_model.fit(training_set[gen_features], training_set['label'])
general_model.save_model('{model_dir}/general_model.json')

