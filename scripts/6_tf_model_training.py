#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 13:48:40 2024

@author: aksu
"""
import sys
import pickle
import pandas as pd
import numpy as np
import xgboost as xgb

#### INPUTS ####
tissue = sys.argv[1]
tf = sys.argv[2]
encode_tissue = sys.argv[3]

ENHANCERS_BED = 'data/all_regions.bed'
outdir = f'/{tissue}-{tf}'
model_dir = 'models'
########

enhancers = pd.read_csv(ENHANCERS_BED, sep='\t', names=['chr', 'start', 'end', 'enh_id'], header=None, dtype=str)
enhancers['coord'] = enhancers['chr'] + ':' + enhancers['start'] + '-' + enhancers['end']
enhancers.index = enhancers.enh_id

general_model = xgb.XGBClassifier()
general_model.load_model('{model_dir}/general_model.json')

# Training sets
with open(f'{outdir}/training_set.pk', 'rb') as f:
    training_set = pickle.load(f)

# Feature sets
emb_cols = [f'emb_{i}' for i in range(1024)]

gen_features = ['maxpwm', 'cons', 'crup', 'crup_mean', 'crup_delta', 'remap', 'tf_exp', 'tf_activity',
            'tfact_crupcor_coef', 'tfact_crupcor_pval', 'trap', 'atac_min', 'atac_max', 'atac_mean', 'atac_mean_mean',
            'atac_delta_min', 'atac_delta_max', 'atac_delta_mean', 'tobias_avg', 'delta_tobias_avg', 'tobias_mean_mean', 'tobias_count',
           'cot_hits_0', 'cot_hits_1', 'cot_hits_2', 'cot_hits_3',
           'cot_maxpwm_0', 'cot_maxpwm_1', 'cot_maxpwm_2', 'cot_maxpwm_3']

tftuned_features = gen_features + ['xgb_general']
t_features = emb_cols + gen_features
tfonly_features = gen_features

# Training new models
params = {
    'booster':'gbtree',
    'objective':'reg:logistic',
    'random_state':1,
    'subsample':0.8,
    'colsample_bytree':0.8,
    'learning_rate': 0.05,
    'device':'cuda',
    'max_depth': 6,
    'min_child_weight':1
}

xgb_gen = general_model.predict(training_set[gen_features])
training_set['xgb_general'] = xgb_gen

def train_model(training, features, name, num_round=500):
    training_set = xgb.DMatrix(np.array(training[features]), label=np.array(training['label']), feature_names = features)
    bst = xgb.train(params, training_set, num_round)
    bst.save_model(f'{outdir}/{name}.json')
    bst.dump_model(f'{outdir}/{name}_readable.txt')

# TF-tuned
train_model(training_set, tftuned_features, 'tf_tuned')
# TF-only
train_model(training_set, tfonly_features, 'tf_only')
# Transformer v2
train_model(training_set, t_features, 'transformer')


