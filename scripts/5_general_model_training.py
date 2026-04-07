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
from config import EMB_COLS, GEN_FEATURES, GENERAL_MODEL_NT, ensure_model_dirs, get_feature_training_files

feature_files = get_feature_training_files()
ensure_model_dirs()

all_features = EMB_COLS + GEN_FEATURES

training_sets = []
for training_file in feature_files:
    with open(training_file, 'rb') as f:
        tmp = pickle.load(f)
        training_sets.append(tmp[all_features + ['label']])
        
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
general_model.fit(training_set[all_features], training_set['label'])
general_model.save_model(GENERAL_MODEL_NT)

