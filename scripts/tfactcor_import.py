#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 16:51:19 2024

@author: ekinda

"""
import pandas as pd
import numpy as np
import qvalue_nfusi as qvalue
from scipy.stats import pearsonr

def get_correlation(tf, max_pval=0.001, max_qval=1):

    
    with open('/project/primate_msa/egrn/encode_samples/chosen_encode_tissues.txt', 'r') as f:
        tissues = f.read().strip().split()
        
    tf_act = pd.read_csv('/project/primate_msa/egrn/tf_enh_links/tf_activities_fixed.tsv', sep='\t', index_col=0)
    crup = pd.read_csv('/project/primate_msa/egrn/enhancers/all_regions_crupscores.bed', sep='\t',
                       names=['chr', 'start', 'end', 'enh_id'] + tissues, header=None, index_col='enh_id').drop(columns=['chr', 'start', 'end']).dropna()

    #max_pval = 0.001
    #max_qval = 1
    
    print(f'Tf {tf} running\n', flush=True)
    
    cors = crup[tissues].apply((lambda x: pearsonr(x, tf_act.loc[:,tf])), axis=1)
    
    pvals = [x[1] for x in cors]
    coefs = [x[0] for x in cors]
    
    # Code from github/nfusi
    result = []
    try:
        qvals = qvalue.estimate(np.array(pvals))
    except AssertionError:
        for i,enh_name in enumerate(crup.index):
            p = pvals[i]
            coef = coefs[i]
            if p <= max_pval:
                result.append([tf, enh_name, coef, p, -1])
        return result
        
    for i,enh_name in enumerate(crup.index):
        
        p = pvals[i]
        q = qvals[i]
        coef = coefs[i]
        
        if p <= max_pval and q <= max_qval:
            result.append([tf, enh_name, str(round(coef,3)), str(round(p,5)), str(round(q,5))])

    return result
