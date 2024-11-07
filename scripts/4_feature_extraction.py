#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 13:48:40 2024

@author: aksu
"""
import pandas as pd
import numpy as np
import qnorm
import sys
import os
import subprocess
from sklearn.preprocessing import StandardScaler
from moods_get_maxscore_import import moods_max
import moods_count_hits
import pickle
import tfactcor_import
from scipy.stats import fisher_exact
import multiprocessing as mp
from functools import partial 
from time import time
from io import StringIO

#### PROGRAM BINARIES ####
TRAP_BIN=''
MOODS_BIN= ''
TOBIAS_BIN = ''
BWTOOL_BIN = ''
BEDTOOLS_BIN = ''

#### INPUTS ####
HG38='' #Human genome in fasta format
BLACKLIST = '' # Human genome blacklist file
ENHANCERS_BED = 'data/all_regions.bed' # Bed file of selected enhancers/genomic regions
ENHANCERS_FASTA = '' # Fasta file of selected enhancers/genomic regions
TRAP_BG = 'data/trap_background.GEVparams'
MOTIF_FILE = '' # Transfac TFBS motifs in fasta format
MOTIF_DIR = ''  # Directory of Transfac TFBS motifs in pfm format

inputs = {
    'tissue':sys.argv[1],
    'tf':sys.argv[2],
    'atac_bigwig':sys.argv[3],
    'atac_bams':[sys.argv[4]],
    'atac_peak':sys.argv[5],
    'chip_peak':sys.argv[6],
    'encode_tissue':sys.argv[7],
    'outdir':sys.argv[8]
}

outdir = inputs['outdir']

def timeit(func):
    def wrapper(*args, **kwargs):
        start = time()
        result = func(*args, **kwargs)
        end = time()
        print('Function', func.__name__, 'time:', round((end-start)/60,1), 'min')
        return result
    return wrapper

# 1. Run TRAP
def run_trap(inputs, outdir, gene_to_transfac):
    f = open(f'{outdir}/trap.txt', 'w')
    cmd = [TRAP_BIN,
           '-s', HG38,
           '-region', ENHANCERS_BED,
           '-matrix', MOTIF_FILE,
           '-norm', TRAP_BG,
           '-name', gene_to_transfac[inputs['tf']]
           ]

    subprocess.call(cmd, stdout=f)
    f.close()
    return
    
# 2. Run MOODS
def run_moods(inputs, outdir, gene_to_transfac):
    motif = 'V$'+gene_to_transfac[inputs["tf"]]
    f = open(f'{outdir}/moods.txt', 'w')
    cmd = [MOODS_BIN,
           '-m', "'" + MOTIF_DIR + motif + ".pfm'",
           '-s', ENHANCERS_FASTA,
           '-t', '0'
           ]
    subprocess.run(' '.join(cmd), stdout = f, shell=True)
    f.close()
    return    

# 3. Run TOBIAS, process
def run_tobias(inputs, outdir, cores=64):
    genome = HG38
    blacklist = BLACKLIST
    peaks = inputs["atac_peak"]
    
    pfm = f'{MOTIF_DIR}/V${gene_to_transfac[inputs["tf"]]}.pfm'
    
    with open(f'{outdir}/tf_pfm.txt', 'w') as f:
        f.write(f'>{inputs["tf"]}\n')
    subprocess.run(f"cat '{pfm}' >> {outdir}/tf_pfm.txt", shell=True)
    
    for i, bam in enumerate(inputs['atac_bams']):
        
        prefix = f'{inputs["tissue"]}_{str(i)}'
        
        cmd = f"{TOBIAS_BIN} BINDetect --motifs {outdir}/tf_pfm.txt --signals {outdir}/{prefix}_footprints.bw --genome {genome} --peaks {peaks} --outdir {outdir} --cores {cores} --skip-excel --naming 'name' --bound-pvalue 0.001"
        subprocess.run(cmd, shell=True)
    
    return
        
def process_tobias(inputs, outdir):

    regions = ENHANCERS_BED
    tf = inputs['tf']
        
    # Intersect bound TFBS counts with all regions
    for i in range(len(inputs['atac_bams'])):
        cmd = f'{BEDTOOLS_BIN} intersect -c -a {regions} -b {outdir}/{tf}/beds/{tf}_{tf}_{str(i)}_footprints_bound.bed | cut -f5 > {outdir}/{tf}_{str(i)}_bound.counts'
        subprocess.run(cmd, shell=True)
        
    cmd = f'paste {outdir}/{tf}*.counts > {outdir}/{tf}.allcounts'
    subprocess.run(cmd, shell=True)
    
    # Mean bound TFBS count over all samples (bio/tech. replicates)
    cmd = '''awk -v n=''' + str(len(inputs['atac_bams'])) + " '{c=0;for(i=1;i<=NF;++i){c+=$i}}{print c / n }'" + f' {outdir}/{tf}.allcounts > {outdir}/{tf}.meancounts'
    subprocess.run(cmd, shell=True)
    return
    
def get_bound_labels(inputs, outdir):
    regions = ENHANCERS_BED
    chip = inputs['chip_peak']
    cmd = f'''{BEDTOOLS_BIN} intersect -loj -F 0.5 -a {regions} -b {chip} | awk '$5 != "."' > {outdir}/bound_regions.bed'''
    subprocess.run(cmd, shell=True)
    
    return pd.read_csv(f'{outdir}/bound_regions.bed', sep='\t', header=None).iloc[:,3].values
    
def get_tfactcors(inputs, outdir, max_pval=0.001, max_qval=1):
    
    cors = []
    for line in tfactcor_import.get_correlation(inputs['tf'], max_pval, max_qval):
        cors.append([str(x) for x in line])
    
    with open(f'{outdir}/tfact_crup_cor.tsv', 'w') as g:
        g.write('tf\tenh_id\tcoef\tpval\tqval\n')
        g.write('\n'.join(['\t'.join(x) for x in cors]))
        
    return

@timeit
def atac_matrices(inputs, outdir):
    # ATAC-seq
    tmp = pd.read_csv(f'{outdir}/atacfeatures.txt', sep='\t')
    tmp['coord'] = tmp['#chrom'] + ':' + tmp['start'].astype(str) + '-' + tmp['end'].astype(str)
    tmp = tmp.merge(enhancers, on='coord', how='right')
    tmp.index = tmp.enh_id
    atac = tmp[['min','max','mean']].copy().fillna(0)
    
    # ATAC mean mean
    atac_mean_mean = pd.read_csv(f'{outdir}/atac_mean_mean.tsv', sep='\t', header=None, names=['enh_id', 'score'])
    atac_mean_mean.index = atac_mean_mean.enh_id
    atac_mean_mean = atac_mean_mean['score']

    # ATAC quantile normalization
    atac_qnorm = pd.DataFrame()
    for typ in ['min', 'mean', 'max']:
        qn_values = np.loadtxt(f'{outdir}/atac_{typ}_qn_values.tsv')
        atac_qnorm[typ] = qnorm.quantile_normalize(atac[[typ]], target=qn_values)
        atac_qnorm[f'delta_{typ}'] = atac_qnorm[typ] - atac_mean_mean
        
    return atac_qnorm, atac_mean_mean

def cotracte_regions(inputs, outdir, enhancers, atac_qnorm, n=5000):
    atac_qnorm_all = pd.read_csv(f'{outdir}/atac_qnorm_mean.tsv', sep='\t', index_col='enh_id')
    atac_qnorm_all = pd.concat([atac_qnorm_all, atac_qnorm['mean']], axis=1)
    atac_qnorm_delta = atac_qnorm_all.sub(atac_qnorm_all.mean(axis=1), axis=0)
    ubiq = atac_qnorm_all.loc[(atac_qnorm_all > 5).sum(axis=1) == 12].head(n)
    ubiq.reindex(ubiq.mean(axis=1).sort_values(ascending=False).index)
    enhancers.loc[ubiq.index, ['chr', 'start', 'end']].to_csv(f'{outdir}/cotracte_ubiq.bed', sep='\t', header=False, index=False)
    
    atac_qnorm_all['delta'] = atac_qnorm_delta['mean']
    atac_qnorm_all['delta_norm'] = atac_qnorm_all['delta'] / atac_qnorm_all['mean']
    specific = atac_qnorm_all[atac_qnorm_all['mean'] > 4].sort_values(by='delta_norm', ascending=False).head(n)
    enhancers.loc[specific.index, ['chr', 'start', 'end']].to_csv(f'{outdir}/cotracte_specific.bed', sep='\t', header=False, index=False)
    return

@timeit
def cotracte_calculation(inputs, outdir, enhancers, n=5000, k=1000):
    trap_bin = TRAP_BIN
    hg = HG38
    pfms = MOTIF_FILE
    trap_bg = TRAP_BG
    
    cmd = f'{trap_bin} -s {hg} -region {outdir}/cotracte_specific.bed -matrix {pfms} -norm {trap_bg} -thread 64 > {outdir}/cotracte_trap_specific.txt'
    subprocess.run(cmd, shell=True)
    
    cmd = f'{trap_bin} -s {hg} -region {outdir}/cotracte_ubiq.bed -matrix {pfms} -norm {trap_bg} -thread 64 > {outdir}/cotracte_trap_ubiq.txt'
    subprocess.run(cmd, shell=True)

    trap_ubiq = pd.read_csv(f'{outdir}/cotracte_trap_ubiq.txt', sep='\t', skiprows=1, header=None, names=['chr', 'start', 'end', 'motif', 'pval'])
    trap_ubiq['coord'] = trap_ubiq.chr + ':' + (trap_ubiq.start -1).astype(str) + '-' + trap_ubiq.end.astype(str)
    trap_ubiq = trap_ubiq.merge(enhancers, on='coord', how='left')
    trap_ubiq = trap_ubiq[['enh_id', 'motif', 'pval']]
    trap_ubiq['ubiq'] = True
    
    trap_sp = pd.read_csv(f'{outdir}/cotracte_trap_specific.txt', sep='\t', skiprows=1, header=None, names=['chr', 'start', 'end', 'motif', 'pval'])
    trap_sp['coord'] = trap_sp.chr + ':' + (trap_sp.start -1).astype(str) + '-' + trap_sp.end.astype(str)
    trap_sp = trap_sp.merge(enhancers, on='coord', how='left')
    trap_sp = trap_sp[['enh_id', 'motif', 'pval']]
    trap_sp['ubiq'] = False
    
    trap_all = pd.concat([trap_sp, trap_ubiq])
    trap_all.to_csv(f'{outdir}/trap_all.csv')
    return list(trap_all.motif.unique())

def corr_mp(tf2, tf1, tf1_bound_specific, tf1_bound_ubiq,
            n, k, bound_together_threshold, L_threshold):
    
    cmd = f"grep {tf2[2:]} {outdir}/trap_all.csv"
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    out = StringIO(p.communicate()[0].decode('utf-8'))
    
    try:
        tf2_df = pd.read_csv(out, header=None, names=['enh_id', 'motif', 'pval', 'ubiq']).sort_values(by='pval').head(k)
    except Exception as error:
        print(error)
        pass
    
    tf2_bound_ubiq = set(tf2_df[tf2_df.ubiq == True].enh_id.values)
    tf2_bound_specific = set(tf2_df[tf2_df.ubiq == False].enh_id.values)

    sp_both = len(tf1_bound_specific.intersection(tf2_bound_specific))
    if sp_both < bound_together_threshold:
        return []
        
    sp_tf1 = len(tf1_bound_specific) - sp_both
    sp_tf2 = len(tf2_bound_specific) - sp_both
    sp_none = n - sp_both - sp_tf1 - sp_tf2
    ubiq_both = len(tf1_bound_ubiq.intersection(tf2_bound_ubiq))
    ubiq_tf1 = len(tf1_bound_ubiq) - ubiq_both
    ubiq_tf2 = len(tf2_bound_ubiq) - ubiq_both
    ubiq_none = n - ubiq_both - ubiq_tf1 - ubiq_tf2
    
    ubiq_table = np.array([[ubiq_both, ubiq_tf1], [ubiq_tf2, ubiq_none]]) + 1
    specific_table = np.array([[sp_both, sp_tf1], [sp_tf2, sp_none]]) + 1
    
    try:
        ubiq_p = fisher_exact(ubiq_table, alternative='greater')[1]
        specific_p = fisher_exact(specific_table, alternative='greater')[1]
    except ValueError:
        print("VALUE ERROR")
        return []
    
    L = -np.log(specific_p) + np.log(ubiq_p)
    if L > L_threshold:
        return [tf1, tf2, str(round(L,2))]
    else:
        return []
    
@timeit
def cotracte_get_tfs(inputs, outdir, gene_to_transfac, tf_list, top_k_tfs = 6, n=5000, k=1000, bound_together_threshold = 10, L_threshold = 2, cores=64):
    tf1 = 'V$'+gene_to_transfac[inputs['tf']]
    
    cmd = f'grep {tf1[2:]} {outdir}/trap_all.csv'
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    out = StringIO(p.communicate()[0].decode('utf-8'))
    
    try:
        tf1_df = pd.read_csv(out, header=None, names=['enh_id', 'motif', 'pval', 'ubiq']).sort_values(by='pval').head(k)
    except Exception as error:
        print(error)
        pass
    
    tf1_bound_ubiq = set(tf1_df[tf1_df.ubiq == True].enh_id.values)
    tf1_bound_specific = set(tf1_df[tf1_df.ubiq == False].enh_id.values)
    
    pool = mp.Pool(cores)
    print("starting pool. n =", len(tf_list))

    results = pool.map(partial(corr_mp, tf1=tf1, tf1_bound_specific=tf1_bound_specific, tf1_bound_ubiq=tf1_bound_ubiq,
                                  n=n, k=k, bound_together_threshold=bound_together_threshold, L_threshold=L_threshold), tf_list)
    pool.close()
    pool.join()
    print("end pool")
    
    with open(f'{outdir}/cotracte_results.txt', 'w') as f:
        f.write('\n'.join(['\t'.join(x) for x in results if len(x) != 0]))
    
    result_df = pd.read_csv(f'{outdir}/cotracte_results.txt', sep='\t', names=['tf1', 'tf2', 'L'])
    
    chosen = result_df.sort_values(by='L', ascending=False).head(6).tf2.values
    
    with open(f'{outdir}/cotracte_chosen_tfs.txt', 'w') as f:
        f.write('\n'.join(chosen))

    # Copy pfms to new directory
    os.makedirs(f'{outdir}/cotracte_pfms', exist_ok=True)
    for motif in chosen:
        cmd = f"cp '{MOTIF_DIR}/{motif}.pfm' {outdir}/cotracte_pfms"
        subprocess.run(cmd, shell=True)
    
    # Run MOODS and process
    
    cmd = f'{MOODS_BIN} -m {outdir}/cotracte_pfms/*.pfm -s {ENHANCERS_FASTA} --batch -p 0.001 > {outdir}/cotracte_moodsout.txt'
    subprocess.run(cmd, shell=True)
    
    moods_count_hits.main(f'{outdir}/cotracte_moodsout.txt', f'{outdir}/cotracte_moods_hits.txt')
    
    return

def construct_datamatrix(inputs, outdir, enh_ids, atac_qnorm, atac_mean_mean, enhancers, cons, crup, crup_mean, remap, bulk_mrna, tf_act):
    
    # ATAC done by seperate function
    
    # TRAP
    trap = pd.read_csv(f'{outdir}/trap.txt', sep='\t', skiprows=1, header=None, names=['chr', 'start', 'end', 'motif', 'pval'])
    trap['coord'] = trap['chr'] + ':' + (trap['start'] - 1).astype(str) + '-' + trap['end'].astype(str)
    trap = trap.merge(enhancers.reset_index(drop=True), on='coord', how='right')
    trap.index = trap['enh_id']
    trap = trap[['motif', 'pval']].fillna(0)
    
    # This sometimes happens for some reason... same row appears twice e.g. in TP53.
    # Maybe duplicated motif in transfac file?
    if len(trap) > len(enhancers):
        trap = trap.drop_duplicates()
        
        if len(trap) != len(enhancers):
            print(trap.head())
            print("Error in trap output.")
            exit(1)
    else:
        print(trap.head())
        print("Error in trap output.")
        exit(1)
    
    # MOODS
    maxpwm = pd.read_csv(f'{outdir}/moods_maxpwm.txt', sep='\t', header=None, names=['coord', 'motif', 'score'])
    maxpwm = maxpwm.merge(enhancers, on='coord', how='right')
    maxpwm.index = maxpwm.enh_id
    maxpwm = maxpwm.loc[enhancers.index, ['motif', 'score']].fillna(0)
    
    # TOBIAS Footprints
    tobias_summaries = []
    for i, bam in enumerate(inputs['atac_bams']):
        prefix = f'{inputs["tf"]}_{str(i)}'
        tobias_summaries.append(pd.read_csv(f'{outdir}/{prefix}_footprints_summary.tsv', sep='\t', header='infer'))
        
    avg = np.nanmean([df['mean'] for df in tobias_summaries], axis=0)
    tobias_mean = pd.DataFrame()
    tobias_mean['coord'] = tobias_summaries[0]['#chrom'] + ':' + tobias_summaries[0]['start'].astype(str) + '-' + tobias_summaries[0]['end'].astype(str)
    tobias_mean['tobias_avg'] = avg
    tobias_mean = tobias_mean.merge(enhancers, on='coord')
    tobias_mean.index = tobias_mean.enh_id
    tobias_mean = tobias_mean[['coord', 'tobias_avg']]
    tobias_mean.columns = ['coord', 'score']
    del avg, tobias_summaries

    # TOBIAS Deltas
    tobias_mean_mean = pd.read_csv(f'{outdir}/tobias_mean_mean.tsv', sep='\t', header=None, names=['enh_id', 'score'])
    tobias_mean_mean.index = tobias_mean_mean.enh_id
    tobias_mean_mean = tobias_mean_mean['score']
    
    tobias_deltas = tobias_mean['score'] - tobias_mean_mean

    # TOBIAS Counts
    tobias_counts = pd.read_csv(f'{outdir}/{inputs["tf"]}.meancounts', sep='\t', header=None, names=['tobias_count'])
    tobias_counts = pd.concat([enhancers,tobias_counts], axis=1)
    
    # Construct final set, start with transformer embeddings
    data = pd.read_csv(ENHANCERS_BED, sep='\t', names=['chrom', 'start', 'end', 'enh_id'])
    data.index = data.enh_id
    data.drop(['enh_id',], axis=1, inplace=True)
    
    tissue = inputs['tissue']
    tf = inputs['tf']

    data['tissue'] = tissue
    data['tf'] = tf    
    data['maxpwm'] = maxpwm.loc[enh_ids, 'score']
    data['cons'] = cons.loc[enh_ids]
    data['crup'] = crup.loc[enh_ids, inputs['encode_tissue']]
    data['crup_mean'] = crup_mean.loc[enh_ids]
    data['crup_delta'] = data['crup'] - data['crup_mean']
    data['remap'] = remap.loc[enh_ids]
    data['tf_exp'] = bulk_mrna.loc[tf, inputs['encode_tissue']]
    data['tf_activity'] = tf_act.loc[inputs['encode_tissue'], tf]
    
    tf_act_cors = pd.read_csv(f'{outdir}/tfact_crup_cor.tsv', sep='\t')
    tf_act_cors.index = tf_act_cors.enh_id
    
    data['tfact_crupcor_coef'] = [tf_act_cors.loc[enh, 'coef'] if enh in tf_act_cors.index else 0 for enh in enh_ids]
    data['tfact_crupcor_pval'] = -np.log(0.000001+np.array([tf_act_cors.loc[enh, 'pval'] if enh in tf_act_cors.index else 1 for enh in enh_ids]))

    #data['tf_crup_cor'] = np.array(-np.log(0.001 + tf_crup_cors[tf].loc[enh_ids].fillna(1).pval))
    #data.loc[data.tf_crup_cor < 0.001, 'tf_crup_cor'] = 0
    
    for typ in ['min', 'max', 'mean']:
        data[f'atac_{typ}'] = atac_qnorm.loc[enh_ids, typ]
        data[f'atac_delta_{typ}'] = atac_qnorm.loc[enh_ids, f'delta_{typ}']
        
    data['atac_mean_mean'] = atac_mean_mean.loc[enh_ids]
    data['trap'] = -np.log(0.0001 + trap.loc[enh_ids, 'pval'])
    data['tobias_avg'] = tobias_mean.loc[enh_ids, 'score']
    data['delta_tobias_avg'] = tobias_deltas.loc[enh_ids]
    data['tobias_mean_mean'] = tobias_mean_mean.loc[enh_ids]
    data['tobias_count'] = tobias_counts.loc[enh_ids, 'tobias_count']
    
    # Labels
    bound_ids = get_bound_labels(inputs, outdir)
    data['label'] = data.index.isin(bound_ids)
    
    # Cotracte
    cotracte = pd.read_csv(f'{outdir}/cotracte_moods_hits.txt', sep='\t', names=['coord', 'motif', 'maxpwm', 'hits'])
    cotracte = cotracte.pivot(index='coord', columns='motif', values=['maxpwm', 'hits']).fillna(0)
    
    try:
        cotracte.columns = [f'cot_maxpwm_{i}' for i in range(6)] + [f'cot_hits_{i}' for i in range(6)]
        cotracte = cotracte.merge(enhancers, on='coord', how='right').fillna(0)
        cotracte.index = cotracte.enh_id
    except Exception as error:
        print(error)
        print("Problem in cotracte output. Assuming no co-binding TFs found.")
        cotracte = np.empty((len(enhancers), 12))
        cotracte[:] = np.nan
        cotracte = pd.DataFrame(cotracte)
        cotracte.columns = [f'cot_maxpwm_{i}' for i in range(6)] + [f'cot_hits_{i}' for i in range(6)]
    
    data[[f'cot_maxpwm_{i}' for i in range(6)] + [f'cot_hits_{i}' for i in range(6)]] = cotracte.loc[enh_ids, [f'cot_maxpwm_{i}' for i in range(6)] +
                                                                                                              [f'cot_hits_{i}' for i in range(6)]]
    return data


# Starting feature extraction

print(inputs)
os.makedirs(outdir, exist_ok=True)

# 1. General data
gene_to_transfac = {}
gene_to_transfacfile = 'data/gene_symbol_to_transfac.txt'
with open(gene_to_transfacfile, 'r') as f:
    for line in f:
        gene, transfac = line.strip().split()
        gene_to_transfac[gene] = transfac[2:]
        
transfac_to_gene = {f'V${val}':key for (key,val) in gene_to_transfac.items()}

with open('data/chosen_encode_tissues.txt', 'r') as f:
    encode_tissues = f.read().strip().split()

bulk_mrna = pd.read_csv('data/tfs_meantpm.tsv', sep='\t', index_col='gene_symbol')
bulk_mrna = bulk_mrna[~bulk_mrna.index.duplicated(keep='first')]
bulk_mrna.drop('gene_id', axis=1, inplace=True)
bulk_mrna.fillna(0, inplace=True)
bulk_mrna = np.log1p(bulk_mrna)
bulk_mrna = bulk_mrna.loc[bulk_mrna.sum(axis=1) > 1,:]

tf_act = pd.read_csv('data/tf_activities_fixed.tsv', sep='\t', index_col=0)
tf_act = tf_act*10000
tf_act = pd.DataFrame(StandardScaler().fit_transform(tf_act), index=tf_act.index, columns=tf_act.columns) # Standardization

enhancers = pd.read_csv(ENHANCERS_BED, sep='\t', names=['chr', 'start', 'end', 'enh_id'], header=None, dtype=str)
enhancers['coord'] = enhancers['chr'] + ':' + enhancers['start'] + '-' + enhancers['end']
enhancers.index = enhancers.enh_id

cons = pd.read_csv('data/all_regions_phastcons.tsv', sep='\t')
cons['coord'] = cons['#chrom'] + ':' + cons['start'].astype(str) + '-' + cons['end'].astype(str)
cons = cons.merge(enhancers, on='coord', how='right')
cons.index = cons.enh_id
cons = cons.loc[enhancers.index, 'mean'].fillna(0)

remap = pd.read_csv('data/all_regions_remap.tsv', sep='\t')
remap['coord'] = remap['#chrom'] + ':' + remap['start'].astype(str) + '-' + remap['end'].astype(str)
remap = remap.merge(enhancers, on='coord', how='right')
remap.index = remap.enh_id
remap = remap.loc[enhancers.index, 'mean'].fillna(0)

crup = pd.read_csv('data/all_regions_crupscores.bed', sep='\t', names=['chr', 'start', 'end', 'enh_id'] + encode_tissues,
                    header=None, index_col='enh_id').drop(columns=['chr', 'start', 'end']).fillna(0)
crup_mean = crup.mean(axis=1)


run_trap(inputs, outdir, gene_to_transfac)
run_moods(inputs, outdir, gene_to_transfac)
moods_max(f'{outdir}/moods.txt', f'{outdir}/moods_maxpwm.txt')
run_tobias(inputs, outdir)
process_tobias(inputs, outdir)
get_tfactcors(inputs, outdir)
atac_qnorm, atac_mean_mean = atac_matrices(inputs, outdir)
cotracte_regions(inputs, outdir, enhancers, atac_qnorm)
tf_list = cotracte_calculation(inputs, outdir, enhancers)
cotracte_get_tfs(inputs, outdir, gene_to_transfac, tf_list)

train_ids = enhancers.loc[~enhancers.chr.isin(['chr1', 'chr8', 'chr21']), 'enh_id'].values
train_data = construct_datamatrix(inputs, outdir, train_ids, atac_qnorm, atac_mean_mean, enhancers, cons, crup, crup_mean,remap, bulk_mrna, tf_act)
with open(f'{outdir}/training_set.pk', 'wb') as f:
    pickle.dump(train_data, f)

test_ids = enhancers.loc[enhancers.chr.isin(['chr1', 'chr8', 'chr21']), 'enh_id'].values
test_data = construct_datamatrix(inputs, outdir, test_ids, atac_qnorm, atac_mean_mean, enhancers, cons, crup, crup_mean,remap, bulk_mrna, tf_act)
with open(f'{outdir}/test_set.pk', 'wb') as f:
    pickle.dump(test_data, f)
    
