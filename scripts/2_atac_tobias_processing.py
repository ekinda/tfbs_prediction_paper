#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 13:48:40 2024

@author: aksu
"""
import pandas as pd
import sys
import os
import subprocess
from time import time

#### PROGRAM BINARIES ####
TOBIAS_BIN = ''
BWTOOL_BIN = ''
BEDTOOLS_BIN = ''

#### INPUTS ####
HG38='' #Human genome in fasta format
BLACKLIST = '' # Human genome blacklist file
ENHANCERS_BED = 'data/all_regions.bed' # Bed file of selected enhancers/genomic regions
ENHANCERS_FASTA = '' # Fasta file of selected enhancers/genomic regions

inputs = {
    'tissue':sys.argv[1],
    'atac_bigwig':sys.argv[3],
    'atac_bams':[sys.argv[4]],
    'atac_peak':sys.argv[5],
    }

outdir = 'data/atac'
########

print(inputs)
os.makedirs(outdir, exist_ok=True)

def timeit(func):
    def wrapper(*args, **kwargs):
        start = time()
        result = func(*args, **kwargs)
        end = time()
        print('Function', func.__name__, 'time:', round((end-start)/60,1), 'min')
        return result
    return wrapper

# Run TOBIAS, process
def run_tobias(inputs, outdir, cores=64):
    genome = HG38
    blacklist = BLACKLIST
    peaks = inputs["atac_peak"]
    
    for i, bam in enumerate(inputs['atac_bams']):
        
        prefix = f'{inputs["tissue"]}_{str(i)}'
        
        cmd = f'{TOBIAS_BIN} ATACorrect --bam {bam} --genome {genome} --peaks {peaks} --blacklist {blacklist} --outdir {outdir} --cores {cores} --prefix {prefix}'
        subprocess.run(cmd, shell=True)
        
        cmd = f'{TOBIAS_BIN} ScoreBigwig --signal {outdir}/{prefix}_corrected.bw --regions {peaks} --output {outdir}/{prefix}_footprints.bw --cores 40'
        subprocess.run(cmd, shell=True)
    
    return
        
def process_tobias(inputs, outdir):

    regions = ENHANCERS_BED
    
    # Summarize the bigwig
    for i, bam in enumerate(inputs['atac_bams']):
        prefix = f'{inputs["tissue"]}_{str(i)}'
        cmd = f'{BWTOOL_BIN} summary {regions} {outdir}/{prefix}_footprints.bw {outdir}/{prefix}_footprints_summary.tsv -header'
        subprocess.run(cmd, shell=True)
        
    return

def get_atac_features(inputs, outdir):
    regions = ENHANCERS_BED
    cmd = f'{BWTOOL_BIN} summary {regions} {inputs["atac_bigwig"]} {outdir}/{inputs["tissue"]}_atacfeatures.txt -header -skip-median'
    subprocess.run(cmd, shell=True)
    
enhancers = pd.read_csv(ENHANCERS_BED, sep='\t', names=['chr', 'start', 'end', 'enh_id'], header=None, dtype=str)
enhancers['coord'] = enhancers['chr'] + ':' + enhancers['start'] + '-' + enhancers['end']
enhancers.index = enhancers.enh_id

run_tobias(inputs, outdir)
process_tobias(inputs, outdir)
get_atac_features(inputs, outdir)

    
