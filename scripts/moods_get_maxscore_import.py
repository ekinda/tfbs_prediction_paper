#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 16:14:45 2023

@author: aksu

Finds the maximum score in each sequence from a MOODS output file.
"""
from collections import defaultdict

#output_sum = f'{input_file}.sums'
#enh_motif_sums = defaultdict(float)

def moods_max(input_file, output_file):
    
    enh_motif_maxscores = defaultdict(lambda: float('-inf'))
    
    with open(input_file, 'r') as f:
        for line in f:
            # lines are of form chr1:3115483-3115493,onecut2_03.pfm,20,+,0.20247693872140649,GGATCAATAT,
            row = line.strip().split(",")
            gene = row[0]
            motif = row[1].split('.')[0]
            score = float(row[4])
            #enh_motif_sums[(gene,motif)] += score
            
            if score > enh_motif_maxscores[(gene, motif)]:
                enh_motif_maxscores[(gene, motif)] = round(score, 2)
                    
    with open(output_file, 'w') as g:
        
        rows = [[enh, mot, str(score)] for (enh, mot), score in enh_motif_maxscores.items()]
        g.write('\n'.join(['\t'.join(x) for x in rows]))

