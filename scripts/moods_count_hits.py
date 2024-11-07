#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 16:14:45 2023

@author: aksu

Finds the max pwm score and the total number of hits in each sequence, using moods output file.
"""

from collections import defaultdict
import sys



def main(input_file, output_file):
    
    enh_motif_maxscores = defaultdict(lambda: float('-inf'))
    region_motif_hits = defaultdict(int)
    
    with open(input_file, 'r') as f:
        for line in f:
            # lines are of form chr1:3115483-3115493,onecut2_03.pfm,20,+,0.20247693872140649,GGATCAATAT,
            row = line.strip().split(",")
            region = row[0]
            motif = row[1].split('.')[0]
            score = float(row[4])
            region_motif_hits[(region,motif)] += 1
            
            if score > enh_motif_maxscores[(region, motif)]:
                enh_motif_maxscores[(region, motif)] = round(score, 2)
                    
    with open(output_file, 'w') as g:
        
        rows = [[enh, mot, str(score), str(region_motif_hits[(enh, mot)])] for (enh, mot), score in enh_motif_maxscores.items()]
        g.write('\n'.join(['\t'.join(x) for x in rows]))
    
    
if __name__ == '__main__':
    input_file = sys.argv[1]
    output_file = f'{input_file}.filtered'
    main(input_file, output_file)
    
#with open(output_sum, 'w') as g:
#    rows = [[enh, mot, str(score)] for (enh, mot), score in enh_motif_sums.items()]
#    g.write('\n'.join(['\t'.join(x) for x in rows]))