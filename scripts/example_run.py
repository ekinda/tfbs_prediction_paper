import subprocess

'''
First, you need to run scripts numbered 1, 2 and 3 to get NT embeddings and mean ATAC features across tissues.
Afterwards, this script can be modified to run feature extraction (in parallel if needed, adjust the command below)
'''

all_inputs = [
{
'tissue':'liver',
'tf':'RXRA',
'atac_bigwig':'encode-atac/bigwig/liver_ENCFF646EEM.bigWig',
'atac_bams':['encode-atac/bams/liver_ENCFF854OTB.bam'],
'atac_peak':'encode-atac/peaks/liver_ENCFF347HUL.bed',
'chip_peak':'encode_chipseq/liver/RXRA_ENCFF537JEP.bed',
'encode_tissue':'liver', # You need to choose 1 tissue from the list of 80
'outdir':'data/liver-RXRA'
},
{
'tissue':'liver',
'tf':'RAD21',
'atac_bigwig':'encode-atac/bigwig/liver_ENCFF646EEM.bigWig',
'atac_bams':['encode-atac/bams/liver_ENCFF854OTB.bam'],
'atac_peak':'encode-atac/peaks/liver_ENCFF347HUL.bed',
'chip_peak':'encode_chipseq/liver/RAD21_ENCFF893YQZ.bed',
'encode_tissue':'liver', # You need to choose 1 tissue from the list of 80
'outdir':'data/liver-RAD21'
}
]

for i, inputs in enumerate(all_inputs):
    a = list(inputs.values())
    cmd = f'python scripts/4_feature_extraction.py {a[0]} {a[1]} {a[2]} {a[3][0]} {a[4]} {a[5]} {a[6]} {a[7]}'
    print(cmd)
    subprocess.run(cmd, shell=True)