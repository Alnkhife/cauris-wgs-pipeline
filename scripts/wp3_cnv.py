#!/usr/bin/env python3
"""
WP3: coverage-based copy-number / amplification screen at virulence and
resistance loci (reviewer: "pathogenicity features were not investigated").

Normalised per-gene coverage (log2 ratio vs sample median) across 19 isolates.
Loci of interest: ALS-family adhesins, SAP proteases, efflux pumps (CDR1/2, MDR1),
HSP genes, plus key resistance genes as internal controls.

Input : analysis/01_alignment/stats/SampleXX.mosdepth.per-gene.bed.gz  (mosdepth --by genes.bed)
        data/raw/reference/genes.bed
Output: analysis/03_virulence_cnv/cnv_table.tsv
        analysis/03_virulence_cnv/cnv_heatmap.pdf
"""
import gzip

import statistics

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = os.path.expanduser('~/Projects/cauris-wgs-analysis')
STATS = os.path.join(ROOT, 'analysis/01_alignment/stats')
GENES_BED = os.path.join(ROOT, 'data/raw/reference/genes.bed')
OUT = os.path.join(ROOT, 'analysis/03_virulence_cnv')
os.makedirs(OUT, exist_ok=True)

SAMPLES = ['Sample01'] + [f'Sample{i:02d}' for i in range(3, 21)]

# loci of interest: gene name patterns to report (substring match on gene attr)
INTEREST = ['CDR1', 'CDR2', 'MDR1', 'ERG11', 'FCA1', 'GSC1', 'FUR1', 'ERG3', 'ERG6',
            'SAP9', 'ABC1', 'SNQ2', 'YCF1', 'ATM1', 'MDL1', 'PXA', 'HSP90', 'HSP70',
            'HSP104', 'HSP60', 'HSP78', 'HSP21']
# B8441 locus tags for named loci (GFF gene= attribute absent for these)
TAGS = ['9J08_02922', '9J08_04280', '9J08_03698', '9J08_02123', '9J08_01838',
        '9J08_01933', '9J08_01595', '9J08_04223', '9J08_03031']


def load_genes_bed(path):
    """genes.bed: chrom, start, end, name (from GFF)."""
    genes = []
    if not os.path.exists(path):
        return genes
    with open(path) as f:
        for line in f:
            p = line.rstrip().split('\t')
            if len(p) >= 4:
                genes.append((p[0], int(p[1]), int(p[2]), p[3]))
    return genes


def load_per_gene(sample):
    path = os.path.join(STATS, f'{sample}.per-gene.tsv')
    df = pd.read_csv(path, sep='\t')
    return df.set_index('gene')['mean_cov']


def main():
    if not os.path.exists(GENES_BED):
        print('genes.bed missing - run mosdepth with --by first')
        return
    data = {}
    for sm in SAMPLES:
        data[sm] = load_per_gene(sm)
    all_genes = sorted(set().union(*[d.index for d in data.values()]))
    table = pd.DataFrame({sm: data[sm].reindex(all_genes) for sm in SAMPLES})
    table = table.fillna(0)
    # normalise: log2(gene_cov / sample_median) where sample median is over genes with cov>0
    norm = pd.DataFrame(index=all_genes)
    for sm in SAMPLES:
        c = table[sm]
        med = c[c > 0].median()
        norm[sm] = np.log2((c + 0.5) / (med + 0.5))
    norm.to_csv(os.path.join(OUT, 'cnv_table.tsv'), sep='\t')
    print(f'cnv_table.tsv: {norm.shape[0]} genes x {norm.shape[1]} samples')

    # report the loci of interest
    pat = '|'.join(INTEREST)
    sel = norm[norm.index.str.upper().str.contains(pat, regex=True, na=False)]
    # add locus-tag matched genes, labelled with their gene names
    tag_names = {'9J08_02922': 'GSC1(FKS1)', '9J08_04280': 'FCA1(FCY1)',
                 '9J08_03698': 'ERG11', '9J08_02123': 'CDR1', '9J08_01838': 'MDR1',
                 '9J08_01933': 'FUR1', '9J08_01595': 'ERG3', '9J08_04223': 'ERG6',
                 '9J08_03031': 'SAP9'}
    for tag, gname in tag_names.items():
        if tag in norm.index:
            sel = pd.concat([sel, norm.loc[[tag]].rename(index={tag: gname})])
    sel = sel[~sel.index.duplicated(keep='first')]
    print(f'interest loci: {sel.shape[0]}')
    print(sel.round(2).to_string())

    # heatmap
    plt.figure(figsize=(14, max(6, 0.3 * sel.shape[0])))
    sns.heatmap(sel, cmap='RdBu_r', vmin=-1.5, vmax=1.5, center=0,
                xticklabels=True, yticklabels=True,
                cbar_kws={'label': 'log2(gene coverage / sample median)'})
    plt.title('Per-gene coverage ratio at virulence/resistance loci (WP3)')
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, 'cnv_heatmap.pdf'))
    plt.savefig(os.path.join(OUT, 'cnv_heatmap.png'), dpi=150)
    print('heatmap written')


if __name__ == '__main__':
    main()
