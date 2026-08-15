#!/usr/bin/env python3
"""WP6: merge context-sample genotypes into the study SNP alignment, run
snp-dists + IQ-TREE2, and report where the 19 study isolates fall relative
to the Saudi context isolates."""
import os
import subprocess
import sys

ROOT = os.path.expanduser('~/Projects/cauris-wgs-analysis')
ENV = os.environ.get('CONDA_PREFIX', '') + '/bin'
OUT = os.path.join(ROOT, 'analysis/07_context')
GENO = os.path.join(OUT, 'genotypes')
SITES = os.path.join(OUT, 'sites.tsv')
STUDY_ALN = os.path.join(ROOT, 'analysis/05_phylogeny/freebayes_snps.fasta')
COMBINED = os.path.join(OUT, 'combined_snps.fasta')

STUDY = ['Sample01', 'Sample03', 'Sample04', 'Sample05', 'Sample06', 'Sample07',
         'Sample08', 'Sample09', 'Sample10', 'Sample11', 'Sample12', 'Sample13',
         'Sample14', 'Sample15', 'Sample16', 'Sample17', 'Sample18', 'Sample19', 'Sample20']


def load_fasta(path):
    seqs, cur = {}, None
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                cur = line[1:].split()[0]
                seqs[cur] = []
            else:
                seqs[cur].append(line)
    return {k: ''.join(v) for k, v in seqs.items()}


def load_sites(path):
    sites = []
    with open(path) as f:
        for line in f:
            p = line.rstrip().split('\t')
            if len(p) == 2:
                sites.append((p[0], int(p[1])))
    return sites


def main():
    sites = load_sites(SITES)
    study = load_fasta(STUDY_ALN)
    n = len(sites)
    print(f'study alignment: {len(study)} seqs x {n} sites')

    # build index: site -> (study position in the alignment is column order)
    combined = {sm: list(seq) for sm, seq in study.items()}
    for sm in sorted(combined):
        assert len(combined[sm]) == n, f'{sm} length mismatch'

    # merge context genotypes at each site (allele-consistency by translation)
    context = {}
    for fn in sorted(os.listdir(GENO)):
        if not fn.endswith('.vcf'):
            continue
        sample = fn[:-4]
        ref2alt = {}
        with open(os.path.join(GENO, fn)) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                p = line.rstrip().split('\t')
                if len(p) < 5:
                    continue
                chrom, pos, refa, alta = p[0], int(p[1]), p[3], p[4].split(',')[0]
                gt = p[9].split(':')[0] if len(p) > 9 else '.'
                ref2alt[(chrom, pos)] = (refa, alta, gt)
        context[sample] = ref2alt
    print(f'context samples: {len(context)}')

    for sample, calls in context.items():
        seq = []
        for chrom, pos in sites:
            if (chrom, pos) in calls:
                refa, alta, gt = calls[(chrom, pos)]
                if len(refa) != 1 or len(alta) != 1:
                    seq.append('N')
                elif gt in ('0', '0/0', '0|0'):
                    seq.append(refa.upper())
                elif gt in ('1', '1/1', '1|1'):
                    seq.append(alta.upper())
                else:
                    seq.append('N')
            else:
                seq.append('N')
        if len(seq) != n:
            print(f'WARN: {sample} seq length {len(seq)} != {n}')
        combined[sample] = seq

    with open(COMBINED, 'w') as f:
        for sm in combined:
            f.write(f'>{sm}\n{"".join(combined[sm])}\n')
    print(f'combined alignment written: {len(combined)} seqs x {n} sites')

    # snp-dists + IQ-TREE2
    subprocess.run([os.path.join(ENV, 'snp-dists'), COMBINED],
                   stdout=open(os.path.join(OUT, 'combined_snpdists.tsv'), 'w'))
    print('snp-dists matrix written')
    subprocess.run([os.path.join(ENV, 'iqtree'), '-s', COMBINED, '-m', 'MFP', '-bb', '1000',
                    '-nt', 'AUTO', '--prefix', os.path.join(OUT, 'tree_combined')],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print('IQ-TREE done')


if __name__ == '__main__':
    main()
