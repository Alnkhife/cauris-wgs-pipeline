#!/usr/bin/env python3
"""Build a SNP-only multi-FASTA alignment from a multi-sample VCF (haploid).
Usage: wp5_build_aln.py <in.vcf> <ref.fasta> <out.fasta>
"""

import sys

vcf, ref, out = sys.argv[1], sys.argv[2], sys.argv[3]

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
    return {k: ''.join(v).upper() for k, v in seqs.items()}

genome = load_fasta(ref)
samples = []
sites = []
with open(vcf) as f:
    for line in f:
        if line.startswith('##'):
            continue
        p = line.rstrip().split('\t')
        if line.startswith('#'):
            samples = p[9:]
            continue
        chrom, pos, refa, alta = p[0], int(p[1]), p[3], p[4]
        if ',' in alta or '*' in alta:
            continue
        if len(refa) != 1 or len(alta) != 1:
            continue
        sites.append((chrom, pos, refa, alta, p[9:]))

aln = {sm: [] for sm in samples}
for chrom, pos, refa, alta, gts in sites:
    for sm, gt in zip(samples, gts):
        a = gt.split(':')[0]
        if a in ('0', '0/0', '0|0'):
            aln[sm].append(refa)
        elif a in ('1', '1/1', '1|1'):
            aln[sm].append(alta)
        else:
            aln[sm].append('N')

with open(out, 'w') as f:
    for sm in samples:
        f.write(f'>{sm}\n{"".join(aln[sm])}\n')
print(f'{len(sites)} sites, {len(samples)} samples -> {out}')
