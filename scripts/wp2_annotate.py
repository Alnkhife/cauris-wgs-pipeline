#!/usr/bin/env python3
"""
WP2 stage C: annotate pileup VCFs with amino-acid changes and build
the concordance table against PathogenWatch (amr-snps-genes.csv).

Input : analysis/02_variant_confirm/vcfs/<gene>.raw.vcf
Output: analysis/02_variant_confirm/variant_table.tsv
        analysis/02_variant_confirm/concordance_table.tsv
"""
import csv
import os
import re

ROOT = os.path.expanduser('~/Projects/cauris-wgs-analysis')
REF = os.path.join(ROOT, 'data/raw/reference/B8441_refseq.fna')
GFF = os.path.join(ROOT, 'data/raw/reference/B8441_refseq.gff')
OUT = os.path.join(ROOT, 'analysis/02_variant_confirm')
PW = os.path.join(ROOT, 'data/raw/auris_final/amr-snps-genes.csv')

SAMPLES = ['Sample01'] + [f'Sample{i:02d}' for i in range(3, 21)]

# gene -> CDS exons (from B8441 V3 RefSeq)
GENES = {
    'ERG11': ('NC_140807.1', [(1474722, 1476296, '-')]),
    'CDR1':  ('NC_140808.1', [(348928, 353454, '-')]),
    'FCY1':  ('NC_140806.1', [(1098591, 1098645, '+'), (1098694, 1099091, '+')]),
    'FKS1':  ('NC_140808.1', [(2088186, 2093852, '-')]),
}

CODON = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S',
    'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C',
    'TGA': '*', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P',
    'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
    'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V',
    'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
    'GGG': 'G',
}
COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def revcomp(s):
    return s.translate(COMP)[::-1]


def load_fasta(path):
    seqs, cur = {}, None
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                cur = line[1:].split()[0]
                seqs[cur] = []
            else:
                seqs[cur].append(line.rstrip())
    return {k: ''.join(v).upper() for k, v in seqs.items()}


def cds_for_gene(genome, chrom, exons):
    order = sorted(exons, key=lambda x: -x[0]) if exons[0][2] == '-' else sorted(exons, key=lambda x: x[0])
    cds = ''.join(genome[chrom][s - 1:e] for s, e, _ in order)
    if exons[0][2] == '-':
        cds = revcomp(cds)
    return cds


def pos_to_codon_idx(exons, chrom, pos, genome):
    """Map a genomic position to (cds index)."""
    strand = exons[0][2]
    for s, e, _ in exons:
        if s <= pos <= e:
            if strand == '-':
                offset = e - pos
            else:
                offset = pos - s
            # accumulate previous exons' lengths
            prev = 0
            order = sorted(exons, key=lambda x: -x[0]) if strand == '-' else sorted(exons, key=lambda x: x[0])
            for ps, pe, _ in order:
                if (ps, pe) == (s, e):
                    break
                prev += pe - ps + 1
            return prev + offset
    return None


def translate_with_var(cds, idx, ref, alt, strand):
    if alt.startswith('<') or ',' in alt:
        return None, None
    alt = alt.split(',')[0]
    seq = list(cds)
    alt_b = revcomp(alt) if strand == '-' else alt
    seq[idx] = alt_b
    seq = ''.join(seq)
    prot_ref = ''.join(CODON[cds[i:i + 3]] for i in range(0, len(cds) - 2, 3))
    prot_alt = ''.join(CODON[seq[i:i + 3]] for i in range(0, len(seq) - 2, 3))
    pos = idx // 3 + 1
    a_ref = prot_ref[pos - 1] if pos <= len(prot_ref) else '?'
    a_alt = prot_alt[pos - 1] if pos <= len(prot_alt) else '?'
    if a_ref == a_alt:
        return None, None
    return f'{a_ref}{pos}{a_alt}', pos


def parse_vcf(path):
    """Return list of (chrom, pos, ref, alt, {sample: (AD_ref, AD_alt, DP)})."""
    records = []
    if not os.path.exists(path):
        return records
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.rstrip().split('\t')
            if len(p) < 10:
                continue
            chrom, pos, ref, alt = p[0], int(p[1]), p[3], p[4]
            if 'DUP' in alt or '*' in alt:
                continue
            fmt = p[8].split(':')
            ad_idx = fmt.index('AD') if 'AD' in fmt else None
            dp_idx = fmt.index('DP') if 'DP' in fmt else None
            samples = {}
            for sm, gt_field in zip(SAMPLES, p[9:]):
                vals = gt_field.split(':')
                adr = ada = dp = 0
                if ad_idx is not None and ad_idx < len(vals):
                    ad = vals[ad_idx].split(',')
                    if len(ad) >= 2:
                        try:
                            adr, ada = int(ad[0]), int(ad[1])
                        except ValueError:
                            pass
                if dp_idx is not None and dp_idx < len(vals):
                    try:
                        dp = int(vals[dp_idx])
                    except ValueError:
                        pass
                samples[sm] = (adr, ada, dp)
            records.append((chrom, pos, ref, alt, samples))
    return records


def main():
    genome = load_fasta(REF)
    all_rows = []
    for gene, (chrom, exons) in GENES.items():
        cds = cds_for_gene(genome, chrom, exons)
        vcf = os.path.join(OUT, 'vcfs', f'{gene}.raw.vcf')
        for chrom_v, pos, ref, alt, samples in parse_vcf(vcf):
            idx = pos_to_codon_idx(exons, chrom, pos, genome)
            if idx is None:
                continue
            if idx >= len(cds):
                continue
            aa_change, aapos = translate_with_var(cds, idx, ref, alt, exons[0][2])
            for sm, (adr, ada, dp) in samples.items():
                if dp == 0:
                    continue
                af = ada / (adr + ada) if (adr + ada) else 0
                all_rows.append({
                    'gene': gene, 'chrom': chrom_v, 'pos': pos, 'ref': ref, 'alt': alt,
                    'aa_change': aa_change or 'synonymous', 'aa_pos': aapos or '',
                    'sample': sm, 'ref_depth': adr, 'alt_depth': ada, 'depth': dp,
                    'alt_fraction': round(af, 3),
                })
    with open(os.path.join(OUT, 'variant_table.tsv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()) if all_rows else ['gene'],
                           delimiter='\t')
        w.writeheader()
        w.writerows(all_rows)
    print(f'variant_table.tsv: {len(all_rows)} rows')

    # ---- key sites summary (the 7 reported sites) ----
    key_sites = {
        'K143R': 'K143R', 'Y132F': 'Y132F',
        'V704L': 'V704L', 'S70R': 'S70R',
        'S639F': 'S639F/P/Y', 'S639P': 'S639F/P/Y', 'S639Y': 'S639F/P/Y',
        'F635L': 'F635L', 'R1354G': 'R1354G',
    }
    summary = {sm: {} for sm in SAMPLES}
    for r in all_rows:
        if r['aa_change'] in key_sites and r['alt_fraction'] >= 0.5:
            summary[r['sample']][key_sites[r['aa_change']]] = f"{r['alt']}({r['alt_fraction']})"
    with open(os.path.join(OUT, 'key_sites_table.tsv'), 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['sample'] + list(dict.fromkeys(key_sites.values())))
        for sm in SAMPLES:
            w.writerow([sm] + [summary[sm].get(k, '-') for k in dict.fromkeys(key_sites.values())])
    print('key_sites_table.tsv written')

    # ---- concordance with PathogenWatch ----
    pw_calls = {}
    if os.path.exists(PW):
        with open(PW) as f:
            for row in csv.DictReader(f):
                name = row['Genome Name'].replace('Candida_Auris_', '').replace('_contigs', '')
                sm = 'Sample' + name[-2:] if name[-2:].isdigit() else name
                pw_calls[sm] = row['Variants']
    with open(os.path.join(OUT, 'concordance_table.tsv'), 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['sample', 'pathogenwatch_variants', 'readlevel_variants', 'concordant'])
        for sm in SAMPLES:
            pv = pw_calls.get(sm, '')
            rv = '+'.join(f"{k}:{v}" for k, v in sorted(summary[sm].items()))
            # normalize names for comparison
            pw_set = set(re.findall(r'[A-Z]+[0-9]+[A-Z]+[0-9]*', pv))
            rv_set = set(re.findall(r'[A-Z]+[0-9]+[A-Z]+[0-9]*', rv))
            concordant = 'YES' if pw_set == rv_set else 'PARTIAL' if pw_set & rv_set else 'NO'
            w.writerow([sm, pv, rv, concordant])
    print('concordance_table.tsv written')


if __name__ == '__main__':
    main()
