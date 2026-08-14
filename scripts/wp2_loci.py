#!/usr/bin/env python3
"""
WP2 stage A: locus table for resistance genes on B8441 (GCF_002759435.1).

For each gene: merge CDS exons, translate, and map key codon positions back to
genomic coordinates (exon-aware). Output: analysis/02_variant_confirm/locus_coordinates.tsv
"""


REF = '$HOME/Projects/cauris-wgs-analysis/data/raw/reference/B8441_refseq.fna'
GFF = '$HOME/Projects/cauris-wgs-analysis/data/raw/reference/B8441_refseq.gff'
OUT = '$HOME/Projects/cauris-wgs-analysis/analysis/02_variant_confirm'

LOCI = {
    'ERG11': ('NC_140807.1', '-', 'lanosterol 14-alpha-demethylase', [('K143', 143), ('Y132', 132)]),
    'CDR1':  ('NC_140808.1', '-', 'ABC efflux transporter', [('V704', 704)]),
    'FCY1':  ('NC_140806.1', '+', 'cytosine deaminase (FCA1)', [('S70', 70)]),
    'FKS1':  ('NC_140808.1', '-', '1,3-beta-glucan synthase (GSC1)', [('S639', 639), ('F635', 635), ('R1354', 1354)]),
}
ALIAS = {'FCY1': 'FCA1', 'FKS1': 'GSC1'}

CODON_TABLE = {
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
            line = line.rstrip()
            if line.startswith('>'):
                cur = line[1:].split()[0]
                seqs[cur] = []
            else:
                seqs[cur].append(line)
    return {k: ''.join(v).upper() for k, v in seqs.items()}


def load_cds_exons(gff, gene_name):
    exons, chrom, strand = [], None, None
    with open(gff) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip().split('\t')
            if len(parts) < 9 or parts[2] != 'CDS':
                continue
            if f'gene={gene_name}' in attrs_check(parts[8]):
                exons.append((int(parts[3]), int(parts[4]), parts[6]))
                chrom, strand = parts[0], parts[6]
    return (chrom, exons, strand) if exons else None


def attrs_check(s):
    return s


def transcript_and_map(cds_seq, strand, exons):
    """Return (translated_protein, list of (codon_idx0, genomic_pos1_of_codon_start))."""
    if strand == '-':
        order = sorted(exons, key=lambda x: -x[0])
    else:
        order = sorted(exons, key=lambda x: x[0])
    # walk exons in transcription order, tracking genomic->transcript offset
    codon_pos = {}
    offset = 0  # transcript base offset (in CDS)
    for s, e, st in order:
        length = e - s + 1
        for i in range(length):
            if offset % 3 == 0:
                if strand == '-':
                    codon_pos[offset // 3] = e - i  # genomic pos of transcript base `offset`
                else:
                    codon_pos[offset // 3] = s + i
            offset += 1
    prot = ''.join(CODON_TABLE[cds_seq[i:i + 3]] for i in range(0, len(cds_seq) - 2, 3))
    return prot, codon_pos


def main():
    os.makedirs(OUT, exist_ok=True)
    genome = load_fasta(REF)
    rows = []
    with open(os.path.join(OUT, 'locus_coordinates.tsv'), 'w') as out:
        out.write('gene\tvariant\tchrom\tpos1\tref_codon\tref_aa\tcds_len\n')
        for gene, (chrom, strand, prod, sites) in LOCI.items():
            feat = load_cds_exons(GFF, ALIAS.get(gene, gene))
            if feat is None:
                print(f'{gene}: CDS NOT FOUND')
                continue
            fchrom, exons, fstrand = feat
            if fchrom != chrom:
                print(f'{gene}: GFF chrom {fchrom} != expected {chrom}')
            # assemble CDS (reverse complement for minus strand)
            order = sorted(exons, key=lambda x: -x[0]) if fstrand == '-' else sorted(exons, key=lambda x: x[0])
            cds = ''.join(genome[fchrom][s - 1:e] for s, e, _ in order)
            if fstrand == '-':
                cds = revcomp(cds)
            prot, codon_pos = transcript_and_map(cds, fstrand, order)
            for label, pos in sites:
                if pos > len(prot):
                    print(f'{gene} {label}: pos {pos} > len {len(prot)}')
                    continue
                idx = pos - 1
                codon = cds[idx * 3:idx * 3 + 3]
                gpos = codon_pos.get(idx)
                row = (gene, label, fchrom, gpos, codon, prot[pos - 1], len(cds))
                rows.append(row)
                out.write('\t'.join(map(str, row)) + '\n')
                print(f'{gene} {label}: chr{fchrom}:{gpos} codon={codon} ref_aa={prot[pos - 1]} (len {len(prot)} aa)')


if __name__ == '__main__':
    main()
