# Methods

This document records the analysis steps, software versions, parameters, and
reference data used in this *Candida auris* whole-genome sequencing project.
Every step is reproducible from the scripts in `scripts/`.

## 1. Data and reference

- **Input**: 19 paired-end FASTQ files (Element Biosciences AVITI, 151 bp PE,
  8+8 dual indexing), demultiplexed with Bases2Fastq v1.8.0.
- **Reference**: *Candida auris* B8441, GCF_002759435.1 (Cand_auris_B8441_V3;
  7 chromosomes, 12.40 Mb), NCBI RefSeq.
- Read counts were verified against the run manifest (exact match for all
  samples).

## 2. Alignment (WP1)

- Tool: bwa 0.7.17 (`bwa mem`, 6 threads, per-sample read groups).
- Post-processing: `samtools sort -n` -> `samtools fixmate -m` ->
  `samtools sort` -> `samtools markdup -r` -> indexed BAM.
- Per-sample mean coverage (median over 5,594 genes): 143x-949x.

## 3. Read-level variant confirmation (WP2)

- Independent of pipeline calls: `bcftools mpileup` (Q20/q20, depth cap 500)
  + `bcftools call --ploidy 1` at four loci:
  ERG11 (NC_140807.1:1474722-1476296), CDR1 (NC_140808.1:348928-353454),
  FCY1/FCA1 (NC_140806.1, two exons), FKS1/GSC1 (NC_140808.1:2088186-2093852).
- Translation: exon-aware, strand-aware (VCF ALT complemented for
  minus-strand genes).
- Result: all 19 isolates concordant with the resistance calls
  (ERG11 K143R 13/19; ERG11 Y132F 6/19; CDR1 V704L 13/19; FCY1 S70R 19/19;
  no FKS1 hot-spot variants).

## 4. Copy-number / virulence screen (WP3)

- Per-gene mean coverage across all annotated genes via `samtools bedcov`;
  normalised as log2(gene/median) per sample.
- No amplification or deletion (ratio beyond +/-1) was observed at any
  resistance- or virulence-associated locus.

## 5. Orthogonal SNP calling (WP4)

- BAMs down-sampled to ~80x (`samtools view -s`, seed 42; fraction from
  per-gene median coverage).
- Caller: freebayes 1.3.10 (haploid; min base Q20, min map Q20, alt fraction
  0.5, alt count 3).
- Filters: biallelic SNPs, QUAL>=50, DP>=10, <10% missing -> 1,301 SNPs.
- Pairwise distances via snp-dists 1.2.0.

## 6. Phylogeny (WP5)

- SNP-only alignment (1,301 sites x 19 isolates) built from the filtered VCF.
- IQ-TREE2 2.1.4: ModelFinder (MFP) + 1,000 ultrafast bootstrap replicates.
- Note: IQ-TREE 3.1.3 (macOS ARM build) rejected the high-ambiguity
  alignment; IQ-TREE 2.1.4 was used instead.
- Three internally clonal lineages resolved (within-lineage medians 38-45
  SNPs, below the 80-SNP clonal threshold; between-lineage 102-233 SNPs).

## 7. Regional context (WP6)

- Public Saudi *C. auris* WGS runs downloaded from NCBI SRA / ENA:
  23 runs (ENA PRJEB84203; 2019 Riyadh outbreak isolates) and 4 runs
  (CDC global collection PRJNA595978).
- Context reads aligned to B8441 and genotyped at the 1,301 study SNP sites
  with `bcftools mpileup/call`.
- Combined alignment (46 isolates x 1,301 sites) -> snp-dists -> IQ-TREE2.
- Result: 26/27 Saudi context isolates within 0-1 SNPs of a study isolate;
  the study isolates belong to the same regional clonal outbreak complex.

## 8. Environment

- Conda environment (see `environment.yml`): bwa 0.7.17, samtools 1.21,
  bcftools 1.21, freebayes 1.3.10, snp-dists 1.2.0, FastQC 0.12.1,
  MultiQC 1.35, Python 3.10 (numpy, pandas, scipy, matplotlib, seaborn).
- Reference mask (nucmer self-mask, 5,849 regions / 2.18 Mb) available for
  re-filtering the SNP set before distance calculation if desired.
