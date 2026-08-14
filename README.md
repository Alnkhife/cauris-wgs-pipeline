# Candida auris WGS Analysis Pipeline

Independent whole-genome sequencing (WGS) analysis pipeline for a *Candida auris*
outbreak investigation


The pipeline performs reference-based variant confirmation, genome-wide SNP
calling, coverage-based copy-number screening, phylogenetic reconstruction, and
placement of study isolates within a national/regional genomic context. It was
built to cross-validate calls made with the CDC MycoSNP-NF pipeline and the
PathogenWatch web platform.

## Workflow

```
WP1  Alignment           bwa mem -> samtools fixmate/markdup -> indexed BAMs
WP2  Variant confirmation  bcftools mpileup/call at resistance loci
                          (ERG11, CDR1, FCA1/FCY1, GSC1/FKS1) + strand-aware
                          protein translation + PathogenWatch concordance
WP3  CNV / virulence screen  per-gene mean coverage (samtools bedcov) across
                          5,594 annotated genes; log2 normalised heatmap
WP4  Orthogonal SNP calling  freebayes (haploid) on ~80x downsampled BAMs;
                          1,301 biallelic SNPs after filtering
WP5  Phylogeny           SNP-only alignment -> snp-dists -> IQ-TREE2 (MFP,
                          1000 UFBoot) -> annotated tree figure
WP6  Regional context    align public Saudi C. auris genomes, genotype them at
                          the study SNP sites, combined 46-isolate phylogeny
```

## Tools and versions

| Tool | Version | Used for |
|---|---|---|
| bwa | 0.7.17 | read alignment to B8441 (GCF_002759435.1) |
| samtools | 1.21 | BAM processing, fixmate, markdup, bedcov, pileup |
| bcftools | 1.21 | variant calling, filtering, site genotyping |
| freebayes | 1.3.10 | independent genome-wide SNP calling |
| snp-dists | 1.2.0 | pairwise SNP distance matrices |
| IQ-TREE2 | 2.1.4 | maximum-likelihood phylogeny + ultrafast bootstrap |
| FastQC | 0.12.1 | read quality control |
| MultiQC | 1.35 | aggregated QC reporting |
| mosdepth | (optional) | per-base/per-window coverage |
| Python | 3.10 | analysis glue, figures (pandas, matplotlib, seaborn) |

## Reference

- *Candida auris* B8441, GCF_002759435.1 (Cand_auris_B8441_V3),
  7 chromosomes, 12.40 Mb, NCBI RefSeq, downloaded via the NCBI datasets API.
- Repeat-mask regions from the MycoSNP nucmer self-mask (5,849 regions,
  2.18 Mb) can be applied to the freebayes SNP set before distance
  calculation (`analysis/06_mycosnp/reference/masked_ref.bed` in the full
  analysis output).

## Usage

```bash
# 1. environment
conda env create -f environment.yml
conda activate cauris

# 2. reference index
bwa index -a bwtsw data/raw/reference/B8441_refseq.fna
samtools faidx data/raw/reference/B8441_refseq.fna

# 3. per-workpackage scripts (see scripts/ dir)
bash scripts/wp1_align.sh          # alignment
bash scripts/wp2_pileup.sh         # read-level variant confirmation
python scripts/wp2_annotate.py     # annotation + concordance tables
bash scripts/wp3_bedcov.sh         # per-gene coverage
python scripts/wp3_cnv.py          # CNV heatmap
bash scripts/wp4_freebayes.sh      # orthogonal SNP calling
python scripts/wp5_build_aln.py ...# SNP alignment
bash scripts/wp5_tree.sh           # distances + IQ-TREE2
bash scripts/wp6_align.sh          # context-genome alignment
bash scripts/wp6_genotype.sh       # site genotyping of context genomes
python scripts/wp6_build_combined.py  # combined alignment + tree
```

## Outputs

- `02_variant_confirm/` — per-locus pileup VCFs, key-site table, PathogenWatch
  concordance table (19/19 concordant)
- `03_virulence_cnv/` — per-gene coverage table and heatmap (no amplifications
  or deletions at resistance/virulence loci)
- `04_freebayes/` — raw and filtered VCFs, downsampled BAMs
- `05_phylogeny/` — SNP alignment, distance matrix, IQ-TREE2 tree, figures
- `07_context/` — combined 46-isolate alignment, matrix, and phylogeny

## Key results

- All 19 isolates: Clade I; ERG11 K143R (13/19), ERG11 Y132F (6/19),
  CDR1 V704L (13/19), FCY1 S70R (19/19); no FKS1 hot-spot mutations.
- 1,301 high-quality biallelic SNPs; three internally clonal lineages
  (within-lineage pairwise distances < 80 SNPs).
- Study isolates 0-1 SNPs from published Saudi outbreak isolates
  (Guan et al., 2025, *Microbiology Spectrum*).

## Notes

- Down-sampling to ~80x for freebayes: `samtools view -s` (seed 42),
  fraction computed from per-gene median coverage.
- IQ-TREE 3.1.3 (macOS ARM) rejects high-N alignments; IQ-TREE 2.1.4 is used
  for the combined alignment (container: quay.io/biocontainers/iqtree:2.1.4_beta).
- Full step-by-step methods and provenance are in `docs/methods.md`.
