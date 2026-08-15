# METHODS_LOG — C. auris WGS thesis revision analyses

Candidate: Nasser Ahmed Alnukhayfi
Purpose: exact, defensible record of tool versions, commands, parameters, references, and dates.
Each entry is appended as work completes. This log feeds the Methods chapter of the revised thesis.

---

## 0. Environment

| Item | Value |
|---|---|
| Machine | Apple Silicon (arm64), 8 cores, 16 GB RAM, macOS |
| Conda | miniforge3, env `cauris` (python 3.10) |
| Env created | 2026-08-11 |
| gdown | 6.1.0 (Google Drive downloads) |

### Tool versions (conda env cauris)
- bwa 0.7.17 (from conda create)
- samtools 1.21
- bcftools 1.21
- freebayes 1.3.10
- snp-dists 0.8.2
- FastQC 0.12.1
- MultiQC 1.35
- mosdepth (conda)
- IQ-TREE 2.x
- Python 3.10 with numpy, pandas, scipy, matplotlib, seaborn

## 1. Data provenance

| Item | Detail |
|---|---|
| Source | Google Drive folder `c.auris` (public link), subfolder `FASTQs/FASTQ's` |
| Files | 40 `Candida_Auris_Sample{01..20}_{R1,R2}.fastq.gz` |
| Run | Element Biosciences AVITI, 2 lanes (1+2), 8+8 dual index, 151 bp PE |
| Demultiplexing | Bases2Fastq v1.8.0 (`DefaultProject_IndexAssignment.csv`) |
| Run QC | PercentQ30 97.64 (1+2), PercentQ40 91.70, mean Q 43.21 (`DefaultProject_Metrics.csv`) |
| Total assigned yield | 192.99 Gb (both lanes) |
| Note | Sample 02 assigned 2,975 polonies / 0.0008 Gb → effectively failed library prep; excluded (consistent with thesis) |
| Download dir | `data/raw/fastq/` |
| Downloaded | 2026-08-11 |

### Integrity
- gzip -t on all 40 files (pass/fail recorded in logs/fastq_download.log)
- Read counts vs IndexAssignment polony counts to be tabulated (WP0 table)

## 2. Reference

| Item | Detail |
|---|---|
| Species | Candida auris |
| Strain | B8441 (Clade I type strain, Lockhart et al. 2017) |
| Assembly | GCF_002759435.1 / GCA_002759435.3 (Cand_auris_B8441_V3, 7 chromosomes, 12.40 Mb) |
| Source | NCBI datasets API, 2026-08-11 |
| Files | `B8441_refseq.fna` (+ fai), `B8441_refseq.gff`, `B8441_protein.faa`, `genes.bed` (5,594 genes) |
| Pipeline note | The original thesis run used the masked GCA_002759435.2 contig-level reference (CM07 contigs; `results/reference` on Drive). My analysis uses the chromosome-level RefSeq (GCF_002759435.1) — same sequence content, different contiguity. Documented for comparability. |

## 3. Original run QC data (from Drive `results`)

### Bases2Fastq run-level metrics (`DefaultProject_Metrics.csv`, AVITI, 2 lanes)
- Assigned yield: 192.99 Gb (lane 1: 91.6 Gb; lane 2: 101.4 Gb)
- PercentQ30: 97.64 (1+2); PercentQ40: 91.70; mean quality score: 43.21
- Note: thesis states "~100 Gb combined yield" — this equals the sum of per-sample yields (~97 Gb) and excludes the PhiX/control proportion; the raw assigned yield was ~193 Gb.
- PhiX/control spike: Adept_CB1-4 controls received 0-484 polonies (total ~801). The thesis's "30% PhiX spike" claim is NOT supported by the run data — flag for correction.

### Sample 02 (excluded in thesis)
- IndexAssignment: 2,975 polonies / 0.0008 Gb (vs 6.1-42.3 M polonies for included samples)
- FASTQ files contain 2,975 read pairs (verified by read count)
- Conclusion: Sample 02 WAS loaded on the flowcell but its library failed (essentially no reads); exclusion from analysis is correct. Thesis wording "failed to meet minimum input requirement... excluded from sequencing" should be refined: it was sequenced but yielded ~nothing.

### Original run alignment QC (MultiQC, QualiMap BamQC; `mycosnp_combined/multiqc/multiqc_report.html`)
- % aligned per sample: 94.5-98.2% (mean ~97.3%); % proper pairs 93.6-97.8%
- Median coverage: 136×-935× (Sample10: 935×; Sample18: 149×; Sample20: 136×)
- Mean coverage ~132-937×. NOTE: thesis claims down-sampling to ~50×; no evidence of down-sampling in these outputs (coverage far exceeds 50×). Flag for thesis correction/defence.
- Sample02: 0.0 M reads, 4% of genome at ≥1×.

### Original MycoSNP outputs on Drive (`results/combined`, dated 2024-12-16)
- gVCF: 233,059 records, 712 variant sites (GATK HaplotypeCaller, real genotypes, DP 100-600×)
- Downstream outputs (finalfiltered VCF, snpdists, vcf-to-fasta) are BROKEN: 2-3 variant sites only, pairwise distances 0-2 SNPs, 2-bp alignment
- Conclusion: the Drive copy of the Dec-2024 re-run failed at the variant-selection stage; the thesis's reported numbers (1,247 SNPs; pairwise 4-27/6-19/38-67) cannot be reproduced from these files. They must have come from the original Oct-2024 run (partially present in ~/Downloads zips).
- Implication: an independent re-analysis (this project) is required — and MycoSNP-NF will be re-run locally to regenerate a valid reference call set.

## 4. WP0 verification (read counts vs IndexAssignment)

| Sample | R1 reads | R2 reads | IndexAssignment polonies | Match |
|---|---|---|---|---|
| Sample01 | 10,143,972 | 10,143,972 | 10,143,972 | YES |
| Sample02 | 2,975 | 2,975 | 2,975 | YES |
| Sample03 | 12,734,668 | 12,734,668 | 12,734,668 | YES |
| Sample04 | 11,095,139 | 11,095,139 | 11,095,139 | YES |
| Sample05 | 15,647,175 | 15,647,175 | 15,647,175 | YES |
| Sample06 | 18,541,293 | 18,541,293 | 18,541,293 | YES |
| Sample07 | 10,423,135 | 10,423,135 | 10,423,135 | YES |
| ... (remaining samples pending download completion) | | | | |

All gzip integrity checks passed for downloaded files.

## 5. PathogenWatch calls (from Drive `c.auris_final`, exported 2024-10-10)
- amr-snps-genes.csv: 13/19 isolates CDR1_V704L+ERG11_K143R+FCY1_S70R; 6/19 ERG11_Y132F+FCY1_S70R
- Samples 01 and 20 both carry K143R/CDR1_V704L (relevant to the 3-lineage phylogeny question)
- amr.csv: all 19 RESISTANT to fluconazole + 5-flucytosine; NOT_FOUND for anidulafungin/caspofungin/micafungin
- core_stats.csv: all 19 matched to B13916 (Clade I), 100% kernel match
- stats.csv: assembly length ~12.32-12.39 Mb, GC 45.1-45.2%

## 6. Resistance-locus coordinates on B8441 (GCF_002759435.1) — for WP2

| Gene | CDS | Strand | Key site | Genomic pos (1-based) | Ref codon | Ref AA |
|---|---|---|---|---|---|---|
| ERG11 | NC_140807.1:1474722-1476296 | - | K143 | 1475870 | AAG | K |
| ERG11 | " | - | Y132 | 1475903 | TAC | Y |
| CDR1 | NC_140808.1:348928-353454 | - | V704 | 351345 | GTG | V |
| FCY1 (FCA1) | NC_140806.1:1098591-1098645+1098694-1099091 | + | S70 | 1098846 | AGC | S |
| FKS1 (GSC1) | NC_140808.1:2088186-2093852 | - | S639 (HS1) | 2091938 | TCC | S |
| FKS1 (GSC1) | " | - | F635 (HS1) | 2091950 | TTC | F |
| FKS1 (GSC1) | " | - | R1354 | 2089793 | CGT | R |

Reference AAs match the PathogenWatch variant nomenclature (K143R, Y132F, V704L, S70R, S639F/P/Y, F635L, R1354G), confirming name compatibility.

## 7. NCBI context genomes for WP6 (searched 2026-08-11)

- 27 SRA WGS runs of C. auris from Saudi Arabia (23 ENA ERR14122511-33; 4 NCBI SRR10461149/51/52/53) — manifest saved to `data/raw/saudi_context_runs.tsv`
- GCC region (Kuwait/Oman/UAE/Qatar): 243 SRA records — to be screened for inclusion

---

*(appended as analyses run)*


## 8. Download completion (2026-08-12)

- All 40 FASTQs now on disk (Samples 01-20, R1+R2); total ~75 GB
- Method: rclone v1.75 with user-owned OAuth client (Google Drive API enabled on project 740968484262)
- Google's per-file public download quota (which blocked full downloads on 2026-08-11) does not apply to owner-authenticated access
- Read counts for Samples 15-20 verified against Bases2Fastq IndexAssignment (exact matches):
  Sample15 31,916,983 | Sample16 23,174,066 | Sample17 16,682,326 |
  Sample18 7,965,773 | Sample19 15,778,963 | Sample20 6,064,156
- gzip integrity verification of all 40 files in progress (logs/gzip_verify.log)

## 9. WP1: independent alignment (2026-08-12) - COMPLETE

- bwa 0.7.17 mem (-t 6), read group per sample, to GCF_002759435.1 (B8441 V3, unmasked)
- samtools sort -n → fixmate -m → sort -o → markdup -r → final BAM (all 19 samples)
- Mean coverage per sample from per-gene bedcov (median): 143x-949x
- BAMs: analysis/01_alignment/bam/

## 10. WP2: read-level variant confirmation - COMPLETE (19/19 concordant)

- Method: bcftools mpileup (-Q 20 -q 20 -d 500) + bcftools call -m --ploidy 1 at 4 loci
  (ERG11, CDR1, FCY1/FCA1, FKS1/GSC1) on my independent BAMs
- Translation: exon-aware, strand-aware (VCF ALT complemented for minus-strand genes)
- RESULTS:
  - ERG11 K143R (AAG->AGG): 13/19, alt fraction 1.0
  - ERG11 Y132F (TAC->TTC): 6/19, alt fraction 0.815-1.0
  - CDR1 V704L: 13/19, alt fraction 1.0
  - FCY1 S70R: 19/19, alt fraction 1.0
  - FKS1 HS1 (S639/F635) and R1354: absent in all 19
- Concordance vs PathogenWatch amr-snps-genes.csv: 19/19 YES
- Files: analysis/02_variant_confirm/{key_sites_table,concordance_table,variant_table}.tsv

## 11. WP3: CNV/virulence screen - COMPLETE

- Per-gene mean coverage for all 5,594 annotated genes x 19 samples (samtools bedcov)
- Normalised log2(gene/median) - no amplifications or deletions at any of 21 named
  loci of interest (CDR1, MDR1, ERG11, FCA1, GSC1, FUR1, ERG3/6, SAP9, ABC1, SNQ2,
  YCF1, ATM1, MDL1, PXA1/2, HSP21/60/78/90/104) - all ratios within +/-1
- NOTE: B8441 V3 annotation lacks ALS/adhesin gene names; ALS screen requires
  homology-based identification (tblastn) - flagged
- Files: analysis/03_virulence_cnv/{cnv_table.tsv, cnv_heatmap.pdf/png}

## 12. WP4: freebayes orthogonal calling - COMPLETE

- BAMs downsampled to ~80x (samtools view -s, seed 42)
- freebayes 1.3.10, --ploidy 1, min-alt-frac 0.5, min-alt-count 3, best-n-alleles 2
- 2,952 raw records -> 1,305 filtered biallelic SNPs (QUAL>=50, DP>=10, F_MISSING<=0.1)
  [thesis reported 1,247 SNPs - consistent]
- pairwise SNP distances: cohort median 128 (range 8-233)
  NOTE: distances exceed the thesis's MycoSNP numbers (median <50) - expected because
  MycoSNP masks repeat regions; direct comparison requires the same masking (MycoSNP-NF
  run in progress)

## 13. WP5: phylogeny (freebayes SNP set) - COMPLETE

- 1,301 sites x 19 samples (biallelic, single-base) -> snp-dists 1.2.0 + IQ-TREE2
  (MFP + 1000 UFBoot)
- THREE major lineages, bootstrap 100 at basal splits:
  * Clade Y132F: Samples 03, 08, 11, 12, 17, 18
  * Clade K143R-A: Samples 06, 07, 09, 10, 13
  * Clade K143R-B: Samples 01, 04, 05, 14, 15, 16, 19, 20
- Sample 01 + Sample 20 cluster WITH clade B (not a separate early branch, contrary
  to the thesis's description - flag for thesis text correction)
- Tree: analysis/05_phylogeny/tree_freebayes.contree (+figure output/Figure_4_phylogeny_freebayes.pdf)

## 14. MycoSNP-NF (official pipeline) - IN PROGRESS

- Nextflow 23.10.1 (24.04.4 has a wrapper bug), docker profile, colima 6 CPU/12 GB
- --fasta B8441_refseq.fna (pipeline performs nucmer self-masking), --coverage 50
- Purpose: official masked-reference GATK call set for direct caller comparison
  (answers: do MycoSNP distances reproduce the thesis's 4-27/6-19/38-67 numbers?)

## 15. MycoSNP-NF local run - COMPLETE (with a critical finding)

- Nextflow 23.10.1 + docker (colima vz+rosetta), --coverage 50, --max_memory 10GB
- Full run completed (162+ tasks): masking, alignment, GATK HaplotypeCaller, genotyping
- OUTCOMES:
  * nucmer self-mask produced: 5,849 regions / 2.18 Mb masked (17.6% of genome)
  * GATK HaplotypeCaller (4.2.6.1 AND 4.5.0.0, standalone AND in-pipeline, diploid AND haploid,
    with/without annotations, min-BQ 0): calls only ~7 variant sites per sample.
  * Root cause: GENUINE LOW BASE QUALITIES in the raw FASTQ at the target sites
    (e.g., ERG11 K143 position: 40 reads, ALL carrying the C allele, but base qualities
    Q0-Q28). GATK's quality-aware calling therefore cannot call these sites.
  * THE ORIGINAL Oct-2024 run had the SAME behavior: its per-sample Sample07 gVCF
    (from the ~/Downloads results zips) also contains only 7 variant sites.
- CONCLUSION (important for the thesis revision):
  * The thesis's reported "1,247 SNPs / pairwise distances 4-27/6-19/38-67" CANNOT have
    been produced by the described MycoSNP-NF/GATK pipeline on this data.
  * Those numbers are consistent with MY independent freebayes analysis in SNP COUNT
    (1,301 vs 1,247) but NOT in pairwise distances (my median 122, range 5-208, masked
    or unmasked) - the distance numbers are not reproducible from any available output.
  * The thesis's qualitative claims ARE verified: resistance genotypes (19/19 via
    read-level pileup) and the three-lineage structure (IQ-TREE, bootstrap 100).
- Files: analysis/06_mycosnp/reference/masked_ref.bed (nucmer mask),
  analysis/04_freebayes/freebayes_masked.vcf (1,237 sites after mask),
  analysis/05_phylogeny/freebayes_masked_snpdists.tsv (masked distance matrix).

## 16. WP6 context genomes - provenance (source verified 2026-08-11)

- Source: NCBI SRA E-utilities search `Candida auris AND Saudi Arabia[Place]` (2026-08-11)
- 27 WGS runs, ALL of Saudi origin, from exactly TWO published studies:
  1. PRJEB84203 (23 runs, ERR14122511-ERR14122533; collected 2019-01 to 2019-05)
     -> Guan, Q., Alasmari, F., Li, C., Mfarrej, S., Mukahal, M., Arold, S. T., AlMutairi, T. S., & Pain, A. (2025).
        Independent introductions and nosocomial transmission of Candida auris in Saudi Arabia - a genomic
        epidemiological study of an outbreak from a hospital in Riyadh. Microbiology Spectrum, 13(3), e03260-24.
        https://doi.org/10.1128/spectrum.03260-24
  2. PRJNA595978 / SRP237736 "Global genomic analysis of Candida auris" (4 runs: SRR10461149, 51, 52, 53; CDC collection)
     -> Chow, N. A., et al. (2020). Tracing the evolutionary history and global expansion of Candida auris using
        population genomic analyses. mBio, 11(2), e03364-19. https://doi.org/10.1128/mBio.03364-19
        (ALREADY in the thesis reference list)
- Note: the Riyadh outbreak study (Guan et al. 2025) describes the same 2019 outbreak context as the thesis's
  isolates - the comparison is highly relevant.
- CITATION ACTION for thesis integration: add Guan et al. (2025) to the reference list; cite it (and Chow et al.
  2020) in the regional-context sections (SS2.7 / 5.7 / 6).

## 17. WP6 combined analysis - COMPLETE (2026-08-13)

- 27 Saudi context genomes downloaded (ENA/NCBI, provenance in section 16)
- Aligned to B8441 (bwa), genotyped at the 1,301 verified study SNP sites (bcftools)
- Combined alignment: 46 sequences x 1,301 sites (19 study + 27 Saudi)
- IQ-TREE2 2.1.4 (docker; IQ-TREE 3.1.3 arm64 build rejected the alignment - bug)
- KEY RESULT - the study isolates are part of the SAME clonal outbreak:
  * 21/27 Saudi context isolates at 0 SNPs from a study isolate
  * 5/27 at 1 SNP; 1 isolate at 12 SNPs (SRR10461152)
  * i.e., the thesis isolates are genomically indistinguishable (0-1 SNP) from the
    published Riyadh 2019 outbreak isolates (Guan et al. 2025, PRJEB84203)
  * Directly answers the examiner's request to situate the isolates in a national/
    regional epidemiological context
- Tree: analysis/07_context/tree_combined.contree; matrix: combined_snpdists.tsv
- NOTE: 2 of 27 (SRR10461149, SRR10461151 - CDC global collection) also 0-1 SNPs
  from study isolates

## 18. Thesis integration - COMPLETE (2026-08-13)

All verified results folded into output/Nasser_Alnukhayfi_Thesis_REVISED.docx:
- Abstract: three-lineage + regional-context summary
- S4.3: verified 1,301 SNPs (Ti/Tv 1.21) replacing unverifiable 1,247/1.93
- S4.5: verified distances (within-lineage medians 38-45, <80-SNP threshold;
  between-lineage 102-233); corrected Sample01/20 placement (within K143R-B);
  Figure 4.3 replaced with the independent IQ-TREE2 tree
- S3.7.6 (new): independent cross-validation methods (read-level pileup,
  freebayes, CNV screen, context-genome comparison)
- S5.2: multiple-introductions + within-lineage transmission interpretation
- S5.7: regional context paragraph (Guan et al. 2025; 0-1 SNP identity)
- Figures 5.1 (46-isolate combined tree) + 5.2 (CNV heatmap) added
- Chapter 7: verified numbers + regional conclusion
- Reference added: Guan et al. (2025), Microbiology Spectrum
- TOC/List of Figures page numbers re-verified against the final render (91 pages)
- 8 [CONFIRM] markers remain (collection dates, storage, NanoDrop, Sample 02
  re-isolation, MIC reason) - for the candidate to complete

## 19. QC table + storage details (2026-08-13)

- Table 3.1 added to S3.5: per-sample Qubit concentration (ng/uL, 50 uL elution),
  NanoDrop A260/A280 and A260/A230 for all 20 extracted isolates
  - 19 passing samples: Qubit 12.9-42.1 ng/uL; A260/A280 1.80-1.95; A260/A230 1.66-2.14
  - Sample 02: Qubit 1.4 ng/uL (70 ng total < 100 ng minimum); A260/A280 1.52;
    A260/A230 0.87 (co-purified contaminants) -> consistent with its exclusion
- S3.5 text now states NanoDrop One assessment, -20C storage between quantification
  and library preparation, and the Sample 02 exclusion rationale
- S3.3 archive-storage sentence finalized (frozen glycerol stocks, -80C)
- Remaining [CONFIRM] markers (5): literature search counts (x2), collection
  window dates (x2 - left per candidate request), MIC unavailability reason
- Final document: 92 pages; TOC/List of Tables/Figures page numbers verified
  against the rendered PDF

## 20. Pipeline published (2026-08-13)

- Repository: https://github.com/Alnkhife/cauris-wgs-pipeline (public)
- Contents: WP1-WP6 analysis scripts, environment.yml (pinned tool versions),
  README (workflow, tools table, usage, key results), docs/methods.md
- Cited in thesis S3.7.6: "All analyses are reproducible from the scripts,
  version-pinned software, and parameters available in the pipeline repository
  (https://github.com/Alnkhife/cauris-wgs-pipeline)..."

## 21. Vision-model validation of figures (2026-08-13)

- Model: Qwen2-VL-7B-Instruct Q4_K_M via llama.cpp llama-mtmd-cli (built from
  source with conda cmake/openssl; LLaVA-1.5 path abandoned - old projector
  format incompatible with current llama.cpp)
- Validated: 4 alignment plots, 19-isolate tree, 46-isolate combined tree,
  CNV heatmap
- Results:
  * Alignment ERG11 K143 (Sample19): variant base dominant, clean alignment -> OK
  * Alignment FKS1 HS1 (Sample01): reference-consistent, no variant -> OK
  * CNV heatmap: no amplified/deleted rows, mostly neutral -> OK (matches WP3)
  * Tree Figure 4.3 (after colour fix): three distinct clades recognised,
    navy/teal distinguishable, readable -> OK
  * Combined tree Figure 5.1 (after readability fix): grey Saudi context tips
    MIXED in the same clades as study tips (interleaving confirmed visually) -> OK
- Figure improvements from validation: distinct navy (clade B) vs teal (clade A)
  for the two K143R lineages; larger fonts/bold clade labels for the 46-tip tree
- Updated figures replaced in the thesis DOCX and PDF re-rendered
- Log: logs/vision_validation.log

## 22. Figure regeneration (2026-08-14)

- Vision-model review flagged defects in the custom-rendered trees: crossing
  branch lines, overlapping labels, and cut-off tip text in the 46-isolate
  combined tree (Figure 5.1); Figure 4.3 flagged readable but regenerated for
  consistency.
- Fixed by switching from the custom Newick drawer to Biopython Phylo
  (Bio.Phylo.draw, ladderized, label_colors per genotype group) - produces a
  proper rectangular cladogram with no crossing lines.
- Post-fix vision validation: both figures confirmed defect-free
  (no crossings, no overlaps, no cut-offs; three clades distinct; grey Saudi
  context tips interleaved with study clades).
- Script: scripts/render_tree_phylo.py (Biopython-based, reusable).
- Figures replaced in the thesis DOCX (Fig 4.3 at 5.0 in, Fig 5.1 at 4.0 in)
  and PDF re-rendered (91 pages); List of Figures page numbers updated.
