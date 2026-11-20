# brca-sv-amp-timing-nf

Nextflow pipeline for timing structural variants and focal amplifications in the BRCA SASH cohort. It runs
on Seqera/AWS Batch.

- **SVCLONE**: SVclone 1.1.4 (`quay.io/biocontainers/svclone:1.1.4--pyr44hdfd78af_0`), one task per sample.
  - Steps: annotate → count (DRAGEN tumour BAM) → filter → cluster → postassign.
  - Filter uses the SVs, the SASH PASS SNVs, PURPLE CN and purity/ploidy.
  - SNVs are subsampled to `snv_subsample` for clustering (seed 20261120).
  - postassign puts SVs on the SNV model.
  - Sex comes from PURPLE (`--XX`/`--XY`).
- **AMPTIMER**: AmplificationTimeR (Wedge-lab, pinned commit) on every focal amplification event in
  `events.tsv`.
  - CN, multiplicity (MutationTimeR MutCN) and SNVs are prepared upstream.
  - Clock-like CpG>TpG mutations are selected by AmplificationTimeR.
  - Runs on the Bioconductor 3.18 image; packages are installed at task start.

## Inputs
- `svclone_manifest` CSV columns:
  - `sample, sex, bam, bai, svs, cnv, purity_ploidy, snv_vcf` (S3 URIs), plus `n_svs`.
  - `svs`: SVclone "simple" format with `dir1/dir2/classification`.
  - `cnv`: Battenberg-style `chr startpos endpos nMaj1_A nMin1_A frac1_A ...`.
- `amptimer_inputs`: a directory with `events.tsv` and `<sample>.{cn,mult,muts}.tsv`.

## Run
```
nextflow run JIAJZHAN-Kairos/brca-sv-amp-timing-nf -r <commit> \
  --svclone_manifest s3://.../svclone/manifest.csv --amptimer_inputs s3://.../amptimer \
  --outdir s3://.../results [--samples L1,L2] [--run_svclone false] [--run_amptimer false]
```
