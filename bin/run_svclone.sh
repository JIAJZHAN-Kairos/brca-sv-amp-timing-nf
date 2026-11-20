#!/usr/bin/env bash
# SVclone 1.1.4 for one sample: annotate -> count -> filter (SVs + SNVs, CNV, purity) -> cluster -> postassign
# (SVs assigned to the SNV model). Each step's exit status is recorded; later steps are skipped after a failure.
# Usage: run_svclone.sh <sample> <sex> <bam> <svs> <cnv> <purity_ploidy> <snv_vcf> <cpus> <snv_subsample>
set -uo pipefail
s=$1; sex=$2; bam=$3; svs=$4; cnv=$5; pp=$6; snv=$7; cpus=$8; nsub=$9
chroms=$(printf 'chr%s,' $(seq 1 22) X | sed 's/,$//')
cat > svclone.ini <<INI
[BamParameters]
read_len: -1
insert_mean: -1
insert_std: -1
mean_cov: 50
max_cn: 10
[SVannotateParameters]
use_dir: True
sv_class_field: classification
trust_sc_position: False
[SVcountParameters]
threshold: 6
norm_overlap: 10
sc_len: 10
[SVclasses]
inversion_class: INV
deletion_class: DEL
dna_gain_class: DUP,INTDUP
dna_loss_class: DEL,INV,TRX
itrx_class: INTRX
[SocratesOpts]
pos1: C1_anchor
dir1: C1_anchor_dir
pos2: C1_realign
dir2: C1_realign_dir
avg_mapq1: C1_avg_realign_mapq
avg_mapq2: C2_avg_realign_mapq
repeat1: repeat1
repeat2: repeat2
filter_repeats: Satellite,Simple_repeat
min_mapq: 20
[DebugParameters]
write_anomalous: False
[FilterParameters]
neutral: False
size_filter: -1
min_dep: 8
min_span: 1
min_split: 1
filter_outliers: False
filter_chroms: True
filter_subclonal_cnvs: False
strict_cnv_filt: True
germline_threshold: 10
support_adjust_factor: 0
sv_offset: 100000
[ValidationParameters]
chroms: ${chroms}
[ClusterParameters]
threads: ${cpus}
repeat: 5
n_iter: 1000
clus_limit: 6
adjusted: True
male: $([ "$sex" = male ] && echo True || echo False)
sv_to_sim: 0
INI
sexflag=$([ "$sex" = male ] && echo --XY || echo --XX)
status=status.tsv; echo -e "sample\tstep\texit" > $status
run() { local step=$1; shift; echo "== $step: $*"; "$@"; local rc=$?; echo -e "$s\t$step\t$rc" >> $status; return $rc; }
run annotate svclone annotate -i "$svs" -b "$bam" -s "$s" --sv_format simple -cfg svclone.ini &&
run count svclone count -i "$s/${s}_svin.txt" -b "$bam" -s "$s" -cfg svclone.ini &&
run filter svclone filter -s "$s" -i "$s/${s}_svinfo.txt" -p "$pp" -c "$cnv" --snvs "$snv" --snv_format mutect -cfg svclone.ini &&
run cluster svclone cluster -s "$s" -cfg svclone.ini $sexflag --subsample "$nsub" --ss_seed 20261120 &&
run postassign svclone postassign -s "$s"
true
