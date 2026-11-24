#!/bin/bash
# PhylogicNDT LeagueModel for one sample group, matching rerun_cohort_league.sh (Results/PhylogicNDT/
# seqera_brca196_20261118/scripts): same image, aggregate step (header-only comp.tsv files accepted, sample
# kept in the list so prevalence denominators are unchanged) and LeagueModel arguments/seed.
# Usage: run_league.sh <group> <samples.txt> <n_jobs>   (comp files staged as ./<s>.comp.tsv)
set -euo pipefail
G=$1; SAMPLES=$2; JOBS=$3
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg MPLCONFIGDIR=$PWD/.mpl
C=$PWD/c; mkdir -p "$C/comparisons" "$C/$G"
while read -r s; do [ -n "$s" ] || continue; mkdir -p "$C/comparisons/$s"; cp "$s.comp.tsv" "$C/comparisons/$s/$s.comp.tsv"; done < "$SAMPLES"
cp "$SAMPLES" "$C/$G/samples.txt"
cp -r /opt/phylogic/workflow wf
python - <<'PY'
p='wf/prepare_comparisons.py'; s=open(p).read()
old="            if count == 0:\n                raise ValueError('Header-only comparison file: '+str(source))\n"
assert s.count(old)==1
s=s.replace(old, "            # rerun 20261119: 0-row files allowed (no comparable events); sample kept in list\n")
open(p,'w').write(s)
PY
( cd wf && python -c "from pathlib import Path; from prepare_comparisons import aggregate; aggregate(Path('$C/$G/samples.txt'), Path('$C/comparisons'), Path('$C/$G/comparisons.tsv'))" )
cd "$C/$G"
python /opt/phylogic/PhylogicNDT.py LeagueModel --cohort "$G" \
  --comparison_fn "$C/$G/comparisons.tsv" --final_sample_list "$C/$G/samples.txt" \
  --n_perms 500 --n_seasons 200 --n_jobs "$JOBS" --perm_chunk_size 10 \
  --percent_subset 0.8 --min_event_prevalence 0.05 --seed 20260813 > league.log 2>&1 || { tail -30 league.log; exit 1; }
echo "LeagueModel done: $G"
