#!/bin/bash --login
# RF parameter estimation under the winning scenario K (manual 6.3): one forest per parameter.
# Output per parameter (moved to results/estimation/): estim_S<K>_<name>.predictions (mean, median, 5%-95% quantiles),
#   estim_S<K>_<name>.oobstats (NMAE, coverage), estim_S<K>_<name>.importance, ...
# PLS axes on with 95% explained variance (manual default; binary default is 0.9).
# --noob must be >= 1000 and <= number of sims of scenario K (~20,000 per scenario).
#
# Usage (from diyabc_rf/; K and PARAMS are required; parameter names come from header.txt):
#   sbatch --export=ALL,K=1,PARAMS="t1 t2 r1" scripts/5_run_estim.sbatch
# Check that every name in PARAMS exists in scenario K (model/priors.txt, "used in scenarios").
#SBATCH --job-name=abcrf_estim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/mnt/ffs24/home/maheymoh/karan_bretfa_GBS/diyabc_rf_v5/logs/estim_%j.log

set -euo pipefail
cd /mnt/ffs24/home/maheymoh/karan_bretfa_GBS/diyabc_rf_v5   # project folder, fixed: works from any submit directory

BIN=./bin/abcranger-linux-v1.16.69
J=$SLURM_CPUS_PER_TASK
K=${K:?set the winning scenario K}
PARAMS=${PARAMS:?set PARAMS explicitly; every name must exist in scenario K (model/priors.txt)}

for p in $PARAMS; do
  name=$(echo "$p" | sed 's#/#_over_#g')     # tE/N1 -> tE_over_N1 (file-name safe)
  "$BIN" -h header.txt -r reftable.bin -b statobs.txt -t 1000 -j $J -s 1 --chosenscen "$K" --noob 2000 --plsmaxvar 0.95 \
         --parameter "$p" -o "estim_S${K}_${name}"
done

mv -f "estim_S${K}_"* results/estimation/
