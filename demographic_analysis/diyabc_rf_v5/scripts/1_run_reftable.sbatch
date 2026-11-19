#!/bin/bash --login
#SBATCH --job-name=diyabc_reftable
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/mnt/ffs24/home/maheymoh/karan_bretfa_GBS/diyabc_rf_v5/logs/reftable_%j.log

set -euo pipefail
cd /mnt/ffs24/home/maheymoh/karan_bretfa_GBS/diyabc_rf_v5   # project folder, fixed: works from any submit directory

BIN=./bin/diyabc-RF-linux-v1.1.54
NSIM=120000
T=$SLURM_CPUS_PER_TASK

# RNG file sets the maximum thread count; make it once (manual section 9)
[ -f RNG_state_0000.bin ] || "$BIN" -p ./ -n "t:${T}"

# -g = sims per batch kept in RAM before writing; larger batches keep all cores busy
"$BIN" -p ./ -m -t "$T" -g 500 -r "$NSIM"
