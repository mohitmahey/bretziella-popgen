#!/bin/bash --login
# RF scenario choice on all scenarios (manual 6.2), using all simulations.
# Scenario numbers S1, S2, ... match the order in header.txt.
# Output in results/model_choice/: modelchoice_all.predictions (votes, posterior prob.), .confusion,
#   .importance, .ooberror (error vs trees, manual 4.3), .lda
# Usage (from diyabc_rf/):  sbatch scripts/3_run_modelchoice.sbatch
#SBATCH --job-name=abcrf_modelchoice
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=2:30:00
#SBATCH --output=/mnt/ffs24/home/maheymoh/karan_bretfa_GBS/diyabc_rf_v5/logs/modelchoice_%j.log

set -euo pipefail
cd /mnt/ffs24/home/maheymoh/karan_bretfa_GBS/diyabc_rf_v5   # project folder, fixed: works from any submit directory

./bin/abcranger-linux-v1.16.69 -h header.txt -r reftable.bin -b statobs.txt -t 2000 -j "$SLURM_CPUS_PER_TASK" -s 1 --allpost -o modelchoice_all

mv -f modelchoice_all* results/model_choice/
