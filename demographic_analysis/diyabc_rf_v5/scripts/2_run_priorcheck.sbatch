#!/bin/bash --login
# Prior-scenario check on the finished reftable (manual 4.1, 5.4.7; same command as the GUI).
# Output (moved to results/prior_check/): pcaloc1_locate.txt (stars = mismatch), pcaloc1_ACP.txt (PCA coordinates).
# Usage (from diyabc_rf/, after 1_run_reftable.sbatch has finished):  sbatch scripts/2_run_priorcheck.sbatch
#SBATCH --job-name=diyabc_priorcheck
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:30:00
#SBATCH --output=/mnt/ffs24/home/maheymoh/karan_bretfa_GBS/diyabc_rf_v5/logs/priorcheck_%j.log

set -euo pipefail
cd /mnt/ffs24/home/maheymoh/karan_bretfa_GBS/diyabc_rf_v5   # project folder, fixed: works from any submit directory

./bin/diyabc-RF-linux-v1.1.54 -p ./ -d a:pl

mv -f pcaloc1_* results/prior_check/
