#!/usr/bin/env bash
set -euo pipefail

# ── Setup ────────────────────────────────────────────────────────────────────
WORKDIR="/mnt/gs21/scratch/maheymoh/karan_project_jgi/AFTER_GATK/04_structure/fastStructure"
INPUT_PREFIX="/mnt/gs21/scratch/maheymoh/karan_project_jgi/AFTER_GATK/03_plink/karan_snps_filtered_plink"
OUTPUT_DIR="$WORKDIR/output"
OUTPUT_PREFIX="$OUTPUT_DIR/karan"
LOG="$WORKDIR/faststructure_run_$(date +%b-%d-%Y_%H-%M).log"

cd "$WORKDIR"

# Redirect all output to log and terminal
exec > >(tee -a "$LOG") 2>&1

echo "======================================================"
echo "Run date   : $(date)"
echo "Working dir: $WORKDIR"
echo "Input prefix: $INPUT_PREFIX"
echo "Output prefix: $OUTPUT_PREFIX"
echo "======================================================"

# Load environment
# Assuming conda is available in the shell
source $(conda info --base)/etc/profile.d/conda.sh || { echo "Error: Could not find conda.sh"; exit 1; }
conda activate faststructure || { echo "Error: Could not activate faststructure environment"; exit 1; }

echo ""
echo "=== Tool versions ==="
# fastStructure doesn't have a standard --version flag that works well in all installs, 
# but we can try to find where it is
which fastStructure
which chooseK.py

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# ── Run fastStructure ───────────────────────────────────────────────────────
echo ""
echo "=== Starting fastStructure runs for K=2 to 20 ==="

for k in $(seq 2 20); do
    echo "[$(date +'%H:%M:%S')] Running for K=${k}..."
    fastStructure --input="${INPUT_PREFIX}" -K ${k} --output="${OUTPUT_PREFIX}" --full --seed=1
done

# ── Choose K ────────────────────────────────────────────────────────────────
echo ""
echo "=== Determining the best K using chooseK.py ==="
chooseK.py --input="${OUTPUT_PREFIX}"

echo ""
echo "=== DONE ==="
echo "Finished: $(date)"
echo "Log saved to: $LOG"
