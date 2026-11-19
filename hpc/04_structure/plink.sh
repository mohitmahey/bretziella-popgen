#!/bin/bash
set -euo pipefail

# ── Setup ────────────────────────────────────────────────────────────────────
WORKDIR="/mnt/gs21/scratch/maheymoh/karan_project_jgi/AFTER_GATK/03_plink"
INPUT_VCF="../02_clone_correction/karan_biallelic_TAGs_snps_miss0.1_MAF_LD_CC.vcf"
OUT_PREFIX="karan_snps_filtered_plink"
LOG="$WORKDIR/plink_run_$(date +%b-%d-%Y_%H-%M).log"

cd "$WORKDIR"

# Redirect all output to log and terminal
exec > >(tee -a "$LOG") 2>&1

echo "======================================================"
echo "Run date   : $(date)"
echo "Working dir: $WORKDIR"
echo "Input VCF  : $INPUT_VCF"
echo "======================================================"

# Load PLINK module
module purge
module load PLINK/2.00a3.7-gfbf-2023a

echo ""
echo "=== PLINK version ==="
plink --version
echo "======================================================"

# ── 1. Convert VCF to Binary (BED/BIM/FAM) ────────────────────────────────────
echo "=== Converting VCF to PLINK binary format ==="
plink --vcf "$INPUT_VCF" \
      --make-bed \
      --out "$OUT_PREFIX" \
      --allow-extra-chr \
      --double-id

# ── 2. Update BIM IDs ────────────────────────────────────────────────────────
# Replaces '.' IDs with 'CHR:POS' format
echo "=== Updating BIM file IDs (CHR:POS) ==="
awk '{if($2==".") $2=$1":"$4; print}' "${OUT_PREFIX}.bim" > temp_binary_with_ids.bim
mv temp_binary_with_ids.bim "${OUT_PREFIX}.bim"

# ── 3. Handle FAM file (Optional Update) ─────────────────────────────────────
# NOTE: The original script expected 'karan_snps_updated_fam.txt' which was not found.
# If you have an updated FAM file with population data, uncomment and adjust the line below.
# mv karan_snps_updated_fam.txt "${OUT_PREFIX}.fam"

# ── 4. Recode to Structure format ───────────────────────────────────────────
echo "=== Recoding to Structure format ==="
plink --bfile "$OUT_PREFIX" \
      --recode structure \
      --out "${OUT_PREFIX}_structure" \
      --allow-extra-chr

# ── 5. Calculate Missingness ────────────────────────────────────────────────
echo "=== Calculating missingness statistics ==="
plink --bfile "$OUT_PREFIX" \
      --double-id \
      --allow-extra-chr \
      --missing \
      --out miss

echo ""
echo "=== DONE ==="
echo "Finished: $(date)"
echo "Log saved to: $LOG"
