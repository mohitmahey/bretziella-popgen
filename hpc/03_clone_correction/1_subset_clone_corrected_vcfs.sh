#!/bin/bash
set -euo pipefail

# ── Setup ────────────────────────────────────────────────────────────────────
WORKDIR="/mnt/gs21/scratch/maheymoh/karan_project_jgi/AFTER_GATK/02_clone_correction"
SAMPLE_LIST="$WORKDIR/clone_corrected_samples_kept.txt"

VCF1="karan_biallelic_TAGs_snps_miss0.1_MAF.vcf.gz"
VCF2="karan_biallelic_TAGs_snps_miss0.1_MAF_LD.vcf.gz"

LOG="$WORKDIR/subset_run_$(date +%b-%d-%Y_%H-%M).log"

cd "$WORKDIR"

# Redirect all output to log and terminal
exec > >(tee -a "$LOG") 2>&1

echo "======================================================"
echo "Run date   : $(date)"
echo "Working dir: $WORKDIR"
echo "Sample list: $SAMPLE_LIST"
echo "Input VCF 1: $VCF1"
echo "Input VCF 2: $VCF2"
echo "======================================================"

module purge
module load BCFtools/1.22-GCC-13.3.0
export BCFTOOLS_PLUGINS=~/softwares/bcftools/plugins/

echo ""
echo "=== Tool versions ==="
bcftools --version | head -2
echo "======================================================"

# ── Subsetting VCF 1 ─────────────────────────────────────────────────────────
OUT1="${VCF1%.vcf.gz}_CC"
echo ""
echo "=== Subsetting $VCF1 ==="
bcftools view -S "$SAMPLE_LIST" -Oz -o "${OUT1}.vcf.gz" "$VCF1"
bcftools index "${OUT1}.vcf.gz"
bcftools view -Ov -o "${OUT1}.vcf" "${OUT1}.vcf.gz"

echo "SNPs in ${OUT1}.vcf.gz:"
bcftools view -H "${OUT1}.vcf.gz" | wc -l

# ── Subsetting VCF 2 ─────────────────────────────────────────────────────────
OUT2="${VCF2%.vcf.gz}_CC"
echo ""
echo "=== Subsetting $VCF2 ==="
bcftools view -S "$SAMPLE_LIST" -Oz -o "${OUT2}.vcf.gz" "$VCF2"
bcftools index "${OUT2}.vcf.gz"
bcftools view -Ov -o "${OUT2}.vcf" "${OUT2}.vcf.gz"

echo "SNPs in ${OUT2}.vcf.gz:"
bcftools view -H "${OUT2}.vcf.gz" | wc -l

echo ""
echo "=== DONE ==="
echo "Final files:"
echo "- ${OUT1}.vcf"
echo "- ${OUT2}.vcf"
echo "Finished: $(date)"
echo "Log saved to: $LOG"
