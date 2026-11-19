#!/bin/bash
set -euo pipefail

# ── Setup ────────────────────────────────────────────────────────────────────
WORKDIR="/mnt/gs21/scratch/maheymoh/karan_project_jgi/AFTER_GATK/01_filter_snps_maf_missing"
INPUT_VCF="/mnt/gs21/scratch/maheymoh/karan_project_jgi/AFTER_GATK/00_input_files/snps_passed_filtered.vcf"
LOG="$WORKDIR/filter_run_$(date +%b-%d-%Y_%H-%M).log"

cd "$WORKDIR"

# Redirect all output to log and terminal
exec > >(tee -a "$LOG") 2>&1

echo "======================================================"
echo "Run date   : $(date)"
echo "Working dir: $WORKDIR"
echo "Input VCF  : $INPUT_VCF"
echo "======================================================"

module purge
module load BCFtools/1.22-GCC-13.3.0
export BCFTOOLS_PLUGINS=~/softwares/bcftools/plugins/

echo ""
echo "=== Tool versions ==="
bcftools --version | head -2
echo "======================================================"

# ── Starting SNP count ───────────────────────────────────────────────────────
echo ""
echo "=== Starting SNP count ==="
bcftools view -H "$INPUT_VCF" | wc -l

# ── Biallelic filter ─────────────────────────────────────────────────────────
echo "=== Keeping biallelic sites only ==="
bcftools view -m2 -M2 -Oz \
    "$INPUT_VCF" \
    -o karan_biallelic_snps.vcf.gz

echo "SNPs after biallelic filter:"
bcftools view -H karan_biallelic_snps.vcf.gz | wc -l

# ── Fill tags ────────────────────────────────────────────────────────────────
echo "=== Filling tags ==="
bcftools +fill-tags karan_biallelic_snps.vcf.gz -Oz -o karan_biallelic_TAGs_snps.vcf.gz -- -t all
bcftools index karan_biallelic_TAGs_snps.vcf.gz

# ── Missing data filter (F_MISSING ≤ 0.1) ───────────────────────────────────
echo "=== Filtering by missingness (F_MISSING > 0.1 removed) ==="
bcftools filter \
    -e 'F_MISSING > 0.1' \
    karan_biallelic_TAGs_snps.vcf.gz \
    -Oz -o karan_biallelic_TAGs_snps_miss0.1.vcf.gz
bcftools index karan_biallelic_TAGs_snps_miss0.1.vcf.gz

echo "SNPs after missingness filter:"
bcftools view -H karan_biallelic_TAGs_snps_miss0.1.vcf.gz | wc -l

# ── MAF filter (MAF ≥ 0.05) ──────────────────────────────────────────────────
echo "=== Filtering by MAF (removing MAF < 0.05) ==="
bcftools view karan_biallelic_TAGs_snps_miss0.1.vcf.gz -Oz -e 'MAF<0.05' -o karan_biallelic_TAGs_snps_miss0.1_MAF.vcf.gz
bcftools index karan_biallelic_TAGs_snps_miss0.1_MAF.vcf.gz

echo "SNPs after MAF filter:"
bcftools view -H karan_biallelic_TAGs_snps_miss0.1_MAF.vcf.gz | wc -l


# ── LD decay (PopLDdecay) ────────────────────────────────────────────────────
#echo "=== Calculating LD decay with PopLDdecay ==="
# Running this BEFORE LD pruning, but AFTER MAF and Missing data filtering to prevent crashes!
#conda run -n popLDdecay \
#    PopLDdecay \
#        -InVCF  karan_biallelic_TAGs_snps_MAF_miss0.1.vcf.gz \
#        -OutStat karan_after_gatk_snps_only

# NOTE: To generate the LD decay plot on a Mac, run locally:
#   perl Plot_OnePop.pl -inFile karan_after_gatk_snps_only.stat.gz -out OUT

# ── LD pruning ───────────────────────────────────────────────────────────────
echo "=== LD pruning (r2 ≤ 0.2, window 50 kb) ==="
bcftools +prune \
    -Oz \
    -o karan_biallelic_TAGs_snps_miss0.1_MAF_LD.vcf.gz \
    -m 0.2 \
    --window 50000bp \
    karan_biallelic_TAGs_snps_miss0.1_MAF.vcf.gz
bcftools index karan_biallelic_TAGs_snps_miss0.1_MAF_LD.vcf.gz

echo "=== Final SNP count (after LD Pruning) ==="
bcftools view -H karan_biallelic_TAGs_snps_miss0.1_MAF_LD.vcf.gz | wc -l

bcftools view -Ov -o karan_biallelic_TAGs_snps_miss0.1_MAF_LD.vcf karan_biallelic_TAGs_snps_miss0.1_MAF_LD.vcf.gz

echo ""
echo "=== DONE ==="
echo "Final file for downstream analysis: karan_biallelic_TAGs_snps_miss0.1_MAF_LD.vcf"
echo "Finished: $(date)"
echo "Log saved to: $LOG"
