#!/bin/bash --login

# If popfile was saved from Excel/Windows, fix it with:
# tr -d '\r' < popfile_original.tsv > tmp && mv tmp popfile_original.tsv
# [ -n "$(tail -c1 popfile_original.tsv)" ] && echo >> popfile_original.tsv

set -e
cd ../input
exec > >(tee 0_build_input.log) 2>&1     # everything below goes to the screen AND input/0_build_input.log

VCF=karan_biallelic_TAGs_snps_miss0.1_MAF_LD_CC.vcf

# 1. Convert VCF to DIYABC format (answers: VCF name, popfile, enter)
#    -> karan_biallelic_TAGs_snps_miss0.1_MAF_LD_CC.DIYABC.snp (all loci coded A = diploid)
module purge
module load Python/2.7.18-GCCcore-12.3.0
printf '%s\npopfile_original.tsv\n\n' $VCF | python2 ../scripts/1_vcf2diyabc.py

# 2. Header line 1 + loci type A -> H (haploid)
sed -e '1s/.*/<NM=1NF> <MAF=0.05> Bretziella fagacearum GBS haploid/' \
    -e '2s/\tA/\tH/g' \
    karan_biallelic_TAGs_snps_miss0.1_MAF_LD_CC.DIYABC.snp > karan_LD_CC_haploid.DIYABC.snp

# 3. Done: summary and tool versions
echo
echo "Done: wrote input/karan_LD_CC_haploid.DIYABC.snp"
echo "Tools used:"
echo "  $(python2 --version 2>&1)  (Python/2.7.18-GCCcore-12.3.0 module)"
echo "  1_vcf2diyabc.py  (converter script, no version number)"
echo "  $(sed --version | head -1)"
echo "  bash $BASH_VERSION"
