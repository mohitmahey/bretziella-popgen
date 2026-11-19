# bretziella-popgen

Scripts used for analysis of population genetic study of the oak wilt fungus *Bretziella fagacearum* in the United States. 

## Workflow

```
 raw GBS reads (SRA)
        |
        v
 hpc/01_variant_calling      trim, align, call variants, quality report
        |
        v
 hpc/02_snp_filtering        filter SNPs (MAF, missing data), prune LD
        |
        v
 hpc/03_clone_correction     find clones, keep one isolate per clone
        |
        +-----------------------+-------------------------+-----------------------+
        |                       |                         |                       |
        v                       v                         v                       v
 hpc/04_structure         rScripts/                  ML tree                 demographic_analysis/
 PLINK + fastSTRUCTURE    PCA, DAPC, diversity,      vcf2phylip +            diyabc_rf_v5/scripts
        |                 differentiation, AMOVA,    IQ-TREE                 DIYABC-RF scenarios
        v                 PERMANOVA
 rScripts/
 ancestry barplots and maps
```

## Scripts

`hpc/` holds the steps that ran on a SLURM cluster. The scripts contain absolute paths of the original cluster. Change these paths before you run them.

- `01_variant_calling/`: trims the reads (fastp), aligns them to the reference genome (BWA-MEM, SAMtools, Picard), calls variants (GATK), and makes a quality report (MultiQC). Run the scripts in number order.
- `02_snp_filtering/`: removes SNPs by minor allele frequency and missing data, and prunes SNPs in linkage disequilibrium (LD).
- `03_clone_correction/`: finds clonal isolates with poppr (`clone_correction.Rmd`), and keeps one isolate for each clone in the VCF files (`1_subset_clone_corrected_vcfs.sh`).
- `04_structure/`: converts the VCF to PLINK format, adds population labels to the `.fam` file, and runs fastSTRUCTURE for a range of K.

`rScripts/` holds the R analyses. Run them from the repository root, for example `Rscript rScripts/stats.R`. Each script makes `figures/` and `results/` if they do not exist.

- `color_palette.R`: the colours and shapes of the populations. The other scripts read it with `source()`.
- `map_oak_wilt_spread.R`: a US county map of oak wilt records with the sample locations.
- `GBS_analysis_postVCF_v2_clone_corrected.R`: PCA, DAPC, and the K-means cluster search.
- `stats.R`: genetic diversity, private alleles, pairwise FST, Nei's distance, Jost's D, and AMOVA.
- `permanova_host.R`: PERMANOVA of population, collection year, and host.
- `structure_barplot_pophelper.R`: fastSTRUCTURE ancestry barplots.
- `maps.R`: fastSTRUCTURE ancestry pies for each isolate on a US map.
- `maps_flow_structure.R`: mean fastSTRUCTURE ancestry pies for each population on a US map.

`demographic_analysis/diyabc_rf_v5/scripts/` holds the DIYABC-RF analysis. Run the scripts in number order. They convert the VCF to the DIYABC format, make the reference table, do a prior check, select the best scenario, and estimate its parameters. `4_make_check_figs.R` draws the check figures. Download DIYABC-RF and abcranger from https://diyabc.github.io/.


## Maximum-likelihood tree

Convert the clone-corrected VCF to PHYLIP format. Then run IQ-TREE with ModelFinder and ultrafast bootstrap:

```
python vcf2phylip.py -i <clone_corrected>.vcf
iqtree3 -s <clone_corrected>.min4.phy -bb 10000 -bnni -nt AUTO
```

## Licence

The scripts are under the MIT licence (`LICENSE`). `demographic_analysis/diyabc_rf_v5/scripts/1_vcf2diyabc.py` is by Etienne Loire and is under GPL-3.0 (`COPYING` in the same folder).
