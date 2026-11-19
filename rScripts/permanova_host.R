# PERMANOVA: does host sub-section explain genetic variation in Bretziella
# fagacearum after population (Clusters) and collection year?
# Standalone; same load idiom and distance as the PERMANOVA in stats.R.
# Run from the repo root: Rscript rScripts/permanova_host.R
rm(list = ls())

library(poppr)       # diss.dist
library(vegan)       # PERMANOVA (adonis2)
library(dplyr)
library(vcfR)

setwd(".")
outdir <- "results"
dir.create(outdir, showWarnings = FALSE)

# Host groups to drop. Empty = keep all four (Red oak, White oak, Chestnut, UA).
# e.g. c("UA", "Chestnut") tests red oak vs white oak only.
drop_hosts <- character(0)

# 1. Load data ----
bret_vcf <- read.vcfR("vcf/karan_biallelic_TAGs_snps_miss0.1_MAF_CC.vcf")

pop_data <- read.csv("metadata/clone_corrected_metadata.csv", header = TRUE,
                     stringsAsFactors = FALSE, fileEncoding = "UTF-8")
pop_data <- pop_data %>% mutate(across(where(is.character), trimws))

# ensuring metadata matches VCF samples
vcf_samples <- colnames(bret_vcf@gt)[-1]
pop_data <- pop_data[match(vcf_samples, pop_data$Sakalidis_lab_code), ]
stopifnot(identical(pop_data$Sakalidis_lab_code, vcf_samples))

gi_bret <- vcfR2genind(bret_vcf, return.alleles = TRUE)
ploidy(gi_bret) <- 1
pop(gi_bret) <- pop_data$Clusters

# 2. Keep isolates with a known year and a kept host group ----
pop_data$Year <- as.numeric(pop_data$Collection.year)
keep     <- !is.na(pop_data$Year) & !pop_data$Sub.section %in% drop_hosts
meta_sub <- pop_data[keep, ]
meta_sub$Host <- factor(meta_sub$Sub.section)

cat("Isolates kept:", sum(keep), "of", nrow(pop_data), "\n")
print(table(meta_sub$Host, meta_sub$Clusters))

# 3. Pairwise genetic distance (same as stats.R) ----
gen_dist    <- diss.dist(gi_bret, percent = FALSE)
dist_matrix <- as.dist(as.matrix(gen_dist)[keep, keep])

# 4. Three-factor model, marginal (type III) effects ----
# Each term is tested after the other two, so SS and R2 do not sum to the total.
set.seed(1)
perm_margin <- adonis2(dist_matrix ~ Clusters + Host + Year, data = meta_sub,
                       permutations = 999, by = "margin")

# 5. Host tested within populations (permutations restricted to Clusters) ----
set.seed(1)
perm_within <- adonis2(dist_matrix ~ Host, data = meta_sub,
                       strata = meta_sub$Clusters, permutations = 999)

print(perm_margin)
print(perm_within)

# 6. Save outputs ----
sink(file.path(outdir, "permanova_clusters_host_year_margin.txt"))
cat("=== PERMANOVA: Clusters + Host sub-section + Year (Marginal Effects) ===\n")
cat("Dropped host groups:", paste(drop_hosts, collapse = ", "), "\n")
cat("Isolates:", sum(keep), "\n\n")
print(perm_margin)
sink()

sink(file.path(outdir, "permanova_host_within_clusters.txt"))
cat("=== PERMANOVA: Host sub-section Tested Within Clusters (Strata) ===\n")
print(perm_within)
sink()

write.csv(as.data.frame(perm_margin),
          file = file.path(outdir, "permanova_clusters_host_year_margin.csv"),
          quote = FALSE)
write.csv(as.data.frame(perm_within),
          file = file.path(outdir, "permanova_host_within_clusters.csv"),
          quote = FALSE)
