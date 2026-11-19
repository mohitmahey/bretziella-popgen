rm(list=ls())
set.seed(1)

#############################################
## 0. Packages & data
#############################################

# ==============================================================================
# Package Summary Table:
# Package   | Key Functions Used                          | What It Calculates / Does
# ----------|---------------------------------------------|-------------------------------------------------------------
# vcfR      | read.vcfR(), vcfR2genlight(), vcfR2genind() | Parses VCF files; converts genomic variant calls into
#           |                                             | population genetic objects.
# adegenet  | genind, ploidy(), pop(), as.matrix()        | Stores haploid genotypes; sets population groupings;
#           |                                             | generates allele dosage matrices.
# hierfstat | allelic.richness(), basic.stats(),          | Calculates Allelic Richness (A_R), Gene Diversity (Hs),
#           | pairwise.WCfst(), genet.dist()              | Pairwise F_ST (Weir & Cockerham 1984), and Nei's genetic
#           |                                             | distance (Ds).
# poppr     | private_alleles(), strata(), poppr.amova()  | Counts population-specific private alleles; sets geographic
#           |                                             | strata; computes AMOVA (SS, MS, sigma^2, Phi_PT).
# ade4      | randtest()                                  | Runs the 999-permutation Monte Carlo significance test for
#           |                                             | AMOVA (Z-score, simulated P-value).
# mmod      | pairwise_D()                                | Calculates Pairwise Jost's D (D_est) allelic
#           |                                             | differentiation between all populations.
# vegan     | adonis2()                                   | Runs PERMANOVA to test the marginal effects of Geography
#           |                                             | (Clusters) vs. Collection Year (Year).
# dplyr     | mutate(), across(), filter(), arrange()     | Data wrangling; cleans whitespace; matches sample IDs
#           |                                             | between metadata and VCF.
# ggplot2   | ggplot(), geom_tile(), geom_col(), ggsave() | Renders heatmaps (F_ST, Jost's D, Nei's Ds) and diversity
#           |                                             | barplots (He, private alleles).
# scales    | rescale(), label_percent()                  | Adjusts color palette gradients and axis label formatting on
#           |                                             | heatmaps and barplots.
# ==============================================================================

library(adegenet)    # genind, PCA
library(hierfstat)   # Na, AR, F-statistics
library(poppr)       # poppr stats, AMOVA wrapper
library(mmod)        # Gst, pairwise Gst
library(ade4)        # randtest for AMOVA
# library(vegan)       # PERMANOVA (adonis2)
library(dplyr)
library(ggplot2)
library(scales)
library(vcfR)
# setting up working directory
setwd(".")

# loading vcf file into the env
bret_vcf <- read.vcfR("vcf/karan_biallelic_TAGs_snps_miss0.1_MAF_CC.vcf")

# loading the population location data file
pop_data <- read.csv("metadata/clone_corrected_metadata.csv", header = TRUE,
                     stringsAsFactors = FALSE, fileEncoding = "UTF-8")
# remove leading/trailing whitespace so population keys match cleanly
pop_data <- pop_data %>% mutate(across(where(is.character), trimws))

# ensuring metadata matches VCF samples
vcf_samples <- colnames(bret_vcf@gt)[-1]
pop_data <- pop_data[match(vcf_samples, pop_data$Sakalidis_lab_code), ]

# converting it to genlight to be able to be used by poppr
gl_bret <- vcfR2genlight(bret_vcf)
gi_bret <- vcfR2genind(bret_vcf,return.alleles = TRUE)

# setting the ploidy level of the species (Bretziella fagacearum is haploid)
ploidy(gl_bret) <- 1
ploidy(gi_bret) <- 1

# population grouping = Clusters (9 US groups; Eastern merges MD/NY/OH/PA/SC/WV)
pop(gl_bret) <- pop_data$Clusters
pop(gi_bret) <- pop_data$Clusters

# checking the data
gl_bret
gi_bret

# dosage matrix
do_bret <- as.matrix(gi_bret)
do_bret[1:5,1:5]
str(do_bret)

# number of populations
nPop(gl_bret)

# name of pop
popNames(gl_bret)

# Output directories
outdir <- "results"
figdir <- "figures"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
if (!dir.exists(figdir)) dir.create(figdir, recursive = TRUE)

# Source palette, to keep colors consistent acorss the board
source("rScripts/color_palette.R")
cols <- get_pop_colors(popNames(gl_bret))

pop_labels <- get_pop_labels(popNames(gl_bret))


#############################################
## 1. Convert to hierfstat format
#############################################

hf_dat <- genind2hierfstat(gi_bret)

pop_names <- levels(hf_dat$pop)
N_pop     <- table(pop(gi_bret))

#############################################
## 2. Na and allelic richness (AR) per population
#############################################

# Na: mean number of alleles per locus, per population (not rarefied -> all pops)
na_mat <- nb.alleles(hf_dat, diploid = FALSE)
colnames(na_mat) <- pop_names
Na_pop <- setNames(colMeans(na_mat, na.rm = TRUE), pop_names)

# AR: allelic richness rarefied to the smallest well-sampled population. A single
# (or near-single) haploid genome can't support rarefied richness, so pops with
# < 3 isolates (AR N=1, IL N=2) are dropped from AR and the rest are rarefied to
# the smallest observed sample among them. Dropped pops appear as NA.
keep      <- names(N_pop)[N_pop >= 3]
hf_ar     <- hf_dat[as.character(hf_dat$pop) %in% keep, , drop = FALSE]
hf_ar$pop <- droplevels(hf_ar$pop)
ar_res    <- allelic.richness(hf_ar,min.n = 3, diploid = FALSE)
AR_mat    <- ar_res$Ar
colnames(AR_mat) <- levels(hf_ar$pop)
AR_pop    <- setNames(rep(NA_real_, length(pop_names)), pop_names)
AR_pop[colnames(AR_mat)] <- colMeans(AR_mat, na.rm = TRUE)

# Save raw Na & AR matrices (optional)
write.csv(na_mat,
          file = file.path(outdir, "Na_per_locus_per_pop.csv"),
          quote = FALSE)
write.csv(AR_mat,
          file = file.path(outdir, "AR_per_locus_per_pop.csv"),
          quote = FALSE)

na_ar_summary <- data.frame(
  Pop = pop_names,
  Na_mean = as.numeric(Na_pop[pop_names]),
  AR_mean = as.numeric(AR_pop[pop_names])
)

write.csv(
  na_ar_summary,
  file = file.path(outdir, "Na_AR_population_means.csv"),
  row.names = FALSE,
  quote = FALSE
)

#############################################
## 2b. Private alleles per population (poppr)
#############################################

# Tabulate private alleles (count.alleles = FALSE counts each unique private allele once)
pa_mat <- private_alleles(gi_bret, report = "table", count.alleles = FALSE)

# Ensure all 9 populations are accounted for in standard order
pa_pop <- setNames(rep(0, length(pop_names)), pop_names)
raw_pa <- rowSums(pa_mat > 0)
pa_pop[names(raw_pa)] <- raw_pa

# Export dedicated private alleles summary CSV
write.csv(
  data.frame(
    Pop = pop_names,
    N   = as.numeric(N_pop[pop_names]),
    Private_Alleles = as.numeric(pa_pop[pop_names])
  ),
  file = file.path(outdir, "private_alleles_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

#############################################
## 3. F-statistics with hierfstat
#############################################

bs <- basic.stats(hf_dat, diploid = FALSE)

sink(file.path(outdir, "basic_stats.txt"))
print(bs$overall)
sink()

# Ho_locus  <- bs$Ho     # haploid: observed heterozygosity undefined (NaN)
He_locus  <- bs$Hs
# Fis_locus <- bs$Fis    # haploid: Fis undefined (derived from Ho)
Fst_locus <- bs$Fst

# Per-population means
# Ho_pop  <- colMeans(Ho_locus,  na.rm = TRUE)   # haploid: Ho undefined
He_pop  <- colMeans(He_locus,  na.rm = TRUE)
# Fis_pop <- colMeans(Fis_locus, na.rm = TRUE)   # haploid: Fis undefined

# Global F-statistics from basic.stats()$overall (Weir & Cockerham).
# Haploid: only FST is defined; FIT and FIS are heterozygosity-based (NaN).
global_Fstat <- data.frame(
  FST = bs$overall["Fst"]
  # FIT = bs$overall["Fit"],   # haploid: FIT undefined
  # FIS = bs$overall["Fis"]    # haploid: FIS undefined
)

write.csv(global_Fstat,
          file = file.path(outdir, "global_Fstats_wc.csv"),
          row.names = FALSE,
          quote = FALSE)


#############################################
## 4. Per-population summary table
#############################################

div_table <- data.frame(
  Pop = pop_names,
  N   = as.numeric(N_pop[pop_names]),
  Na  = Na_pop[pop_names],
  AR  = AR_pop[pop_names],
  Pa  = as.numeric(pa_pop[pop_names]),
  He  = He_pop[pop_names],
  # Ho  = Ho_pop[pop_names],    # haploid: Ho undefined (no heterozygotes)
  # Fis = Fis_pop[pop_names],   # haploid: Fis undefined (derived from Ho)
  row.names = NULL
)

write.csv(div_table,
          file = file.path(outdir, "population_diversity_summary.csv"),
          row.names = FALSE,
          quote = FALSE)

#############################################
## 5. Gene diversity (He) per population — haploid plot
#############################################
# Haploid: observed heterozygosity (Ho) is undefined, so plot gene diversity
# (He = Nei's Hs) alone. Old Ho-vs-He barplot kept below, commented out.

he_df <- data.frame(Pop = pop_names, He = He_pop[pop_names])

p_he <- ggplot(he_df, aes(x = Pop, y = He)) +
  geom_bar(stat = "identity", fill = "#d6604d", width = 0.6) +
  scale_x_discrete(labels = pop_labels) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Population", y = "Gene diversity (He)")

ggsave(
  filename = file.path(figdir, "gene_diversity_He_barplot.png"),
  plot     = p_he,
  width    = 8,
  height   = 5,
  dpi      = 300
)

#############################################
## 5b. Private alleles per population — barplot
#############################################

pa_df <- data.frame(
  Pop = pop_names,
  Private_Alleles = as.numeric(pa_pop[pop_names])
)
pa_df <- pa_df[order(-pa_df$Private_Alleles), ]
pa_df$Pop <- factor(pa_df$Pop, levels = pa_df$Pop)

p_pa <- ggplot(pa_df, aes(x = Pop, y = Private_Alleles, fill = Pop)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(aes(label = Private_Alleles), vjust = -0.4, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = cols, guide = "none") +
  scale_x_discrete(labels = pop_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.ticks   = element_line(color = "black"),
    panel.grid   = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) +
  labs(
    x = "States",
    y = "Number of Private Alleles"
  )

ggsave(
  filename = file.path(figdir, "private_alleles_barplot.png"),
  plot     = p_pa,
  width    = 8,
  height   = 5,
  dpi      = 300
)

# --- Old Ho vs He barplot (diploid; Ho is undefined for a haploid) ---
# ho_he_long <- data.frame(
#   Pop   = rep(pop_names, 2),
#   Value = c(Ho_pop[pop_names], He_pop[pop_names]),
#   Stat  = rep(c("Ho", "He"), each = length(pop_names))
# )
# p_hohe <- ggplot(ho_he_long, aes(x = Pop, y = Value, fill = Stat)) +
#   geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
#   scale_x_discrete(labels = pop_labels) +
#   scale_fill_manual(values = c("Ho" = "#4393c3", "He" = "#d6604d")) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   labs(x = "Population", y = "Heterozygosity", fill = NULL)
# ggsave(
#   filename = file.path(figdir, "Ho_vs_He_barplot.png"),
#   plot     = p_hohe, width = 8, height = 5, dpi = 300
# )

#############################################
## 7b. Pairwise FST (Weir & Cockerham 1984) + heatmap
#############################################

pairwise_fst_wc <- pairwise.WCfst(hf_dat, diploid = F)

# Ensure diagonal is 0 for filled white square display
diag(pairwise_fst_wc) <- 0

write.csv(as.matrix(pairwise_fst_wc),
          file = file.path(outdir, "pairwise_FST_WC84.csv"),
          quote = FALSE)

fst_wc_df <- as.data.frame(as.table(as.matrix(pairwise_fst_wc)))
colnames(fst_wc_df) <- c("Pop1", "Pop2", "FST")
# Floor negative sampling noise to 0 for color mapping so 0 is pure white
fst_wc_df$FST_fill <- pmax(0, fst_wc_df$FST)
fst_max <- max(fst_wc_df$FST_fill, na.rm = TRUE)
fst_mid <- median(fst_wc_df$FST_fill[fst_wc_df$FST_fill > 0], na.rm = TRUE)

p_fst_wc <- ggplot(fst_wc_df, aes(x = Pop1, y = Pop2, fill = FST_fill)) +
  geom_tile(color = "grey80", linewidth = 0.5) +
  scale_x_discrete(labels = pop_labels) +
  scale_y_discrete(labels = pop_labels) +
  geom_text(aes(label = sprintf("%.3f", FST)),
            size = 4.6, color = "black") +
  scale_fill_gradientn(
    colors = c("white", "orange", "red"),
    values = rescale(c(0, fst_mid, fst_max)),
    limits = c(0, fst_max),
    name   = expression(F[ST])
  ) +
  coord_equal() +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    axis.title    = element_text(size = 12),
    panel.grid    = element_blank()
  ) +
  labs(
    x    = "Population",
    y    = "Population",
    fill = expression(F[ST])
  )

ggsave(
  filename = file.path(figdir, "pairwise_FST_heatmap_WC84.png"),
  plot     = p_fst_wc,
  width    = 10,
  height   = 8,
  dpi      = 300
)


#############################################
## 8. AMOVA (poppr + ade4) — two one-level designs
#############################################

# County is deliberately NOT a stratum: 53 of the 69 State_County units hold a
# single isolate, so that level measures differences between individuals rather
# than between geographic units. Each design is therefore one level deep.
strata(gi_bret) <- data.frame(
  State    = pop_data$State,
  Clusters = pop_data$Clusters
)

## 8a. AMOVA: Clusters only
amova_cl      <- poppr.amova(gi_bret, ~Clusters)
set.seed(1)   # randtest permutes; seed it so the p-values reproduce exactly
amova_cl_rand <- ade4::randtest(amova_cl, nrepet = 999)

sink(file.path(outdir, "amova_clusters_poppr.txt"))
print(amova_cl)
sink()

sink(file.path(outdir, "amova_clusters_randtest.txt"))
print(amova_cl_rand)
sink()

write.csv(as.data.frame(amova_cl$componentsofcovariance),
          file = file.path(outdir, "AMOVA_clusters_components.csv"),
          quote = FALSE)
write.csv(data.frame(obs = amova_cl_rand$obs, pval = amova_cl_rand$pvalue),
          file = file.path(outdir, "AMOVA_clusters_randtest.csv"),
          row.names = FALSE, quote = FALSE)

# PERMANOVA disabled in this release: rScripts/permanova_host.R runs the PERMANOVA (Clusters + Year + host).
# #############################################
# ## 9. PERMANOVA: Clusters & Collection Year (vegan)
# #############################################
# 
# # 1. Parse numeric Year and filter out missing ('-') samples
# pop_data$Year <- as.numeric(pop_data$Collection.year)
# keep          <- !is.na(pop_data$Year)
# meta_sub      <- pop_data[keep, ]
# 
# # 2. Compute individual pairwise genetic distance matrix
# gen_dist    <- diss.dist(gi_bret, percent = FALSE)
# dist_matrix <- as.dist(as.matrix(gen_dist)[keep, keep])
# 
# # 3. Two-factor model: Clusters vs Year (Marginal / Type III effects)
# set.seed(1)
# perm_margin <- adonis2(dist_matrix ~ Clusters + Year, data = meta_sub, permutations = 999, by = "margin")
# 
# # 4. Year tested strictly within Clusters (constrained permutation)
# set.seed(1)
# perm_within <- adonis2(dist_matrix ~ Year, data = meta_sub, strata = meta_sub$Clusters, permutations = 999)
# 
# # 5. Save text outputs to results/
# sink(file.path(outdir, "permanova_clusters_year_margin.txt"))
# cat("=== PERMANOVA: Clusters + Collection Year (Marginal Effects) ===\n")
# print(perm_margin)
# sink()
# 
# sink(file.path(outdir, "permanova_year_within_clusters.txt"))
# cat("=== PERMANOVA: Collection Year Tested Within Clusters (Strata) ===\n")
# print(perm_within)
# sink()
# 
# # 6. Save CSV summaries to results/
# write.csv(as.data.frame(perm_margin),
#           file = file.path(outdir, "permanova_clusters_year_margin.csv"),
#           quote = FALSE)
# write.csv(as.data.frame(perm_within),
#           file = file.path(outdir, "permanova_year_within_clusters.csv"),
#           quote = FALSE)

#############################################
## 10. Nei genetic distance Dst
#############################################

nei_genetic <- genet.dist(hf_dat, method = "Ds", diploid = FALSE)

write.csv(as.matrix(nei_genetic),
          file = file.path(outdir, "nei_genetic.csv"),
          quote = FALSE)

# Long-format for ggplot heatmap
nei_df <- as.data.frame(as.table(as.matrix(nei_genetic)))
colnames(nei_df) <- c("Pop1", "Pop2", "Ds")

p_nei <- ggplot(nei_df, aes(x = Pop1, y = Pop2, fill = Ds)) +
  geom_tile(color = "grey80", linewidth = 0.5) +
  scale_x_discrete(labels = pop_labels) +
  scale_y_discrete(labels = pop_labels) +
  geom_text(aes(label = sprintf("%.3f", Ds)),
            size = 4.6, color = "black") +
  scale_fill_gradient(
    low    = "white",
    high   = "#2166AC",
    limits = c(0, max(nei_df$Ds, na.rm = TRUE)),
    name   = "Nei's Ds"
  ) +
  coord_equal() +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    axis.title   = element_text(size = 12),
    panel.grid   = element_blank()
  ) +
  labs(
    x     = "Population",
    y     = "Population"
  )

ggsave(
  filename = file.path(figdir, "nei_genetic_distance_heatmap.png"),
  plot     = p_nei,
  width    = 10,
  height   = 8,
  dpi      = 300
)

#############################################
## 11. Pairwise FST (Nei) + heatmap
#############################################

# pairwise_fst_nei <- pairwise.neifst(hf_dat, diploid = F)
# 
# # Permutation-based CIs / p-values via boot.ppfst — DISABLED for haploid data.
# # hierfstat::boot.ppfst() errors with diploid = FALSE:
# #   "number of items to replace is not a multiple of replacement length".
# # Its output was never written to any file, so it's commented out to let the run
# # finish. (For haploid pairwise-FST significance, permute individuals among pops.)
# # fst_perm <- boot.ppfst(hf_dat, nboot = 1000, diploid = FALSE)
# # fst_perm$ll   # lower limit of CI
# # fst_perm$ul   # upper limit
# # fst_perm$estimate
# # pval <- ifelse(fst_perm$ll > 0, "<0.05", ">0.05")
# 
# # Save matrix
# write.csv(as.matrix(pairwise_fst_nei),
#           file = file.path(outdir, "pairwise_FST_nei.csv"),
#           quote = FALSE)
# 
# # Long-format for ggplot heatmap
# fst_df_nei <- as.data.frame(as.table(as.matrix(pairwise_fst_nei)))
# colnames(fst_df_nei) <- c("Pop1", "Pop2", "FST")
# 
# # Heatmap
# p_fst_nei <- ggplot(fst_df_nei, aes(x = Pop1, y = Pop2, fill = FST)) +
#   geom_tile(color = "white") +
#   scale_x_discrete(labels = pop_labels) +
#   scale_y_discrete(labels = pop_labels) +
#   geom_text(aes(label = ifelse(is.na(FST), "", sprintf("%.3f", FST))),
#             size = 2.5, color = "black") +
#   scale_fill_gradient2(
#     low  = "white",
#     mid  = "orange",
#     high = "red",
#     midpoint = median(fst_df_nei$FST, na.rm = TRUE),
#     na.value = "white"
#   ) +
#   coord_equal() +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     axis.title.x = element_text(size = 12),
#     axis.title.y = element_text(size = 12)
#   ) +
#   labs(
#     x    = "Population",
#     y    = "Population",
#     fill = expression(F[ST])
#   )
# 
# ggsave(
#   filename = file.path(figdir, "pairwise_FST_heatmap_nei.png"),
#   plot     = p_fst_nei,
#   width    = 10,
#   height   = 8,
#   dpi      = 300
# )

#############################################
## 12. Pairwise Jost's D (mmod) + heatmap
#############################################

message("Calculating pairwise Jost's D...")
pairwise_d_jost <- suppressWarnings(pairwise_D(gi_bret))

# Save matrix to CSV
write.csv(as.matrix(pairwise_d_jost),
          file = file.path(outdir, "pairwise_Jost_D.csv"),
          quote = FALSE)

# Long-format for ggplot heatmap
d_jost_df <- as.data.frame(as.table(as.matrix(pairwise_d_jost)))
colnames(d_jost_df) <- c("Pop1", "Pop2", "Jost_D")
jost_max <- max(d_jost_df$Jost_D, na.rm = TRUE)
jost_mid <- median(d_jost_df$Jost_D[d_jost_df$Jost_D > 0], na.rm = TRUE)

# Heatmap matching the WC84 and Nei aesthetics
p_jost_d <- ggplot(d_jost_df, aes(x = Pop1, y = Pop2, fill = Jost_D)) +
  geom_tile(color = "grey80", linewidth = 0.5) +
  scale_x_discrete(labels = pop_labels) +
  scale_y_discrete(labels = pop_labels) +
  geom_text(aes(label = sprintf("%.3f", Jost_D)),
            size = 4.6, color = "black") +
  scale_fill_gradientn(
    colors = c("#FFFFFF", "#E1BEE7", "#BA68C8"),
    values = rescale(c(0, jost_mid, jost_max)),
    limits = c(0, jost_max),
    name   = expression("Jost's " * italic(D))
  ) +
  coord_equal() +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    axis.title   = element_text(size = 12),
    panel.grid   = element_blank()
  ) +
  labs(
    x    = "Population",
    y    = "Population",
    fill = expression("Jost's " * italic(D))
  )

ggsave(
  filename = file.path(figdir, "pairwise_Jost_D_heatmap.png"),
  plot     = p_jost_d,
  width    = 10,
  height   = 8,
  dpi      = 300
)

save.image("~/Library/CloudStorage/OneDrive-MichiganStateUniversity/SIDE_PROJECTS/karan_bai_popgen/results/stats.RData")
