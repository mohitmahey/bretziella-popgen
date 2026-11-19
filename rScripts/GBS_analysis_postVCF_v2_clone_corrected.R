# GBS analysis of the Bretziella fagacerium
rm(list = ls())
set.seed(1)

# loading required libraries
library(vcfR)
library(tidyverse)
library(poppr)
library(ape)
library(RColorBrewer)
library(igraph)
library(ggplot2)
library(adegenet)

# setting up working directory
setwd(".")

# 1. Load Data ----
message("Loading VCF and Metadata...")
bret_vcf <- read.vcfR("vcf/karan_biallelic_TAGs_snps_miss0.1_MAF_LD_CC.vcf",
                      verbose = FALSE)
pop_data <- read.csv("metadata/clone_corrected_metadata.csv", header = TRUE,
                     stringsAsFactors = FALSE)
# removing trailing white spaces from front and back
pop_data <- pop_data %>% mutate(across(where(is.character), trimws))

# Match metadata to VCF
vcf_samples <- colnames(bret_vcf@gt)[-1]
pop_data <- pop_data[match(vcf_samples, pop_data$Sakalidis_lab_code), ]

# Convert to Genlight
gl_bret <- vcfR2genlight(bret_vcf)
gi_bret <- vcfR2genind(bret_vcf,return.alleles = TRUE)

ploidy(gl_bret) <- 1
ploidy(gi_bret) <- 1

pop(gi_bret) <- pop_data$Clusters

gl_bret
# 3. Dynamic Colour & Shape Palette (Sourced from central palette) ----
source("rScripts/color_palette.R")

# Group by the `Clusters` column (9 groups; the six eastern states MD/NY/OH/
# PA/SC/WV are merged into "Eastern").
pops_ordered <- sort(unique(pop_data$Clusters))
pop(gl_bret) <- factor(pop_data$Clusters, levels = pops_ordered)
cols <- get_pop_colors(pops_ordered)
pop_shapes <- get_pop_shapes(pops_ordered)
pop_labels <- get_pop_labels(pops_ordered)

# 4. Publication Theme ----
pub_theme <- theme_minimal(base_family = "Times New Roman", base_size = 12) +
  theme(
    panel.border     = element_rect(colour = "grey30", fill = NA, linewidth = 0.5),
    axis.line        = element_blank(),
    panel.background = element_blank(),
    legend.title     = element_text(size = 15, face = "bold"),
    legend.text      = element_text(size = 12),
    plot.title       = element_text(face = "bold", size = 14)
  )

# 5. Analysis 1: UPGMA Tree ----
message("Building UPGMA Tree...")
tree <- aboot(gl_bret, tree = "upgma", 
              distance = bitwise.dist, sample = 100, 
              showtree = F, 
              cutoff = 50)

# writing tree as output for figtree

# to block NA
tree$node.label[is.na(tree$node.label)] <- ""
ape::write.tree(tree, "results/genetic_tree.nwk")

# annotation for fig tree
ann <- data.frame(taxa = tree$tip.label,
                  Cluster = pop_data$Clusters[match(tree$tip.label, pop_data$Sakalidis_lab_code)])
# to get exact colors
ann$Color <- cols[match(ann$Cluster, pops_ordered)]
write.table(ann, "results/genetic_tree_annotations.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
# colour every branch: a node takes its cluster colour if all descendants agree
n_tip   <- length(tree$tip.label)
node_cl <- rep(NA_character_, n_tip + tree$Nnode)
node_cl[1:n_tip] <- ann$Cluster

po <- reorder(tree, "postorder")
for (i in seq_len(nrow(po$edge))) {
  p <- po$edge[i, 1]; ch <- po$edge[i, 2]
  v <- node_cl[ch]
  node_cl[p] <- if (is.na(node_cl[p])) v else if (identical(node_cl[p], v)) v else "MIXED"
}
node_col <- ifelse(node_cl == "MIXED", "#BBBBBB", cols[match(node_cl, pops_ordered)])

t2 <- tree
t2$tip.label  <- paste0(tree$tip.label,  "[&!color=", node_col[1:n_tip],  "]")
t2$node.label <- paste0(tree$node.label, "[&!color=", node_col[-(1:n_tip)], "]")

writeLines(c("#NEXUS", "begin trees;",
             paste("\ttree t1 = [&R]", write.tree(t2)), "end;"),
           "results/genetic_tree_coloured.nex")

# Tree tips keep the VCF sample IDs (BF codes); no separate display-name column.
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
png(file = "figures/genetic_tree.png", width = 7, height = 10, units = "in",
    res = 600, type = "cairo")
par(family = "sans", mar = c(4, 0.5, 0.5, 0.5), lwd = 0.5)
plot.phylo(tree, cex = 0.45, tip.color = cols[pop(gl_bret)], 
           label.offset = 0.003, edge.width = 0.6)

legend("topleft", legend = pop_labels[pops_ordered],
       fill = cols, border = NA, bty = "n"
       , cex = 0.6, title = "Population")
dev.off()

# 6. Analysis 2: PCA ----
message("Running PCA...")
# nf keeps every axis: section 6b needs the full set to form the denominator of
# its variance percentages. Eigenvalues and the PC1/PC2 coordinates are unchanged.
rubi.pca <- glPca(gl_bret, parallel = FALSE, nf = nInd(gl_bret) - 1)
pca_scores <- as.data.frame(rubi.pca$scores)
pca_scores$pop <- pop(gl_bret)

var_exp <- rubi.pca$eig / sum(rubi.pca$eig) * 100
pc1_lab <- paste0("PC1 (", round(var_exp[1], 1), "%)")
pc2_lab <- paste0("PC2 (", round(var_exp[2], 1), "%)")

# Save master PCA
p_pca_master <- ggplot(pca_scores, aes(x = PC1, y = PC2, colour = pop, shape = pop)) +
  stat_ellipse(aes(fill = pop), geom = "polygon", level = 0.95, alpha = 0.08,
               color = NA, show.legend = FALSE) +
  stat_ellipse(level = 0.95, linewidth = 0.9, show.legend = FALSE) +  # highlighted outline (per-cluster colour)
  geom_point(size = 5, alpha = 0.80) +
  scale_colour_manual(values = cols, name = "Population", labels = pop_labels) +
  scale_fill_manual(values = cols) +
  scale_shape_manual(values = pop_shapes, name = "Population", 
                     labels = pop_labels) +
  labs(x = pc1_lab, y = pc2_lab) +
  pub_theme +
  theme(axis.title   = element_text(size = 36),
         axis.text    = element_text(size = 26),
         legend.title = element_text(size = 24, face = "bold"),   # "Population" title
         legend.text  = element_text(size = 20)) +                # state names
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1)))  # symbols, was 2.5
p_pca_master

ggsave("figures/PCA_Master_v2.png", p_pca_master, width = 9, height = 6, dpi = 600)

  # 6b. Analysis 2b: Zoomed PCA — MI, TX and Eastern dropped ----
# variation among the six that remain.
message("Plotting zoomed PCA...")
zoom_drop <- c("MI", "TX", "Eastern")
zoom_keep <- !pca_scores$pop %in% zoom_drop

zoom_var <- apply(pca_scores[zoom_keep, grep("^PC", names(pca_scores))], 2, var)
zoom_pct <- 100 * zoom_var / sum(zoom_var)

# Adding a data frame swaps the data into the master plot, so theme, palette,
# shapes and both ellipse layers stay in sync with it. droplevels() keeps the
# dropped groups out of the legend.
p_pca_zoom <- p_pca_master + droplevels(pca_scores[zoom_keep, ]) +
  labs(x = paste0("PC1 (", round(zoom_pct[1], 1), "%)"),
       y = paste0("PC2 (", round(zoom_pct[2], 1), "%)"))

# AR (n=1) and IL (n=2) are too small for a 95% ellipse; their points plot bare
# and stat_ellipse prints a harmless "Too few points" message for each.
ggsave("figures/PCA_Zoom_v2.png", p_pca_zoom, width = 8, height = 6, dpi = 600)

 # 6c. PCA coloured by host sub-section ----
# Same PCA, same coordinates and look as the master; ellipses and colour = host,
# shape = population. Chestnut (n = 2) is too small for an ellipse.
pca_scores$host <- pop_data$Sub.section
p_pca_host <- ggplot(pca_scores, aes(x = PC1, y = PC2, colour = host)) +
  stat_ellipse(aes(fill = host), geom = "polygon", level = 0.95, alpha = 0.08,
               colour = NA, show.legend = FALSE) +
  stat_ellipse(level = 0.95, linewidth = 0.9, show.legend = FALSE) +
  geom_point(aes(shape = pop), size = 5, alpha = 0.80) +
  scale_colour_manual(values = host_colors, name = "Host",
                      labels = c(UA = "Unknown")) +
  scale_fill_manual(values = host_colors) +
  scale_shape_manual(values = pop_shapes, name = "Population", labels = pop_labels) +
  labs(x = pc1_lab, y = pc2_lab) +
  pub_theme +
  theme(axis.title   = element_text(size = 36),
        axis.text    = element_text(size = 26),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text  = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1)),
         shape  = guide_legend(override.aes = list(size = 5)))
p_pca_host

ggsave("figures/PCA_Host_v2.png", p_pca_host, width = 8, height = 6, dpi = 600)

# 7. Analysis 3: DAPC ----
message("Running DAPC...")

# Retention based on previous optimization (xvalDapc)
pnw.dapc <- dapc(gl_bret, n.pca = 8, n.da = 4)

# # optimum number of PCs to retain, run this just first time
# png(file = "figures/dapc_optimum_score.png",
#     width = 7, height = 7, units = "in",
#     res = 600, type = "cairo")
# optim.a.score(pnw.dapc)
# dev.off()

# scatter(pnw.dapc, col = cols, cex = 2, legend = TRUE,
#         clabel = F, posi.leg = "bottomleft", scree.pca = TRUE,
#         posi.pca = "topleft", cleg = 0.75)

# Get DAPC scores for ggplot
dapc_scores <- as.data.frame(pnw.dapc$ind.coord)
dapc_scores$pop <- pop(gl_bret)

# Calculate centroids for star-plot look (emulating adegenet::scatter)
centroids_dapc <- dapc_scores %>%
  group_by(pop) %>%
  summarize(LD1 = mean(LD1), LD2 = mean(LD2))

dapc_scores <- dapc_scores %>%
  left_join(centroids_dapc, by = "pop", suffix = c("", ".cent"))

# Variance explained by LDs
dapc_var <- pnw.dapc$eig / sum(pnw.dapc$eig) * 100
ld1_lab <- paste0("LD1 (", round(dapc_var[1], 1), "%)")
ld2_lab <- paste0("LD2 (", round(dapc_var[2], 1), "%)")

# Helper: Scree plot inset
dapc_scree <- function(dapc_obj) {
  eig_df <- data.frame(LD = 1:length(dapc_obj$eig), Eig = dapc_obj$eig)
  ggplot(eig_df, aes(x = LD, y = Eig)) +
    geom_bar(stat = "identity", fill = "grey70") +
    theme_void() +
    theme(panel.background = element_rect(fill = "white", color = "grey80"))
}

# Master DAPC (with star-lines)
p_dapc_master <- ggplot(dapc_scores, aes(x = LD1, y = LD2, color = pop, shape = pop)) +
  geom_segment(aes(xend = LD1.cent, yend = LD2.cent), alpha = 0.15, linewidth = 0.2) +
  stat_ellipse(aes(fill = pop), geom = "polygon", level = 0.95, alpha = 0.08, color = NA, show.legend = FALSE) +
  stat_ellipse(level = 0.95, linewidth = 0.9, show.legend = FALSE) +  # highlighted outline (per-cluster colour)
  geom_point(size = 5, alpha = 0.80) +
  scale_color_manual(values = cols, name = "Population", labels = pop_labels) +
  scale_fill_manual(values = cols) +
  scale_shape_manual(values = pop_shapes, name = "Population", labels = pop_labels) +
  labs(x = ld1_lab, y = ld2_lab) +
  pub_theme +
  theme(axis.title   = element_text(size = 36),
        axis.text    = element_text(size = 26),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text  = element_text(size = 20)) +
  guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))
p_dapc_master
ggsave("figures/DAPC_Master_v2.png", p_dapc_master, width = 7, height = 5, dpi = 600)

# 7b. DAPC by host sub-section ----
# Same SNPs and n.pca as the master DAPC; groups = host instead of population.
# 4 host groups, so at most 3 discriminant axes.
host_grp  <- factor(pop_data$Sub.section)
host.dapc <- dapc(gl_bret, pop = host_grp, n.pca = 8, n.da = 3)

dapc_host <- as.data.frame(host.dapc$ind.coord)
dapc_host$host <- host_grp
dapc_host$pop  <- pop(gl_bret)

host_var <- host.dapc$eig / sum(host.dapc$eig) * 100

p_dapc_host <- ggplot(dapc_host, aes(x = LD1, y = LD2, colour = host)) +
  stat_ellipse(aes(fill = host), geom = "polygon", level = 0.95, alpha = 0.08,
               colour = NA, show.legend = FALSE) +
  stat_ellipse(level = 0.95, linewidth = 0.9, show.legend = FALSE) +
  geom_point(aes(shape = pop), size = 5, alpha = 0.80) +
  scale_colour_manual(values = host_colors, name = "Host",
                      labels = c(UA = "Unknown")) +
  scale_fill_manual(values = host_colors) +
  scale_shape_manual(values = pop_shapes, name = "Population", labels = pop_labels) +
  labs(x = paste0("LD1 (", round(host_var[1], 1), "%)"),
       y = paste0("LD2 (", round(host_var[2], 1), "%)")) +
  pub_theme +
  theme(axis.title   = element_text(size = 36),
        axis.text    = element_text(size = 26),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text  = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1)),
         shape  = guide_legend(override.aes = list(size = 5)))
p_dapc_host

ggsave("figures/DAPC_Host_v2.png", p_dapc_host, width = 7, height = 6, dpi = 600)

# 8. DAPC Membership Barplot ----
message("Generating DAPC membership barplot...")
dapc.results <- as.data.frame(pnw.dapc$posterior)
dapc.results$pop <- pop(gl_bret)
dapc.results$Sample <- rownames(dapc.results)
dapc.results <- pivot_longer(dapc.results, -c(pop, Sample), names_to = "Assigned_Pop", values_to = "Probability")

p_mem <- ggplot(dapc.results, aes(x=Sample, y=Probability, fill=Assigned_Pop)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values = cols, labels = pop_labels) +
  facet_grid(~pop, scales = "free", space = "free",
             labeller = as_labeller(pop_labels)) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        panel.spacing = unit(0.1, "lines"), strip.text.x = element_text(angle = 90, size = 6)) +
  labs(y = "Posterior Probability", x = "Individual Samples")

ggsave("figures/DAPC_Membership_v2.png", p_mem, width = 12, height = 5, dpi = 600)
saveRDS(p_mem, "figures/DAPC_Membership_v2.rds")

# Save results
dir.create("results", recursive = TRUE, showWarnings = FALSE)
saveRDS(pnw.dapc, "results/dapc_object_v2.rds")
message("DAPC Analysis Complete. Plots in figures/")

find.clusters(gl_bret)

# Inspect cluster sizes
save.image("~/Library/CloudStorage/OneDrive-MichiganStateUniversity/SIDE_PROJECTS/karan_bai_popgen/GBS_analysis_postVCF_v2.RData")

#############################################
## 9. poppr diversity indices
#############################################

poppr_res <- poppr(gi_bret)
poppr_df  <- as.data.frame(poppr_res)
outdir <- "results"
write.csv(poppr_df,
          file = file.path(outdir, "poppr_diversity_summary.csv"),
          row.names = FALSE,
          quote = FALSE)

