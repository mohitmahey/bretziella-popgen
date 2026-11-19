rm(list=ls())

# loading required libraries
suppressMessages({ library(pophelper); library(dplyr); library(ggplot2); library(RColorBrewer) })
# Force the cairo raster device: macOS' default quartz png/tiff device colour-manages
# fills (e.g. #4DAF4A renders as #40A33A), dulling the palette. cairo writes true sRGB
# so the bars match the reference / map-pie colours exactly.
options(bitmapType = "cairo")
# Sourcing standard palette to ensure clusters match the map and PCA colours
source("rScripts/color_palette.R")

# ── presets ─────────────────────────────────────────────────────────────────────
k_all     <- 2:20
meanq_dir <- "structure_results/fastStructure/output"
fam_file  <- "structure_results/03_plink/karan_snps_filtered_plink.fam"
meta_file <- "metadata/metadata.csv"
out_dir   <- "figures/structure_plot/pophelper"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── data wrangling ────────────────────────────────────
# V1 = Family, V2 = sample id
fam  <- read.table(fam_file, stringsAsFactors = FALSE)   
meta <- read.csv(meta_file, check.names = FALSE)
names(meta) <- trimws(sub("^﻿", "", names(meta)))
meta <- meta %>% mutate(across(where(is.character), trimws))

# Join meta to fam to keep exact order and get Clusters
fam_meta <- left_join(fam, meta, by = c("V2" = "Sakalidis_lab_code"))

# west -> east state order (by mean longitude)
state_levels <- fam_meta %>%
  group_by(Clusters) %>%
  summarise(lon = mean(as.numeric(Longitude), na.rm = TRUE), .groups = "drop") %>%
  arrange(lon) %>% pull(Clusters)
grp_order <- state_levels

# one-column group-label data.frame (Clusters names, in .fam/.meanQ order)
grp <- data.frame(Clusters = fam_meta$Clusters, stringsAsFactors = FALSE)

# ── pophelper: read + align ─────────────────────────────────────────────────────
qfiles <- sprintf("%s/karan.%d.meanQ", meanq_dir, k_all)             # ascending K
qlist  <- readQ(qfiles, filetype = "basic")
qlist  <- alignK(qlist)                                               # cross-K continuity
names(qlist) <- paste0("K", k_all)

# ── post-alignment: match cluster colours chronologically ─────────────────────
# plotQ uses a single palette array sequentially (col 1 -> pal[1], etc).
# To keep cluster colours stable across K, we forward-propagate column indices
# so inherited clusters stay in the same column, then build the palette at highest K.

# 1. Forward-propagate column indices K->K-1 (greedy correlation match). alignK alone
#    is unreliable across K on this data (sign flips, NAs at high K), so THIS loop is
#    what actually keeps a cluster's colour stable as K increases. Do NOT remove it.
for (i in 2:length(qlist)) {
  cormat <- cor(as.matrix(qlist[[i - 1]]), as.matrix(qlist[[i]]))
  cormat[is.na(cormat)] <- 0
  new_order <- integer(ncol(qlist[[i]]))
  for (j in seq_len(nrow(cormat))) {
    best_col <- which.max(cormat[j, ])
    new_order[j] <- best_col
    cormat[, best_col] <- -Inf    # prevent reuse
  }
  new_order[new_order == 0] <- setdiff(seq_len(ncol(qlist[[i]])), new_order)
  qlist[[i]] <- as.data.frame(qlist[[i]][, new_order, drop = FALSE])
  colnames(qlist[[i]]) <- paste0("Cluster", seq_len(ncol(qlist[[i]])))
}

# 2. Colour each aligned cluster by the population that dominates it at the OPTIMAL K=9,

.map_colors <- .base_colors
# Distinct colours for the higher-K columns (cols >9), assigned positionally.
# Padding heavily so it never runs out of colours for high K values.
.hi_k_colors <- c("#FFE119", "#42D4F4", "#F032E6", "#BFEF45", "#000075",
                  "#800000", "#808000", "#A9A9A9", "#DCBEFF", "#000000", "#C71585",
                  brewer.pal(12, "Set3"), brewer.pal(8, "Pastel2"))
cluster_grp <- fam_meta$Clusters
k_ref <- as.matrix(qlist[["K9"]])                                  # anchor colours at optimal K=9
dom   <- apply(k_ref, 2, function(q) names(which.max(tapply(q, cluster_grp, mean))))
pal   <- rep(NA_character_, ncol(qlist[[length(qlist)]]))          # length = highest K
pal[seq_along(dom)] <- ifelse(duplicated(dom), NA, .map_colors[dom])  # cols 1-9 -> your K=9 palette
pal[is.na(pal)] <- .hi_k_colors[seq_len(sum(is.na(pal)))]         # cols >9 -> distinct hi-K colours
pal   <- unname(pal)

# Set sample IDs as rownames so pophelper's useindlab=TRUE displays them
for (i in seq_along(qlist)) {
  rownames(qlist[[i]]) <- fam$V2
}

# ── pophelper: plot ─────────────────────────────────────────────────────────────
graphics.off()   # macOS quartz device reset before pophelper draws

cat("Drawing joined plots in splits ...\n")
splits <- list(
  list(k_range = 2:10, q_idx = 1:9, name = "K2-10", show_headers = TRUE),
  list(k_range = 11:20, q_idx = 10:19, name = "K11-20", show_headers = TRUE),
  list(k_range = 3:6, q_idx = 2:5, name = "K3-6", show_headers = FALSE),
  list(k_range = 4:8, q_idx = 3:7, name = "K4-8", show_headers = FALSE),
  list(k_range = c(4, 6, 8), q_idx = c(3, 5, 7), name = "K4_6_8", show_headers = FALSE)
)

formats <- c("png", "tiff", "pdf", "rds")

for (s in splits) {
  cat(sprintf("Drawing %s joined plot ...\n", s$name))
  plot_height <- 2.0 * (length(s$k_range) / 9)
  if (plot_height < 2.0) plot_height <- 2.0
  
  for (fmt in formats) {
    fmt_dir <- file.path(out_dir, "sweeps", fmt)
    dir.create(fmt_dir, recursive = TRUE, showWarnings = FALSE)
    
    if (fmt == "rds") {
      p <- plotQ(qlist[s$q_idx],
            exportplot = FALSE, returnplot = TRUE,
            imgoutput = "join",
            grplab = grp, subsetgrp = grp_order, ordergrp = TRUE,
            grplabsize = 6, grplabangle = 45, grplabface = "bold",
            grplabheight = 10, grplabspacer = 0.5, grplabpos = 0.65,
            showdiv = TRUE, divgrp = "Clusters", divcol = "black", divtype = 1, divsize = 0.7,
            barbordersize = 0,
            showsp = TRUE, sppos = "left", splab = paste0("K=", s$k_range),
            splabsize = 9, splabface = "bold", splabcol = "black", spbgcol = "white",
            clustercol = unname(pal),
            showindlab = FALSE, sharedindlab = FALSE, useindlab = FALSE,
            indlabsize = 3, indlabangle = 90, indlabvjust = 0.5,
            indlabhjust = 1, indlabcol = "grey30", indlabspacer = 0.5,
            indlabheight = 1.5,
            showlegend = FALSE, showyaxis = FALSE,
            showtitle = s$show_headers,
            titlelab = if(s$show_headers) paste0("STRUCTURE ancestry sweep, ", sub("-", "..", s$name)) else "",
            titlesize = 12, titleface = "bold",
            showsubtitle = FALSE,
            subtitlelab = "",
            subtitlesize = 8,
            panelspacer = 0.1, panelratio = c(3, 5))
            
      saveRDS(p, file = file.path(fmt_dir, paste0("structure_sweep_", s$name, "_pophelper.rds")))
    } else {
      plotQ(qlist[s$q_idx],
            imgoutput = "join",
            grplab = grp, subsetgrp = grp_order, ordergrp = TRUE,
            grplabsize = 6, grplabangle = 45, grplabface = "bold",
            grplabheight = 10, grplabspacer = 0.5, grplabpos = 0.65,
            showdiv = TRUE, divgrp = "Clusters", divcol = "black", divtype = 1, divsize = 0.7,
            barbordersize = 0,
            showsp = TRUE, sppos = "left", splab = paste0("K=", s$k_range),
            splabsize = 9, splabface = "bold", splabcol = "black", spbgcol = "white",
            clustercol = unname(pal),
            showindlab = FALSE, sharedindlab = FALSE, useindlab = FALSE,
            indlabsize = 3, indlabangle = 90, indlabvjust = 0.5,
            indlabhjust = 1, indlabcol = "grey30", indlabspacer = 0.5,
            indlabheight = 1.5,
            showlegend = FALSE, showyaxis = FALSE,
            showtitle = s$show_headers,
            titlelab = if(s$show_headers) paste0("STRUCTURE ancestry sweep, ", 
                                                 sub("-", "..", s$name)) else "",
            titlesize = 12, titleface = "bold",
            showsubtitle = FALSE,
            subtitlelab = "",
            subtitlesize = 8,
            panelspacer = 0.1, panelratio = c(3, 5),
            height = plot_height, width = 35, dpi = 300, imgtype = fmt,
            exportpath = fmt_dir,
            outputfilename = paste0("structure_sweep_", s$name, "_pophelper"))
    }
  }
}

cat("Drawing separate K plots ...\n")
for (fmt in formats) {
  sep_out_dir <- file.path(out_dir, "individual_K", fmt)
  dir.create(sep_out_dir, showWarnings = FALSE, recursive = TRUE)
  
  if (fmt == "rds") {
    plist <- plotQ(qlist,
          exportplot = FALSE, returnplot = TRUE,
          imgoutput = "sep",
          grplab = grp, subsetgrp = grp_order, ordergrp = TRUE,
          grplabsize = 3.5, grplabangle = 90, grplabface = "bold",
          grplabheight = 28, grplabspacer = 0.3, grplabpos = 0.5,
          showdiv = TRUE, divgrp = "Clusters", divcol = "black", divtype = 1, divsize = 0.7,
          barbordersize = 0,
          clustercol = unname(pal),
          showindlab = FALSE, useindlab = FALSE,
          indlabsize = 3, indlabangle = 90, indlabvjust = 0.5,
          indlabhjust = 1, indlabcol = "grey30", indlabspacer = 0.5,
          indlabheight = 1.5,
          showlegend = FALSE, showyaxis = TRUE, showticks = FALSE, basesize = 14,
          showtitle = TRUE, titlesize = 12, titleface = "bold")$plot   # $plot = one grob per K (K2..K20)

    for(i in seq_along(plist)) {
      saveRDS(plist[[i]], file = file.path(sep_out_dir, paste0("K", i+1, ".rds")))
    }
  } else {
    plotQ(qlist,
          imgoutput = "sep",
          grplab = grp, subsetgrp = grp_order, ordergrp = TRUE,
          grplabsize = 3.5, grplabangle = 90, grplabface = "bold",
          grplabheight = 28, grplabspacer = 0.3, grplabpos = 0.5,
          showdiv = TRUE, divgrp = "Clusters", divcol = "black", divtype = 1, divsize = 0.7,
          barbordersize = 0,
          clustercol = unname(pal),
          showindlab = FALSE, useindlab = FALSE,
          indlabsize = 3, indlabangle = 90, indlabvjust = 0.5,
          indlabhjust = 1, indlabcol = "grey30", indlabspacer = 0.5,
          indlabheight = 1.5,
          showlegend = FALSE, showyaxis = TRUE, showticks = FALSE, basesize = 14,
          showtitle = TRUE, titlesize = 12, titleface = "bold",
          height = 4, width = 80, dpi = 300, imgtype = fmt,
          exportpath = sep_out_dir)
  }
}

cat("Done!\n")

# ══ pophelper-brewer DIAGNOSTIC figures ══════════════════════════════════════════
# A second, self-contained rendering of the same sweeps that intentionally uses
# pophelper's OWN behaviour instead of the custom pipeline above: no forward-propagation
# re-order (l.55-67), no map-matched palette. It re-reads the Q files fresh, uses
# pophelper's alignK for cross-K continuity, and colours clusters with an RColorBrewer
# Set1 + Dark2 + Set3 palette (brewer only -- the same source palettes .map_colors was
# drawn from). Written to its own figures/.../pophelper_brewer_diagnostic/ folder so it
# never touches the production figures above.
cat("Drawing pophelper-brewer diagnostic figures ...\n")
diag_dir <- file.path(out_dir, "pophelper_brewer_diagnostic")
dir.create(diag_dir, showWarnings = FALSE, recursive = TRUE)

nat <- alignK(readQ(qfiles, filetype = "basic"))            # pophelper's own cross-K alignment
names(nat) <- paste0("K", k_all)
brewer_cols <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Dark2"),
                 brewer.pal(12, "Set3"))[seq_len(max(k_all))]   # 29 distinct -> first max(K)

diag_splits <- list(
  list(name = "K2-10",  ks = 2:10),
  list(name = "K10-20", ks = 10:20),
  list(name = "K3-6",   ks = 3:6),
  list(name = "K4-8",   ks = 4:8),
  list(name = "K4_6_8", ks = c(4, 6, 8))
)
for (fmt in formats) {                       # same formats vector as production (png/tiff/pdf/rds)
  fmt_dir <- file.path(diag_dir, fmt)
  dir.create(fmt_dir, showWarnings = FALSE, recursive = TRUE)

  for (s in diag_splits) {
    if (fmt == "rds") {
      p <- plotQ(nat[match(s$ks, k_all)],
            exportplot = FALSE, returnplot = TRUE,
            imgoutput = "join",
            grplab = grp, subsetgrp = grp_order, ordergrp = TRUE,
            grplabsize = 6, grplabangle = 45, grplabface = "bold",
            showsp = TRUE, sppos = "left", splab = paste0("K=", s$ks),
            clustercol = brewer_cols,
            showlegend = FALSE, showyaxis = FALSE)
      saveRDS(p, file = file.path(fmt_dir, paste0("structure_sweep_", s$name, "_brewer.rds")))
    } else {
      plotQ(nat[match(s$ks, k_all)],
            imgoutput = "join",
            grplab = grp, subsetgrp = grp_order, ordergrp = TRUE,
            grplabsize = 6, grplabangle = 45, grplabface = "bold",
            showsp = TRUE, sppos = "left", splab = paste0("K=", s$ks),
            clustercol = brewer_cols,
            showlegend = FALSE, showyaxis = FALSE,
            height = 0.6 * length(s$ks), width = 30, dpi = 150, imgtype = fmt,
            exportpath = fmt_dir,
            outputfilename = paste0("structure_sweep_", s$name, "_brewer"))
    }
  }
}
cat("pophelper-brewer diagnostic figures written to", diag_dir, "\n")

# ── convert all generated PDFs to vector SVG via pdftocairo ────────────────────
if (nzchar(Sys.which("pdftocairo"))) {
  cat("Converting PDF plots to SVG via pdftocairo ...\n")
  all_pdfs <- list.files(out_dir, pattern = "\\.pdf$", full.names = TRUE, recursive = TRUE)
  for (pdf_path in all_pdfs) {
    svg_path <- sub("/pdf/", "/svg/", sub("\\.pdf$", ".svg", pdf_path))
    dir.create(dirname(svg_path), showWarnings = FALSE, recursive = TRUE)
    system2("pdftocairo", c("-svg", shQuote(pdf_path), shQuote(svg_path)))
  }
  cat("SVG conversion complete! SVG files written to:\n",
      file.path(out_dir, "sweeps/svg/"), "\n",
      file.path(out_dir, "individual_K/svg/"), "\n",
      file.path(diag_dir, "svg/"), "\n")
} else {
  warning("pdftocairo not found on PATH; skipped SVG conversion.")
}
