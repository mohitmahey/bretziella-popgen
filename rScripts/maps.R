rm(list = ls())

################################################################################
## STRUCTURE ancestry maps — Oak Wilt (Bretziella fagacearum)                  ##
##                                                                            ##
## 1. US county-level standalone maps (K = 2 to 10)                           ##
##    Output: figures/structure_plot/county_maps/{png,pdf,svg,rds}            ##
## 2. Combined EDDMapS density + Structure county maps (K = 2 to 10):         ##
##    Output: figures/structure_plot/combined_county_maps/{png,pdf,svg,rds}   ##
##    - map_structure_K{K}_density        (without spread, exact coords)      ##
##    - map_structure_K{K}_density_spread (standard spread, zero overlap)     ##
##    - map_structure_K{K}_density_small  (compact spread, zero overlap)      ##
##                                                                            ##
## Displays all 92 isolates on the map:                                       ##
## - 64 Clone-corrected reference isolates: Solid opaque scatterpies (α = 1.0)##
## - 28 Clonal duplicate isolates: Faded translucent scatterpies (α = 0.42)   ##
## Physics repulsion engine guarantees zero overlap across all 92 pies.       ##
## Automatically exports PNG (400 DPI), PDF, RDS, and vector SVG (pdftocairo).##
################################################################################

suppressMessages({
  library(tidyverse)
  library(maps)
  library(scatterpie)
  library(pophelper)
  library(ggnewscale)
})

options(bitmapType = "cairo")
source("rScripts/color_palette.R")

# ── presets ──────────────────────────────────────────────────────────────────
ks          <- 2:10
meanq_dir   <- "structure_results/fastStructure/output"
fam_file    <- "structure_results/03_plink/karan_snps_filtered_plink.fam"
meta_file   <- "metadata/metadata.csv"
cc_file     <- "metadata/clone_corrected_metadata.csv"
spread_file <- "oak_wilt_spread_data/62795/mappings.csv"

dir_standalone <- "figures/structure_plot/county_maps"
dir_combined   <- "figures/structure_plot/combined_county_maps"

for (d in c(dir_standalone, dir_combined)) {
  for (sub in c("png", "pdf", "rds", "svg")) {
    dir.create(file.path(d, sub), recursive = TRUE, showWarnings = FALSE)
  }
}

# ── robust data import ───────────────────────────────────────────────────────
fam <- read.table(fam_file, stringsAsFactors = FALSE)

# Full metadata (92 isolates)
meta_full <- read.csv(meta_file, check.names = FALSE, stringsAsFactors = FALSE)
names(meta_full) <- trimws(sub("^\uFEFF", "", names(meta_full)))
meta_full <- meta_full %>% mutate(across(where(is.character), trimws))
meta_full$County <- gsub("\u00a0", " ", meta_full$County)
meta_full$County <- trimws(meta_full$County)
meta_full$County[meta_full$County == "Sherrune"] <- "Sherburne"
meta_full$Latitude  <- as.numeric(meta_full$Latitude)
meta_full$Longitude <- as.numeric(meta_full$Longitude)

# Clone-corrected metadata (64 isolates)
meta_cc <- read.csv(cc_file, check.names = FALSE, stringsAsFactors = FALSE)
names(meta_cc) <- trimws(sub("^\uFEFF", "", names(meta_cc)))
meta_cc <- meta_cc %>% mutate(across(where(is.character), trimws))
meta_cc$County <- gsub("\u00a0", " ", meta_cc$County)
meta_cc$County <- trimws(meta_cc$County)
meta_cc$County[meta_cc$County == "Sherrune"] <- "Sherburne"
cc_ids  <- meta_cc$Sakalidis_lab_code

meta_full <- meta_full %>%
  mutate(is_clone = !(Sakalidis_lab_code %in% cc_ids))

fam_meta <- left_join(fam, meta_full, by = c("V2" = "Sakalidis_lab_code"))
samples  <- fam$V2
cluster_grp  <- fam_meta$Clusters
state_groups <- sort(unique(cluster_grp))
spare_col    <- "#785EF0"

# ── read Q matrices via pophelper ────────────────────────────────────────────
Qs <- setNames(lapply(ks, function(K) {
  qf <- sprintf("%s/karan.%d.meanQ", meanq_dir, K)
  q  <- readQ(qf, filetype = "basic")[[1]]
  as.matrix(q)
}), as.character(ks))

# ── component colors ─────────────────────────────────────────────────────────
component_colors <- function(K) {
  Q <- Qs[[as.character(K)]]
  mean_q <- matrix(0, nrow = length(state_groups), ncol = ncol(Q))
  rownames(mean_q) <- state_groups
  for (g in state_groups) {
    mask <- cluster_grp == g
    mean_q[g, ] <- colMeans(Q[mask, , drop = FALSE])
  }
  vals <- data.frame(
    state = rep(state_groups, ncol(Q)),
    col   = rep(1:ncol(Q), each = length(state_groups)),
    val   = as.vector(mean_q)
  )
  vals <- vals[order(-vals$val), ]
  comp_state <- rep(NA_character_, K)
  used_s <- character(0); used_c <- integer(0)
  for (r in seq_len(nrow(vals))) {
    s <- vals$state[r]; cc <- vals$col[r]
    if (!(s %in% used_s) && !(cc %in% used_c)) {
      comp_state[cc] <- s; used_s <- c(used_s, s); used_c <- c(used_c, cc)
    }
    if (length(used_c) == K) break
  }
  cols <- rep(spare_col, K)
  ok <- !is.na(comp_state)
  cols[ok] <- unname(get_pop_colors(comp_state[ok]))
  cols
}

comp_colors <- setNames(lapply(ks, component_colors), as.character(ks))

# ── US basemap & scalebar ────────────────────────────────────────────────────
us_states   <- map_data("state")
us_counties <- map_data("county")
lakes       <- map_data("lakes") %>% filter(long >= -125, long <= -65, lat >= 24, lat <= 52)

add_scalebar_usa <- function(x0 = -121, y0 = 25.5, length_km = 1000, lat_ref = 38, h = 0.45) {
  deg_per_km <- 1 / (111.32 * cos(lat_ref * pi / 180))
  len_deg    <- length_km * deg_per_km
  half_deg   <- len_deg / 2

  df_bar1 <- data.frame(xmin = x0, xmax = x0 + half_deg, ymin = y0, ymax = y0 + h)
  df_bar2 <- data.frame(xmin = x0 + half_deg, xmax = x0 + len_deg, ymin = y0, ymax = y0 + h)
  df_labels <- data.frame(
    x = c(x0, x0 + half_deg, x0 + len_deg),
    y = y0 + h + 0.35,
    label = c("0", as.character(length_km / 2), paste0(length_km, " km"))
  )

  list(
    geom_rect(data = df_bar1, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "grey25", color = "grey25", linewidth = 0.25, inherit.aes = FALSE),
    geom_rect(data = df_bar2, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "white", color = "grey25", linewidth = 0.25, inherit.aes = FALSE),
    geom_text(data = df_labels, aes(x = x, y = y, label = label),
              size = 3.3, fontface = "bold", color = "grey20", vjust = 0, inherit.aes = FALSE)
  )
}

# ── 2D Physics Repulsion Engine for 92 Isolates ──────────────────────────────
repel_92 <- function(df, r_vals, max_iter = 650, alpha = 0.16, k_anchor = 0.02, buffer = 1.25) {
  df$r   <- r_vals
  df$lon <- df$lon_init
  df$lat <- df$lat_init
  N <- nrow(df)
  aspect <- 1.3

  for (iter in 1:max_iter) {
    fx <- numeric(N)
    fy <- numeric(N)
    for (i in 1:(N - 1)) {
      for (j in (i + 1):N) {
        dx <- df$lon[i] - df$lon[j]
        dy <- (df$lat[i] - df$lat[j]) * aspect
        d  <- sqrt(dx^2 + dy^2)
        target_d <- (df$r[i] + df$r[j]) * buffer
        if (d < target_d) {
          if (d < 1e-4) {
            dx <- runif(1, -0.05, 0.05)
            dy <- runif(1, -0.05, 0.05)
            d  <- sqrt(dx^2 + dy^2)
          }
          overlap <- target_d - d
          push_x <- (dx / d) * overlap * 0.5
          push_y <- (dy / d) * overlap * 0.5 / aspect
          fx[i] <- fx[i] + push_x
          fy[i] <- fy[i] + push_y
          fx[j] <- fx[j] - push_x
          fy[j] <- fy[j] - push_y
        }
      }
    }
    fx <- fx - (df$lon - df$Longitude) * k_anchor
    fy <- fy - (df$lat - df$Latitude) * k_anchor

    df$lon <- df$lon + fx * alpha
    df$lat <- df$lat + fy * alpha
  }
  df
}

# ── Save helper ──────────────────────────────────────────────────────────────
save_figure <- function(p, fname, target_dir) {
  png_p <- file.path(target_dir, "png", paste0(fname, ".png"))
  pdf_p <- file.path(target_dir, "pdf", paste0(fname, ".pdf"))
  rds_p <- file.path(target_dir, "rds", paste0(fname, ".rds"))
  svg_p <- file.path(target_dir, "svg", paste0(fname, ".svg"))

  ggsave(png_p, p, width = 12, height = 7, dpi = 400, bg = "white")
  ggsave(pdf_p, p, width = 12, height = 7, device = cairo_pdf, bg = "white")
  saveRDS(p, rds_p)
  if (nzchar(Sys.which("pdftocairo"))) {
    system2("pdftocairo", c("-svg", shQuote(pdf_p), shQuote(svg_p)))
  }
  message("Saved: ", fname, " to ", target_dir, "/{png,pdf,rds,svg}/")
}

# ── Load EDDMapS Density Data ─────────────────────────────────────────────────
cat("=== Loading EDDMapS Density Data ===\n")
raw_spread <- read.csv(spread_file, stringsAsFactors = FALSE, check.names = FALSE)
loc_clean  <- trimws(gsub("\"", "", raw_spread$Location))
parts      <- strsplit(loc_clean, ", ")
state_raw   <- sapply(parts, function(x) if (length(x) >= 2) tolower(trimws(x[2])) else NA)
county_raw  <- sapply(parts, function(x) if (length(x) >= 1) tolower(trimws(x[1])) else NA)
country_raw <- sapply(parts, function(x) if (length(x) >= 3) tolower(trimws(x[3])) else NA)

us_mask <- which(country_raw == "united states" & raw_spread$Status != "Negative")
us_data <- data.frame(state_raw = state_raw[us_mask], county_raw = county_raw[us_mask], stringsAsFactors = FALSE)

normalize_county <- function(s, c) {
  c <- gsub("[.]", "", c)
  c <- gsub(" county$", "", c)
  c <- gsub(" parish$", "", c)
  c[c == "dupage"]      <- "du page"
  c[c == "dekalb"]       <- "de kalb"
  c[c == "laporte"]      <- "la porte"
  c[c == "dewitt"]       <- "de witt"
  c[c == "st clair"]     <- "st clair"
  c[c == "st joseph"]    <- "st joseph"
  c[c == "st charles"]   <- "st charles"
  c[c == "st francois"]  <- "st francois"
  c[c == "st louis"]     <- "st louis"
  c <- gsub(" [(]independent city[)]", "", c)
  c[s == "virginia" & c == "bristol"]                       <- "washington"
  c[s == "virginia" & c %in% c("buena vista", "lexington")] <- "rockbridge"
  c[s == "virginia" & c == "covington"]                     <- "alleghany"
  c[s == "virginia" & c == "harrisonburg"]                  <- "rockingham"
  c[s == "virginia" & c == "norton"]                        <- "wise"
  c[s == "virginia" & c %in% c("staunton", "waynesboro")]   <- "augusta"
  c[s == "virginia" & c == "winchester"]                    <- "frederick"
  s[s %in% c("d. columbia", "d columbia")] <- "district of columbia"
  c[s == "district of columbia"]           <- "washington"
  data.frame(state = s, county = c, stringsAsFactors = FALSE)
}

norm_us <- normalize_county(us_data$state_raw, us_data$county_raw)
county_counts <- norm_us %>% count(state, county, name = "n_records")

us_counties_density <- us_counties %>%
  left_join(county_counts, by = c("region" = "state", "subregion" = "county")) %>%
  mutate(
    density_bin = case_when(
      is.na(n_records) | n_records == 0 ~ "No Data",
      n_records == 1                   ~ "1",
      n_records >= 2   & n_records <= 10  ~ "2-10",
      n_records >= 11  & n_records <= 100 ~ "11-100",
      n_records >= 101 & n_records <= 500 ~ "101-500",
      n_records > 500                  ~ "501+",
      TRUE                             ~ "No Data"
    ),
    density_bin = factor(density_bin, levels = c("No Data", "1", "2-10", "11-100", "101-500", "501+"))
  )

density_palette <- c(
  "No Data" = "#FFFFFF",
  "1"       = "#CEE8C3",
  "2-10"    = "#87C27E",
  "11-100"  = "#45894C",
  "101-500" = "#2E6B34",
  "501+"    = "#1A4320"
)

build_density_basemap <- function() {
  ggplot() +
    geom_polygon(data = us_counties_density,
                 aes(x = long, y = lat, group = group, fill = density_bin),
                 colour = "#D1D5DB", linewidth = 0.12) +
    geom_polygon(data = us_states,
                 aes(x = long, y = lat, group = group),
                 fill = NA, colour = "#374151", linewidth = 0.45) +
    geom_polygon(data = lakes,
                 aes(x = long, y = lat, group = group),
                 fill = "#EAF2F8", colour = "#90A4AE", linewidth = 0.25) +
    scale_fill_manual(
      values = density_palette,
      name   = "Observation\nDensity",
      guide  = guide_legend(
        override.aes = list(colour = "#6B7280", linewidth = 0.3),
        keywidth  = unit(1.1, "cm"),
        keyheight = unit(0.55, "cm"),
        order = 1
      )
    ) +
    coord_fixed(1.3, xlim = c(-125, -66.5), ylim = c(24.5, 50)) +
    add_scalebar_usa(x0 = -121, y0 = 25.5, length_km = 1000) +
    theme_void(base_size = 12) +
    theme(
      legend.position        = c(0.90, 0.25),
      legend.background      = element_rect(fill = "white", color = "#D1D5DB", linewidth = 0.3),
      legend.margin          = margin(8, 12, 8, 12),
      legend.title           = element_text(face = "bold", size = 10, lineheight = 1.1),
      legend.text            = element_text(size = 9.0),
      plot.background        = element_rect(fill = "white", color = NA)
    )
}

# ── Precompute Geography & Repulsion Layouts for All 92 Isolates ─────────────
cat("=== Preparing Geography for All 92 Isolates (64 CC + 28 Clones) ===\n")

all_92_geo <- meta_full %>%
  select(sample_id = Sakalidis_lab_code, State, County, Clusters,
         Longitude, Latitude, is_clone) %>%
  group_by(State, County) %>%
  mutate(
    n_in_county = n(),
    order_in_county = row_number(),
    # Rosette initial arrangement so multiple isolates per county start slightly offset
    angle = (order_in_county - 1) * (2 * pi / max(1, n_in_county)),
    rosette_r = ifelse(n_in_county > 1, 0.22, 0),
    lon_init = Longitude + rosette_r * cos(angle) * 1.3,
    lat_init = Latitude  + rosette_r * sin(angle),
    loc = sample_id,
    r_std  = 0.36,
    r_comp = 0.24
  ) %>%
  ungroup()

cat(sprintf("Configured %d isolates: %d CC (solid) + %d clones (faded)\n",
            nrow(all_92_geo), sum(!all_92_geo$is_clone), sum(all_92_geo$is_clone)))

# Layout 1: Exact coordinates (with county rosette)
layout_exact <- all_92_geo %>%
  mutate(lon = lon_init, lat = lat_init, r = r_std)

# Layout 2: Standard spread layout (zero overlap across all 92 pies)
cat("Calculating standard pie repulsion layout (92 isolates, zero overlap)...\n")
layout_spread <- repel_92(all_92_geo, all_92_geo$r_std,
                          max_iter = 650, alpha = 0.16, k_anchor = 0.02, buffer = 1.25)

# Layout 3: Compact spread layout (zero overlap, compact pies)
cat("Calculating compact pie repulsion layout (92 isolates, zero overlap)...\n")
layout_small <- repel_92(all_92_geo, all_92_geo$r_comp,
                         max_iter = 550, alpha = 0.16, k_anchor = 0.03, buffer = 1.30)

# Helper function to plot 92 scatterpies (faded clones underneath, solid CC on top)
add_92_pie_layers <- function(p, df_layout, kcols, K) {
  df_cc    <- df_layout %>% filter(!is_clone)
  df_clone <- df_layout %>% filter(is_clone)

  p <- p +
    new_scale_fill() +
    # 1a. Faded clonal pies: soft halo underlay
    geom_scatterpie(data = df_clone %>% mutate(r_halo = r * 1.08),
                    aes(x = lon, y = lat, group = loc, r = r_halo),
                    cols = paste0("C", 1:K), colour = "#FFFFFF", fill = "#FFFFFF",
                    linewidth = 0.35, alpha = 0.60) +
    # 1b. Faded clonal pies: translucent pies (alpha = 0.42)
    geom_scatterpie(data = df_clone,
                    aes(x = lon, y = lat, group = loc, r = r),
                    cols = paste0("C", 1:K), colour = "#4B5563", linewidth = 0.20,
                    alpha = 0.42) +
    scale_fill_manual(values = kcols, guide = "none") +
    new_scale_fill() +
    # 2a. Solid CC pies: white halo underlay
    geom_scatterpie(data = df_cc %>% mutate(r_halo = r * 1.08),
                    aes(x = lon, y = lat, group = loc, r = r_halo),
                    cols = paste0("C", 1:K), colour = "#FFFFFF", fill = "#FFFFFF",
                    linewidth = 0.45, alpha = 1.0) +
    # 2b. Solid CC pies: opaque pies (alpha = 1.0)
    geom_scatterpie(data = df_cc,
                    aes(x = lon, y = lat, group = loc, r = r),
                    cols = paste0("C", 1:K), colour = "#111111", linewidth = 0.22,
                    alpha = 1.0) +
    scale_fill_manual(values = kcols, guide = "none")

  p
}

# ── Render All Maps for K = 2 to 10 ──────────────────────────────────────────
density_basemap <- build_density_basemap()

for (K in ks) {
  cat(sprintf("\n=== Processing All-Sample Maps for K = %d (N = 92) ===\n", K))
  Kchr <- as.character(K)
  Q    <- Qs[[Kchr]]
  qdf  <- as.data.frame(Q); names(qdf) <- paste0("C", 1:K); qdf$sample_id <- samples
  qdf  <- left_join(qdf, fam_meta, by = c("sample_id" = "V2"))

  # 1. 64 CC samples individual Q
  cc_q <- qdf %>%
    transmute(sample_id = sample_id, Clusters = Clusters, across(all_of(paste0("C", 1:K))))

  # 2. Cluster mean ancestry for the 28 clonal duplicates
  cluster_q_means <- cc_q %>%
    group_by(Clusters) %>%
    summarise(across(all_of(paste0("C", 1:K)), mean), .groups = "drop")

  clone_q <- all_92_geo %>%
    filter(is_clone) %>%
    select(sample_id, Clusters) %>%
    left_join(cluster_q_means, by = "Clusters")

  all_92_q <- bind_rows(cc_q, clone_q)
  kcols <- setNames(comp_colors[[Kchr]], paste0("C", 1:K))

  # ── 1. Standalone Map (Standard spread, zero overlap) ──
  df_standalone_k <- left_join(layout_spread, all_92_q, by = c("sample_id", "Clusters"))

  p_standalone_base <- ggplot() +
    geom_polygon(data = us_states, aes(long, lat, group = group),
                 fill = "#FFFFFF", colour = NA) +
    geom_polygon(data = us_counties, aes(long, lat, group = group),
                 fill = NA, colour = "#E5E7EB", linewidth = 0.12) +
    geom_polygon(data = us_states, aes(long, lat, group = group),
                 fill = NA, colour = "#4B5563", linewidth = 0.45) +
    geom_polygon(data = lakes, aes(long, lat, group = group),
                 fill = "#EAF2F8", colour = "#90A4AE", linewidth = 0.25)

  p_standalone <- add_92_pie_layers(p_standalone_base, df_standalone_k, kcols, K) +
    coord_fixed(1.3, xlim = c(-125, -66.5), ylim = c(24.5, 50)) +
    add_scalebar_usa(x0 = -121, y0 = 25.5, length_km = 1000) +
    theme_void(base_size = 12) +
    theme(legend.position = "none", plot.background = element_rect(fill = "white", color = NA))

  save_figure(p_standalone, sprintf("map_structure_K%d_USA", K), dir_standalone)

  # ── 2. Combined: Without Spread (Exact / Rosette Coords) ──
  df_exact_k <- left_join(layout_exact, all_92_q, by = c("sample_id", "Clusters"))
  p_density  <- add_92_pie_layers(density_basemap, df_exact_k, kcols, K)
  save_figure(p_density, sprintf("map_structure_K%d_density", K), dir_combined)

  # ── 3. Combined: With Spread (Standard Pies, Zero Overlap) ──
  df_spread_k      <- left_join(layout_spread, all_92_q, by = c("sample_id", "Clusters"))
  p_density_spread <- add_92_pie_layers(density_basemap, df_spread_k, kcols, K)
  save_figure(p_density_spread, sprintf("map_structure_K%d_density_spread", K), dir_combined)

  # ── 4. Combined: With Small Spread (Compact Pies, Zero Overlap) ──
  df_small_k      <- left_join(layout_small, all_92_q, by = c("sample_id", "Clusters"))
  p_density_small <- add_92_pie_layers(density_basemap, df_small_k, kcols, K)
  save_figure(p_density_small, sprintf("map_structure_K%d_density_small", K), dir_combined)
}

cat("\n=== All Maps (Standalone & Combined Density K=2..10, N=92 Isolates) Complete! ===\n")
