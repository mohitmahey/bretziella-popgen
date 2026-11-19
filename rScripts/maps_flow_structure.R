rm(list = ls())

################################################################################
## STRUCTURE ancestry map — Oak Wilt (Bretziella fagacearum)                   ##
##                                                                            ##
## 1. US cluster-average standalone maps (K = 2 to 10)                        ##
##    Output: figures/structure_plot/average_maps/{png,pdf,svg,rds}           ##
## 2. Combined EDDMapS density + Structure cluster average maps (K = 2 to 10): ##
##    Output: figures/structure_plot/combined_average_maps/{png,pdf,svg,rds}  ##
##    - average_map_structure_K{K}_density        (without spread)            ##
##    - average_map_structure_K{K}_density_spread (with spread, standard pies)##
##    - average_map_structure_K{K}_density_small  (with spread, smaller pies) ##
##                                                                            ##
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
spread_file <- "oak_wilt_spread_data/62795/mappings.csv"

dir_standalone <- "figures/structure_plot/average_maps"
dir_combined   <- "figures/structure_plot/combined_average_maps"

for (d in c(dir_standalone, dir_combined)) {
  for (sub in c("png", "pdf", "rds", "svg")) {
    dir.create(file.path(d, sub), recursive = TRUE, showWarnings = FALSE)
  }
}

# ── Scale bar helper for Full USA (1000 km, 0 - 500 - 1000 km) ────────────────
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

# ── robust data import ──────────────────────────────────────────────────────────
fam  <- read.table(fam_file, stringsAsFactors = FALSE)
meta <- read.csv(meta_file, check.names = FALSE)
names(meta) <- trimws(sub("^\uFEFF", "", names(meta)))
names(meta) <- trimws(names(meta))
meta <- meta %>% mutate(across(where(is.character), trimws))
meta$County[meta$County == "Sherrune"] <- "Sherburne"

fam_meta <- left_join(fam, meta, by = c("V2" = "Sakalidis_lab_code"))
samples  <- fam$V2
cluster_grp  <- fam_meta$Clusters
state_groups <- sort(unique(cluster_grp))
spare_col    <- "#785EF0"

# ── Manual declutter nudges for cluster average maps ───────────────────────────
nudge <- tibble::tribble(
  ~grp,       ~dlon, ~dlat,  ~lx,   ~ly,  ~hj, ~vj,
  "MN",         0.0,   0.0,  0.0,   1.8,  0.5,  0,
  "WI",         0.0,   0.0,  0.0,   1.8,  0.5,  0,
  "MI",         0.0,   0.0,  0.0,   1.9,  0.5,  0,
  "IA",        -0.3,   0.2, -1.8,   0.0,  1.0,  0.5,
  "MO",        -1.0,  -0.9, -1.8,   0.0,  1.0,  0.5,
  "IL",         0.0,   0.0,  0.0,  -1.4,  0.5,  1,
  "AR",         0.0,   0.0,  0.0,  -1.3,  0.5,  1,
  "TX",         0.0,   0.0,  0.0,   1.8,  0.5,  0,
  "Eastern",    0.0,   0.0,  0.0,   2.2,  0.5,  0
)

us          <- map_data("state")
us_counties <- map_data("county")
lakes       <- map_data("lakes") %>% filter(long >= -125, long <= -65, lat >= 24, lat <= 52)

# ── Repulsion engine for cluster average pies ─────────────────────────────────
repel_cluster_pies <- function(df, r_vals, max_iter = 350, alpha = 0.16, k_anchor = 0.05, buffer = 1.10) {
  df$r   <- r_vals
  df$cur_lon <- df$lon
  df$cur_lat <- df$lat
  N <- nrow(df)
  aspect <- 1.3

  for (iter in 1:max_iter) {
    fx <- numeric(N)
    fy <- numeric(N)
    for (i in 1:(N - 1)) {
      for (j in (i + 1):N) {
        dx <- df$cur_lon[i] - df$cur_lon[j]
        dy <- (df$cur_lat[i] - df$cur_lat[j]) * aspect
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
    fx <- fx - (df$cur_lon - df$lon) * k_anchor
    fy <- fy - (df$cur_lat - df$lat) * k_anchor

    df$cur_lon <- df$cur_lon + fx * alpha
    df$cur_lat <- df$cur_lat + fy * alpha
  }
  df$lon <- df$cur_lon
  df$lat <- df$cur_lat
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
    geom_polygon(data = us,
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

# ── Precompute Cluster Baseline & Repulsion Layouts ───────────────────────────
info_base <- fam_meta %>%
  transmute(sample_id = V2, grp = Clusters,
            lon = as.numeric(Longitude), lat = as.numeric(Latitude))

grp_base <- info_base %>%
  group_by(grp) %>%
  summarise(lon = mean(lon), lat = mean(lat), n = dplyr::n(), .groups = "drop") %>%
  left_join(nudge, by = "grp") %>%
  mutate(across(c(dlon, dlat, lx, ly), ~ coalesce(., 0)),
         hj = coalesce(hj, 0.5), vj = coalesce(vj, 0.5),
         lon = lon + dlon, lat = lat + dlat,
         label = get_pop_labels(grp),
         lab_lon = lon + lx, lab_lat = lat + ly)

# 1. Exact coords (standard sizing)
r_std_avg <- 0.55 + 0.15 * sqrt(grp_base$n)
layout_exact <- grp_base %>% mutate(r = r_std_avg)

# 2. Standard spread layout
cat("Calculating standard cluster pie repulsion layout...\n")
layout_spread <- repel_cluster_pies(grp_base, r_std_avg, max_iter = 350, alpha = 0.16, k_anchor = 0.05)

# 3. Small spread layout
cat("Calculating small cluster pie repulsion layout...\n")
r_small_avg <- 0.35 + 0.09 * sqrt(grp_base$n)
layout_small <- repel_cluster_pies(grp_base, r_small_avg, max_iter = 300, alpha = 0.14, k_anchor = 0.06)

# ── Render All Maps for K = 2 to 10 ──────────────────────────────────────────
density_basemap <- build_density_basemap()

for (K in ks) {
  qf <- sprintf("%s/karan.%d.meanQ", meanq_dir, K)
  Q  <- as.matrix(readQ(qf, filetype = "basic")[[1]])

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
  ccols <- rep(spare_col, K)
  ok <- !is.na(comp_state)
  ccols[ok] <- unname(get_pop_colors(comp_state[ok]))
  ccols <- setNames(ccols, paste0("C", 1:K))

  qdf  <- as.data.frame(Q); names(qdf) <- paste0("C", 1:K); qdf$sample_id <- samples
  qdf  <- left_join(qdf, info_base, by = "sample_id")

  # Calculate cluster mean ancestry proportions
  grp_q <- qdf %>%
    group_by(grp) %>%
    summarise(across(all_of(paste0("C", 1:K)), mean), .groups = "drop")

  # ── 1. Standalone Cluster Average Map ──
  df_standalone <- left_join(layout_exact, grp_q, by = "grp")

  p_standalone <- ggplot() +
    geom_polygon(data = us, aes(long, lat, group = group),
                 fill = "#FFFFFF", colour = NA) +
    geom_polygon(data = us_counties, aes(long, lat, group = group),
                 fill = NA, colour = "#E5E7EB", linewidth = 0.12) +
    geom_polygon(data = us, aes(long, lat, group = group),
                 fill = NA, colour = "#4B5563", linewidth = 0.45) +
    geom_scatterpie(data = df_standalone,
                    aes(x = lon, y = lat, r = r, group = grp),
                    cols = paste0("C", 1:K), colour = "#1F2937", linewidth = 0.25) +
    scale_fill_manual(values = ccols, guide = "none") +
    coord_fixed(1.3, xlim = c(-125, -66.5), ylim = c(24.5, 50)) +
    add_scalebar_usa(x0 = -121, y0 = 25.5, length_km = 1000) +
    theme_void(base_size = 12) +
    theme(legend.position = "none")

  save_figure(p_standalone, sprintf("average_map_structure_K%d_USA", K), dir_standalone)

  # ── 2. Combined: Without Spread (Exact Coords + Halo) ──
  df_exact_k <- left_join(layout_exact, grp_q, by = "grp")

  p_avg_density <- density_basemap +
    new_scale_fill() +
    geom_scatterpie(data = df_exact_k %>% mutate(r_halo = r * 1.06),
                    aes(x = lon, y = lat, r = r_halo, group = grp),
                    cols = paste0("C", 1:K), colour = "#FFFFFF", fill = "#FFFFFF", linewidth = 0.50) +
    geom_scatterpie(data = df_exact_k,
                    aes(x = lon, y = lat, r = r, group = grp),
                    cols = paste0("C", 1:K), colour = "#111111", linewidth = 0.30) +
    scale_fill_manual(values = ccols, guide = "none")

  save_figure(p_avg_density, sprintf("average_map_structure_K%d_density", K), dir_combined)

  # ── 3. Combined: With Spread (Standard Pies + Halo) ──
  df_spread_k <- left_join(layout_spread, grp_q, by = "grp")

  p_avg_density_spread <- density_basemap +
    new_scale_fill() +
    geom_scatterpie(data = df_spread_k %>% mutate(r_halo = r * 1.06),
                    aes(x = lon, y = lat, r = r_halo, group = grp),
                    cols = paste0("C", 1:K), colour = "#FFFFFF", fill = "#FFFFFF", linewidth = 0.50) +
    geom_scatterpie(data = df_spread_k,
                    aes(x = lon, y = lat, r = r, group = grp),
                    cols = paste0("C", 1:K), colour = "#111111", linewidth = 0.30) +
    scale_fill_manual(values = ccols, guide = "none")

  save_figure(p_avg_density_spread, sprintf("average_map_structure_K%d_density_spread", K), dir_combined)

  # ── 4. Combined: With Small Spread (Compact Pies + Halo) ──
  df_small_k <- left_join(layout_small, grp_q, by = "grp")

  p_avg_density_small <- density_basemap +
    new_scale_fill() +
    geom_scatterpie(data = df_small_k %>% mutate(r_halo = r * 1.08),
                    aes(x = lon, y = lat, r = r_halo, group = grp),
                    cols = paste0("C", 1:K), colour = "#FFFFFF", fill = "#FFFFFF", linewidth = 0.40) +
    geom_scatterpie(data = df_small_k,
                    aes(x = lon, y = lat, r = r, group = grp),
                    cols = paste0("C", 1:K), colour = "#111111", linewidth = 0.25) +
    scale_fill_manual(values = ccols, guide = "none")

  save_figure(p_avg_density_small, sprintf("average_map_structure_K%d_density_small", K), dir_combined)
}

cat("=== All Cluster Average Maps (Standalone & Combined Density K=2..10) Complete! ===\n")
