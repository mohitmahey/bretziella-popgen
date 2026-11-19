rm(list = ls())

################################################################################
## Oak Wilt (Bretziella fagacearum) US County Spread & Density Map            ##
##                                                                            ##
## Reconstructs the EDDMapS county density map in publication-quality R       ##
## focusing strictly on the Continental United States (USA only).            ##
## Data source: oak_wilt_spread_data/62795/mappings.csv                       ##
## Outputs: Publication-grade PNG (400 DPI), vector PDF (Cairo), and SVG.     ##
################################################################################

suppressMessages({
  library(tidyverse)
  library(maps)
  library(ggnewscale)
})

# Cairo graphics device for true sRGB colour rendering
options(bitmapType = "cairo")

# Source central project color palette if available
if (file.exists("rScripts/color_palette.R")) {
  source("rScripts/color_palette.R")
}

# ── Paths & Output Directories ────────────────────────────────────────────────
data_file <- "oak_wilt_spread_data/62795/mappings.csv"
meta_file <- "metadata/metadata.csv"
out_dir   <- "figures/oak_wilt_spread"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

png_dir <- file.path(out_dir, "png")
pdf_dir <- file.path(out_dir, "pdf")
svg_dir <- file.path(out_dir, "svg")
dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(svg_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== 1. Loading and Parsing EDDMapS Occurrence Data ===\n")
stopifnot(file.exists(data_file))

raw_df <- read.csv(data_file, stringsAsFactors = FALSE, check.names = FALSE)
cat("Loaded", nrow(raw_df), "total records from", data_file, "\n")

# Parse Location ("County, State, Country")
loc_clean <- trimws(gsub("\"", "", raw_df$Location))
parts     <- strsplit(loc_clean, ", ")

state_raw   <- sapply(parts, function(x) if (length(x) >= 2) tolower(trimws(x[2])) else NA)
county_raw  <- sapply(parts, function(x) if (length(x) >= 1) tolower(trimws(x[1])) else NA)
country_raw <- sapply(parts, function(x) if (length(x) >= 3) tolower(trimws(x[3])) else NA)

# Filter strictly to United States positive / detected records
us_mask <- which(country_raw == "united states" & raw_df$Status != "Negative")
us_data <- data.frame(
  state_raw  = state_raw[us_mask],
  county_raw = county_raw[us_mask],
  stringsAsFactors = FALSE
)
cat("Filtered to", nrow(us_data), "confirmed US records.\n")

# ── County Normalization for Seamless Join with map_data('county') ─────────────
normalize_county <- function(s, c) {
  c <- gsub("[.]", "", c)               # Remove periods (e.g. st. -> st)
  c <- gsub(" county$", "", c)          # Strip trailing 'county'
  c <- gsub(" parish$", "", c)          # Strip trailing 'parish'

  # Standard county name variants in maps package
  c[c == "dupage"]      <- "du page"
  c[c == "dekalb"]       <- "de kalb"
  c[c == "laporte"]      <- "la porte"
  c[c == "dewitt"]       <- "de witt"
  c[c == "st clair"]     <- "st clair"
  c[c == "st joseph"]    <- "st joseph"
  c[c == "st charles"]   <- "st charles"
  c[c == "st francois"]  <- "st francois"
  c[c == "st louis"]     <- "st louis"

  # Virginia independent cities -> mapped to encompassing counties
  c <- gsub(" [(]independent city[)]", "", c)
  c[s == "virginia" & c == "bristol"]                       <- "washington"
  c[s == "virginia" & c %in% c("buena vista", "lexington")] <- "rockbridge"
  c[s == "virginia" & c == "covington"]                     <- "alleghany"
  c[s == "virginia" & c == "harrisonburg"]                  <- "rockingham"
  c[s == "virginia" & c == "norton"]                        <- "wise"
  c[s == "virginia" & c %in% c("staunton", "waynesboro")]   <- "augusta"
  c[s == "virginia" & c == "winchester"]                    <- "frederick"

  # District of Columbia
  s[s %in% c("d. columbia", "d columbia")] <- "district of columbia"
  c[s == "district of columbia"]           <- "washington"

  data.frame(state = s, county = c, stringsAsFactors = FALSE)
}

norm_us <- normalize_county(us_data$state_raw, us_data$county_raw)

# County aggregation: count reports per county
county_counts <- norm_us %>%
  count(state, county, name = "n_records")

cat("Aggregated into", nrow(county_counts), "unique reporting US counties.\n")

# ── Basemap Layers & Joining ──────────────────────────────────────────────────
cat("=== 2. Preparing US Basemap Layers ===\n")
us_states   <- map_data("state")
us_counties <- map_data("county")
lakes       <- map_data("lakes") %>% filter(long >= -125, long <= -65, lat >= 24, lat <= 52)

# Join counts to all US counties
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

# ── Exact EDDMapS Colour Palette ──────────────────────────────────────────────
density_palette <- c(
  "No Data" = "#FFFFFF",  # Unreported / white
  "1"       = "#CEE8C3",  # Lightest sage green
  "2-10"    = "#87C27E",  # Soft spring green
  "11-100"  = "#45894C",  # Forest green
  "101-500" = "#2E6B34",  # Deep pine green
  "501+"    = "#1A4320"   # Hyper-dense dark evergreen
)

# ── Scale bar helper for Full Continental USA (1000 km) ──────────────────────
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

# ══════════════════════════════════════════════════════════════════════════════
# MAP 1: Standalone EDDMapS US County Density Map (USA Only)
# ══════════════════════════════════════════════════════════════════════════════
cat("=== 3. Rendering MAP 1: Standalone EDDMapS US County Density Map ===\n")

p1 <- ggplot() +
  # 1. County polygons with density choropleth fill
  geom_polygon(data = us_counties_density,
               aes(x = long, y = lat, group = group, fill = density_bin),
               colour = "#D1D5DB", linewidth = 0.12) +
  # 2. State boundary outlines
  geom_polygon(data = us_states,
               aes(x = long, y = lat, group = group),
               fill = NA, colour = "#374151", linewidth = 0.45) +
  # 3. Great Lakes water bodies
  geom_polygon(data = lakes,
               aes(x = long, y = lat, group = group),
               fill = "#EAF2F8", colour = "#90A4AE", linewidth = 0.25) +
  # 4. Color scale & Legend
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
  # 5. Fixed Continental USA projection
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

# Export MAP 1
fname1 <- "map_oak_wilt_eddmaps_density_USA"
ggsave(file.path(png_dir, paste0(fname1, ".png")), p1, width = 12, height = 7, dpi = 400, bg = "white")
ggsave(file.path(pdf_dir, paste0(fname1, ".pdf")), p1, width = 12, height = 7, device = cairo_pdf, bg = "white")
saveRDS(p1, file.path(out_dir, paste0(fname1, ".rds")))
cat("Saved:", fname1, "in PNG, PDF, and RDS.\n")

# ══════════════════════════════════════════════════════════════════════════════
# MAP 2: Overlay with PopGen Sample Collection Sites (Sized Obsidian Dots)
# ══════════════════════════════════════════════════════════════════════════════
if (file.exists(meta_file)) {
  cat("=== 4. Rendering MAP 2: Spread Map Overlaid with PopGen Isolates ===\n")
  meta <- read.csv(meta_file, check.names = FALSE, stringsAsFactors = FALSE)
  names(meta) <- trimws(sub("^\uFEFF", "", names(meta)))
  meta <- meta %>% mutate(across(where(is.character), trimws))

  # Aggregated county points for popgen samples
  sample_pts <- meta %>%
    filter(!is.na(Longitude) & !is.na(Latitude)) %>%
    group_by(Clusters, State, County) %>%
    summarise(
      Longitude = mean(as.numeric(Longitude)),
      Latitude  = mean(as.numeric(Latitude)),
      n         = dplyr::n(),
      .groups   = "drop"
    )

  p2 <- p1 +
    # White halo backing for high contrast over dark greens
    geom_point(data = sample_pts,
               aes(x = Longitude, y = Latitude, size = n),
               shape = 16, colour = "#FFFFFF", alpha = 0.95) +
    # Obsidian core with crisp white border
    geom_point(data = sample_pts,
               aes(x = Longitude, y = Latitude, size = n),
               shape = 21, fill = "#111111", colour = "#FFFFFF", stroke = 0.45) +
    scale_size_continuous(
      range  = c(1.4, 3.4),
      breaks = c(1, 2, 3, 4),
      name   = "PopGen\nSamples",
      guide  = guide_legend(
        override.aes = list(shape = 21, fill = "#111111", colour = "#FFFFFF", stroke = 0.5),
        keywidth  = unit(1.1, "cm"),
        keyheight = unit(0.55, "cm"),
        order = 2
      )
    ) +
    theme(
      legend.position  = c(0.90, 0.27),
      legend.box       = "vertical",
      legend.box.just  = "left",
      legend.spacing.y = unit(0.25, "cm")
    )

  fname2 <- "map_oak_wilt_density_with_popgen_samples"
  ggsave(file.path(png_dir, paste0(fname2, ".png")), p2, width = 12, height = 7, dpi = 400, bg = "white")
  ggsave(file.path(pdf_dir, paste0(fname2, ".pdf")), p2, width = 12, height = 7, device = cairo_pdf, bg = "white")
  saveRDS(p2, file.path(out_dir, paste0(fname2, ".rds")))
  cat("Saved:", fname2, "in PNG, PDF, and RDS.\n")
}

# ── Convert PDFs to Vector SVG via pdftocairo (if installed) ───────────────────
if (nzchar(Sys.which("pdftocairo"))) {
  cat("Converting PDF maps to vector SVG via pdftocairo ...\n")
  all_pdfs <- list.files(pdf_dir, pattern = "\\.pdf$", full.names = TRUE)
  for (pdf_path in all_pdfs) {
    svg_path <- file.path(svg_dir, sub("\\.pdf$", ".svg", basename(pdf_path)))
    system2("pdftocairo", c("-svg", shQuote(pdf_path), shQuote(svg_path)))
  }
  cat("SVG conversion complete.\n")
}

cat("=== All Done! Maps saved to:", out_dir, "===\n")
