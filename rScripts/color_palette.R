# Central Color Palette for Oak Wilt (Bretziella fagacearum) GBS Analysis
#
# Keyed by the `Clusters` column in metadata/metadata.csv: 9 US groups (the six
# eastern states MD/NY/OH/PA/SC/WV are merged into "Eastern"). Sourced by the
# analysis scripts, which call get_pop_colors()/get_pop_shapes() on
# popNames(gl_bret). Colours are the 9-colour Okabe-Ito colourblind-safe palette;
# shapes are provided as a secondary channel.

# Okabe-Ito colour per Clusters group (colourblind-safe)
.base_colors <- c(
  "AR"      = "#E69F00",  # orange
  "Eastern" = "#56B4E9",  # sky blue
  "IA"      = "#009E73",  # bluish green
  "IL"      = "#F0E442",  # yellow
  "MI"      = "#0072B2",  # blue
  "MN"      = "#D55E00",  # vermillion
  "MO"      = "#CC79A7",  # reddish purple
  "TX"      = "#000000",  # black
  "WI"      = "#999999"   # grey
)

# Distinct point shape per Clusters group (secondary channel to colour)
.base_shapes <- c(
  "AR"      = 16,  # filled circle
  "Eastern" = 10,  # circle-plus
  "IA"      = 17,  # filled triangle
  "IL"      = 15,  # filled square
  "MI"      = 8,   # asterisk
  "MN"      = 3,   # plus
  "MO"      = 4,   # cross
  "TX"      = 6,   # open inverted triangle
  "WI"      = 7    # square-cross
)

# Okabe-Ito colour per host `Sub-section` (UA = unknown host)
host_colors <- c(
  "Red oak"   = "#D55E00",  # vermillion
  "White oak" = "#0072B2",  # blue
  "Chestnut"  = "#009E73",  # bluish green
  "UA"        = "#999999"   # grey
)

# Full names for pretty legends: get_pop_labels(popNames(gl_bret))
.state_names <- c(
  "AR" = "Arkansas",   "Eastern" = "Eastern", "IA" = "Iowa",
  "IL" = "Illinois",   "MI" = "Michigan",     "MN" = "Minnesota",
  "MO" = "Missouri",   "TX" = "Texas",        "WI" = "Wisconsin"
)

# Function to get colors for a set of names
get_pop_colors <- function(names) {
  missing <- setdiff(names, names(.base_colors))
  if (length(missing) > 0) {
    fallback_cols <- c("#332288","#88CCEE","#44AA99","#117733","#999933",
                       "#DDCC77","#CC6677","#882255","#AA4499","#D55E00")
    new_cols <- colorRampPalette(fallback_cols)(length(missing))
    names(new_cols) <- missing
    return(c(.base_colors, new_cols)[names])
  }
  return(.base_colors[names])
}

# Function to get shapes for a set of names
get_pop_shapes <- function(names) {
  missing <- setdiff(names, names(.base_shapes))
  if (length(missing) > 0) {
    new_shapes <- rep(c(0, 1, 2, 5, 6), length.out = length(missing))
    names(new_shapes) <- missing
    return(c(.base_shapes, new_shapes)[names])
  }
  return(.base_shapes[names])
}

# Function to get pretty (full-name) labels; falls back to the code itself
get_pop_labels <- function(names) {
  out <- .state_names[names]
  out[is.na(out)] <- names[is.na(out)]
  setNames(out, names)
}
