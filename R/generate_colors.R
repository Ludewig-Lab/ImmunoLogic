#' Generate Automatic Color Palettes for Plots
#'
#' @param x An integer (number of colors), a character/factor vector (categories),
#'   or a data.frame (uses the unique values of the first column).
#' @param palette_type Character, either `"default"` (uses heuristic grouping palettes)
#'   or `"condition"` (uses distinct, high-contrast condition colors).
#'
#' @return A character vector of hex color codes. If `x` contains names/categories,
#'   the returned vector will be named accordingly.
#' @export
generate_colors <- function(x, palette_type = "default") {
  palette_type <- match.arg(palette_type, c("default", "condition"))

  # 1. Parse Input to determine categories and 'n'
  if (is.data.frame(x)) {
    cats <- unique(as.character(x[[1]]))
    n <- length(cats)
  } else if (is.numeric(x) && length(x) == 1) {
    cats <- NULL
    n <- x
  } else {
    cats <- unique(as.character(x))
    n <- length(cats)
  }

  if (n == 0) return(character(0))

  # 2. Define Palettes
  disease_colors <- c("#dfc27d", "#BE3144", "#202547", "#355C7D", "#779d8d")
  fibroblast_colors <- c("#D53E4F", "#f4a582", "#ff7b7b", "#8e0b00", "#FEE08B",
                         "#42090D", "#FF7B00", "#FFF4DF")
  nice_colors <- c("#67001f", "#D53E4F", "#f4a582", "#FEE08B", "#003c30", "#01665e",
                   "#66C2A5", "#3288BD", "#BEAED4", "#c7eae5", "#355C7D", "#202547",
                   "#B45B5C", "#8c510a")
  extended_colors <- c("#fde0dd", "#fa9fb5", "#d95f0e", "#dd1c77", "#D53E4F",
                       "#f4a582", "#FEE08B", "#f03b20", "#ffffcc", "#43a2ca",
                       "#1c9099", "#355C7D", "#3288BD", "#BEAED4", "#756bb1",
                       "#c7eae5")
  large_palette <- c(
    "#fde0dd", "#fa9fb5", "#f768a1", "#dd1c77", "#980043",
    "#f4a582", "#fdae61", "#f46d43", "#d73027", "#a50026",
    "#fee08b", "#ffffbf", "#e6f598", "#99d594", "#66c2a5",
    "#43a2ca", "#1c9099", "#016c59", "#3288bd", "#5e4fa2",
    "#beaed4", "#9e9ac8", "#756bb1", "#542788", "#3f007d",
    "#c7eae5", "#a6bddb", "#74a9cf", "#3690c0", "#045a8d"
  )
  cond_palette <- c(
    "steelblue", "orange", "purple", "forestgreen", "firebrick",
    "goldenrod", "turquoise", "violet", "darkolivegreen", "coral",
    "slateblue", "tomato", "mediumorchid", "darkgoldenrod", "cadetblue",
    "deeppink", "darkseagreen", "dodgerblue", "sienna", "darkcyan",
    "rosybrown", "lightblue", "limegreen", "maroon", "peru"
  )

  # 3. Assign Colors based on Type and Heuristics
  out_colors <- NULL

  if (palette_type == "condition") {
    if (n <= length(cond_palette)) {
      out_colors <- cond_palette[1:n]
    } else {
      out_colors <- c(cond_palette, grDevices::rainbow(n - length(cond_palette)))
    }
  } else {
    # Default Group Logic
    is_disease <- !is.null(cats) && any(grepl("healthy|explant|visit", cats, ignore.case = TRUE))
    is_fibro <- !is.null(cats) && any(grepl("Fb|Periv|VSMC", cats, ignore.case = TRUE))

    if (n <= 5 && is_disease) {
      out_colors <- disease_colors[1:n]
    } else if (n <= 8 && is_fibro) {
      out_colors <- fibroblast_colors[1:n]
    } else if (n <= 14) {
      out_colors <- nice_colors[1:n]
    } else if (n <= 16) {
      out_colors <- extended_colors[1:n]
    } else if (n <= 30) {
      out_colors <- large_palette[1:n]
    } else {
      out_colors <- c(large_palette, grDevices::rainbow(n - length(large_palette)))
    }
  }

  # 4. Attach names if categorical input was provided
  if (!is.null(cats)) {
    names(out_colors) <- cats
  }

  return(out_colors)
}
