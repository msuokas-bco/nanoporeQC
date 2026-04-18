# Theme helper functions for generating theme configurations from minimal inputs

#' Convert hex color to RGB
#' @keywords internal
hex_to_rgb <- function(hex) {
  hex <- sub("^#", "", hex)
  r <- strtoi(substr(hex, 1, 2), 16)
  g <- strtoi(substr(hex, 3, 4), 16)
  b <- strtoi(substr(hex, 5, 6), 16)
  c(r, g, b)
}

#' Convert RGB to hex color
#' @keywords internal
rgb_to_hex <- function(r, g, b) {
  sprintf("#%02X%02X%02X", round(r), round(g), round(b))
}

#' Lighten a color by a factor
#' @keywords internal
lighten_color <- function(hex, factor = 0.2) {
  rgb <- hex_to_rgb(hex)
  # Move towards white (255, 255, 255)
  rgb <- rgb + (255 - rgb) * factor
  rgb_to_hex(rgb[1], rgb[2], rgb[3])
}

#' Darken a color by a factor
#' @keywords internal
darken_color <- function(hex, factor = 0.2) {
  rgb <- hex_to_rgb(hex)
  # Move towards black (0, 0, 0)
  rgb <- rgb * (1 - factor)
  rgb_to_hex(rgb[1], rgb[2], rgb[3])
}

#' Get contrasting text color (dark or light) for a background
#' @keywords internal
get_contrast_color <- function(bg_hex) {
  rgb <- hex_to_rgb(bg_hex)
  # Calculate luminance using WCAG formula
  luminance <- (0.299 * rgb[1] + 0.587 * rgb[2] + 0.114 * rgb[3]) / 255
  if (luminance > 0.5) "#333333" else "#FFFFFF"
}

#' Generate a complete theme configuration from minimal inputs
#' @keywords internal
generate_theme_config <- function(primary_color, secondary_color = NULL,
                                   background = "light", font_family = "sans") {
  # Determine secondary color if not provided
  if (is.null(secondary_color)) {
    secondary_color <- lighten_color(primary_color, 0.2)
  }

  # Determine background color
  bg_color <- switch(background,
                     "light" = lighten_color(primary_color, 0.85),
                     "white" = "#FFFFFF",
                     "subtle" = lighten_color(primary_color, 0.9),
                     "#FFFFFF")

  # Determine header colors with contrast
  header_color <- get_contrast_color(primary_color)

  # Determine font family CSS stack
  font_stack <- switch(font_family,
                       "sans" = "'Helvetica Neue', Helvetica, Arial, sans-serif",
                       "serif" = "'Georgia', 'Times New Roman', Times, serif",
                       "ubuntu" = "'Ubuntu', 'Helvetica Neue', Helvetica, Arial, sans-serif",
                       font_family)

  list(
    header_bg = primary_color,
    header_color = header_color,
    stripe_color = lighten_color(primary_color, 0.80),
    hover_color = lighten_color(primary_color, 0.60),
    font_family = font_stack,
    title_color = darken_color(primary_color, 0.15),
    value_color = primary_color,
    fill_color = primary_color,
    mean_line_color = darken_color(primary_color, 0.30),
    percentile_line_color = secondary_color,
    bg_color = bg_color,
    text_color = "#333333"
  )
}

# Theme specifications: minimal inputs needed to generate complete theme configs.
# Primary colors match Quarto/Bootswatch themes. Secondary colors are harmonious
# lighter variants of the primary, used for the percentile line in plots.
# When secondary is NULL, it is auto-derived as a lightened primary.
THEME_SPECS <- list(
  # Existing themes (preserve original visual identity)
  sandstone = list(primary = "#DFA878", secondary = "#C4956C", bg = "light", font = "sans"),
  simplex   = list(primary = "#D9230F", secondary = "#FF6B5E", bg = "white", font = "sans"),
  flatly    = list(primary = "#2C3E50", secondary = "#18BC9C", bg = "light", font = "sans"),
  journal   = list(primary = "#EB6864", secondary = "#F7CACA", bg = "white", font = "serif"),
  lumen     = list(primary = "#4D79B6", secondary = "#A3BEE8", bg = "white", font = "sans"),
  spacelab  = list(primary = "#446E9B", secondary = "#8FAAC7", bg = "light", font = "sans"),
  united    = list(primary = "#FF6F21", secondary = "#FFB289", bg = "white", font = "ubuntu"),

  # New light themes from Quarto/Bootswatch (primary colors from Bootswatch SCSS).
  # Secondary colors use harmonious lighter shades rather than Bootswatch's
  # "success" (often green) to avoid visual clashing in plots.
  default   = list(primary = "#0d6efd", secondary = "#6ea8fe", bg = "light", font = "sans"),
  cerulean  = list(primary = "#2fa4e7", secondary = "#7DC7F2", bg = "light", font = "sans"),
  cosmo     = list(primary = "#2780e3", secondary = "#7BADEE", bg = "light", font = "sans"),
  litera    = list(primary = "#4582ec", secondary = "#89ADF4", bg = "light", font = "serif"),
  lux       = list(primary = "#1a1a1a", secondary = "#8C8C8C", bg = "subtle", font = "sans"),
  materia   = list(primary = "#2196f3", secondary = "#64B5F6", bg = "light", font = "sans"),
  minty     = list(primary = "#78c2ad", secondary = "#A6D8C7", bg = "light", font = "sans"),
  zephyr    = list(primary = "#3459e6", secondary = "#7D95EE", bg = "light", font = "sans")
)

#' Generate theme config by name
#' @keywords internal
get_theme_config <- function(theme) {
  if (!theme %in% names(THEME_SPECS)) {
    stop("Unknown theme: ", theme)
  }

  spec <- THEME_SPECS[[theme]]
  generate_theme_config(
    primary_color = spec$primary,
    secondary_color = spec$secondary,
    background = spec$bg,
    font_family = spec$font
  )
}
