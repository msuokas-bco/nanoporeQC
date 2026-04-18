# tests/testthat/test-theme-helpers.R
# Tests for theme helper functions (no minimap2/samtools needed)

library(testthat)

# ============================================================================
# Color Utility Tests
# ============================================================================

test_that("hex_to_rgb parses hex colors correctly", {
  expect_equal(nanoporeQC:::hex_to_rgb("#FF0000"), c(255, 0, 0))
  expect_equal(nanoporeQC:::hex_to_rgb("#00FF00"), c(0, 255, 0))
  expect_equal(nanoporeQC:::hex_to_rgb("#0000FF"), c(0, 0, 255))
  expect_equal(nanoporeQC:::hex_to_rgb("#FFFFFF"), c(255, 255, 255))
  expect_equal(nanoporeQC:::hex_to_rgb("#000000"), c(0, 0, 0))
})

test_that("hex_to_rgb works with and without leading hash", {
  expect_equal(nanoporeQC:::hex_to_rgb("FF0000"), c(255, 0, 0))
  expect_equal(nanoporeQC:::hex_to_rgb("#FF0000"), c(255, 0, 0))
})

test_that("rgb_to_hex formats correctly", {
  expect_equal(nanoporeQC:::rgb_to_hex(255, 0, 0), "#FF0000")
  expect_equal(nanoporeQC:::rgb_to_hex(0, 255, 0), "#00FF00")
  expect_equal(nanoporeQC:::rgb_to_hex(0, 0, 255), "#0000FF")
})

test_that("lighten_color moves color towards white", {
  result <- nanoporeQC:::lighten_color("#000000", 0.5)
  expect_equal(result, "#808080")

  result <- nanoporeQC:::lighten_color("#FF0000", 1.0)
  expect_equal(result, "#FFFFFF")
})

test_that("darken_color moves color towards black", {
  result <- nanoporeQC:::darken_color("#FFFFFF", 0.5)
  expect_equal(result, "#808080")

  result <- nanoporeQC:::darken_color("#FF0000", 1.0)
  expect_equal(result, "#000000")
})

test_that("get_contrast_color returns appropriate text color", {
  # Dark backgrounds should get white text
  expect_equal(nanoporeQC:::get_contrast_color("#000000"), "#FFFFFF")
  expect_equal(nanoporeQC:::get_contrast_color("#1a1a1a"), "#FFFFFF")

  # Light backgrounds should get dark text
  expect_equal(nanoporeQC:::get_contrast_color("#FFFFFF"), "#333333")
  expect_equal(nanoporeQC:::get_contrast_color("#F0F0F0"), "#333333")
})

# ============================================================================
# Theme Generation Tests
# ============================================================================

test_that("generate_theme_config returns all required properties", {
  config <- nanoporeQC:::generate_theme_config("#DFA878")
  required_props <- c("header_bg", "header_color", "stripe_color", "hover_color",
                      "font_family", "title_color", "value_color", "fill_color",
                      "mean_line_color", "percentile_line_color", "bg_color",
                      "text_color")
  expect_true(all(required_props %in% names(config)))
})

test_that("generate_theme_config uses primary color appropriately", {
  config <- nanoporeQC:::generate_theme_config("#2780e3")
  expect_equal(config$header_bg, "#2780e3")
  expect_equal(config$fill_color, "#2780e3")
  expect_equal(config$value_color, "#2780e3")
})

test_that("generate_theme_config derives secondary color when not provided", {
  config <- nanoporeQC:::generate_theme_config("#2780e3")
  # Secondary should be a lighter variant of primary
  expect_false(config$percentile_line_color == "#2780e3")
  # Should still be a valid hex
  expect_match(config$percentile_line_color, "^#[0-9A-F]{6}$")
})

test_that("generate_theme_config uses provided secondary color", {
  config <- nanoporeQC:::generate_theme_config("#2780e3", secondary_color = "#7BADEE")
  expect_equal(config$percentile_line_color, "#7BADEE")
})

test_that("generate_theme_config font_family handles presets", {
  config_sans <- nanoporeQC:::generate_theme_config("#2780e3", font_family = "sans")
  expect_match(config_sans$font_family, "Helvetica")

  config_serif <- nanoporeQC:::generate_theme_config("#2780e3", font_family = "serif")
  expect_match(config_serif$font_family, "Georgia")

  config_ubuntu <- nanoporeQC:::generate_theme_config("#2780e3", font_family = "ubuntu")
  expect_match(config_ubuntu$font_family, "Ubuntu")
})

test_that("generate_theme_config font_family passes through custom strings", {
  custom_font <- "'My Custom Font', sans-serif"
  config <- nanoporeQC:::generate_theme_config("#2780e3", font_family = custom_font)
  expect_equal(config$font_family, custom_font)
})

# ============================================================================
# THEME_SPECS Tests
# ============================================================================

test_that("THEME_SPECS contains all 15 expected themes", {
  expected <- c(
    "sandstone", "simplex", "flatly", "journal", "lumen", "spacelab", "united",
    "default", "cerulean", "cosmo", "litera", "lux", "materia", "minty", "zephyr"
  )
  expect_equal(sort(names(nanoporeQC:::THEME_SPECS)), sort(expected))
})

test_that("each THEME_SPECS entry has required fields", {
  for (theme_name in names(nanoporeQC:::THEME_SPECS)) {
    spec <- nanoporeQC:::THEME_SPECS[[theme_name]]
    expect_true(all(c("primary", "secondary", "bg", "font") %in% names(spec)),
                info = paste("Theme:", theme_name))
    expect_match(spec$primary, "^#[0-9a-fA-F]{6}$",
                 info = paste("Theme:", theme_name))
    expect_true(spec$bg %in% c("light", "white", "subtle"),
                info = paste("Theme:", theme_name))
  }
})

test_that("get_theme_config works for all themes", {
  for (theme_name in names(nanoporeQC:::THEME_SPECS)) {
    config <- nanoporeQC:::get_theme_config(theme_name)
    expect_type(config, "list")
    expect_true("header_bg" %in% names(config),
                info = paste("Theme:", theme_name))
  }
})

test_that("get_theme_config errors on unknown theme", {
  expect_error(nanoporeQC:::get_theme_config("nonexistent"), "Unknown theme")
})
