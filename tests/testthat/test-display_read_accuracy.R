# tests/testthat/test-display_read_accuracy.R

library(testthat)
library(ggplot2)

# Helper function to create a minimal test FASTQ file
create_test_fastq <- function(n_reads = 10) {
  # DNA CS standard sequence (first 500bp for testing)
  ref_seq <- "GCCATCAGATTGTGTTTGTTAGTCGCTTTTTTTTTTTGGAATTTTTTTTTTGGAATTTTTTTTTTGCGCTAACAACCTCCTGCCGTTTTGCCCGTGCATATCGGTCACGAACAAATCTGATTACTAAACACAGTAGCCTGGATTTGTTCTATCAGTAATCGACCTTATTCCTAATTAAATAGAGCAAATCCCCTTATTGGGGGTAAGACATGAAGATGCCAGAAAAACATGACCTGTTGGCCGCCATTCTCGCGGCAAAGGAACAAGGCATCGGGGCAATCCTTGCGTTTGCAATGGCGTACCTTCGCGGCAGATATAATGGCGGTGCGTTTACAAAAACAGTAATCGACGCAACGATGTGCGCCATTATCGCCTAGTTCATTCGTGACCTTCTCGACTTCGCCGGACTAAGTAGCAATCTCGCTTATATAACGAGCGTGTTTATCGGCTACATCGGTACTGACTCGATTGGTTCGCTTATCAAACGCTTCGCTGCTAAAAAAGCCGGAGTAGAAGATGGTAGAAATCAATAATCAACGTAAGGCGTTCCTCGA"

  temp_fastq <- tempfile(fileext = ".fastq")
  fastq_lines <- character()

  for (i in 1:n_reads) {
    # Create reads with slight variations (simulate sequencing errors)
    read_seq <- ref_seq
    # Introduce 1-2 random substitutions
    if (runif(1) > 0.3) {
      pos <- sample(1:nchar(read_seq), 1)
      substr(read_seq, pos, pos) <- sample(c("A", "T", "G", "C"), 1)
    }

    qual_string <- paste(rep("I", nchar(read_seq)), collapse = "")

    fastq_lines <- c(fastq_lines,
                     sprintf("@read_%d", i),
                     read_seq,
                     "+",
                     qual_string)
  }

  writeLines(fastq_lines, temp_fastq)
  return(temp_fastq)
}

# Check if system dependencies are available
has_minimap2 <- system("which minimap2", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
has_samtools <- system("which samtools", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
has_dependencies <- has_minimap2 && has_samtools

# ============================================================================
# Input Validation Tests
# ============================================================================

test_that("display_read_accuracy validates FASTQ input", {
  expect_error(
    display_read_accuracy(NULL),
    "FASTQ file must be a single file path string"
  )

  expect_error(
    display_read_accuracy(c("file1.fastq", "file2.fastq")),
    "FASTQ file must be a single file path string"
  )

  expect_error(
    display_read_accuracy(123),
    "FASTQ file must be a single file path string"
  )
})

test_that("display_read_accuracy checks file existence", {
  expect_error(
    display_read_accuracy("nonexistent_file.fastq"),
    "FASTQ file does not exist"
  )
})

test_that("display_read_accuracy validates output_format parameter", {
  skip_if_not(has_dependencies, "minimap2 or samtools not available")

  temp_fastq <- create_test_fastq(5)
  on.exit(unlink(temp_fastq))

  expect_error(
    display_read_accuracy(temp_fastq, output_format = "invalid"),
    "'arg' should be one of"
  )

  expect_error(
    display_read_accuracy(temp_fastq, output_format = "PDF"),
    "'arg' should be one of"
  )
})

test_that("display_read_accuracy validates table_format parameter", {
  skip_if_not(has_dependencies, "minimap2 or samtools not available")

  temp_fastq <- create_test_fastq(5)
  on.exit(unlink(temp_fastq))

  expect_error(
    display_read_accuracy(temp_fastq, table_format = "invalid"),
    "'arg' should be one of"
  )
})

test_that("display_read_accuracy validates theme parameter", {
  skip_if_not(has_dependencies, "minimap2 or samtools not available")

  temp_fastq <- create_test_fastq(5)
  on.exit(unlink(temp_fastq))

  expect_error(
    display_read_accuracy(temp_fastq, theme = "invalid"),
    "'arg' should be one of"
  )
})

# ============================================================================
# System Dependencies Tests
# ============================================================================

test_that("display_read_accuracy checks for minimap2", {
  skip_if(has_minimap2, "minimap2 is available")

  temp_fastq <- create_test_fastq(5)
  on.exit(unlink(temp_fastq))

  expect_error(
    display_read_accuracy(temp_fastq),
    "minimap2 is not installed or not in PATH"
  )
})

test_that("display_read_accuracy checks for samtools", {
  skip_if(has_samtools, "samtools is available")
  skip_if_not(has_minimap2, "minimap2 not available")

  temp_fastq <- create_test_fastq(5)
  on.exit(unlink(temp_fastq))

  expect_error(
    display_read_accuracy(temp_fastq),
    "samtools is not installed or not in PATH"
  )
})

# ============================================================================
# Output Format Tests
# ============================================================================

test_that("display_read_accuracy returns ggplot object with output_format='plot'", {
  skip_if_not(has_dependencies, "minimap2 or samtools not available")

  temp_fastq <- create_test_fastq(10)
  on.exit(unlink(temp_fastq))

  result <- display_read_accuracy(temp_fastq, output_format = "plot")

  expect_s3_class(result, "ggplot")
  expect_true("accuracy" %in% names(result$data))
})

test_that("display_read_accuracy returns table with output_format='table'", {
  skip_if_not(has_dependencies, "minimap2 or samtools not available")
  skip_if_not_installed("htmltools")

  temp_fastq <- create_test_fastq(10)
  on.exit(unlink(temp_fastq))

  result <- display_read_accuracy(temp_fastq, output_format = "table", table_format = "kable")

  # htmltools::browsable returns an object with class "shiny.tag.list" or similar
  # Check that it's an HTML object that can be displayed
  expect_true(inherits(result, "html") || inherits(result, "shiny.tag.list"))
  expect_true("html" %in% class(result) || "shiny.tag" %in% class(result))
})

test_that("display_read_accuracy returns data.frame with table_format='data.frame'", {
  skip_if_not(has_dependencies, "minimap2 or samtools not available")

  temp_fastq <- create_test_fastq(10)
  on.exit(unlink(temp_fastq))

  result <- display_read_accuracy(temp_fastq, output_format = "table", table_format = "data.frame")

  expect_s3_class(result, "data.frame")
  expect_true(all(c("Metric", "Value") %in% names(result)))
  expect_true("Mean Accuracy" %in% result$Metric)
  expect_true("Median Accuracy" %in% result$Metric)
})

test_that("display_read_accuracy returns list with output_format='both'", {
  skip_if_not(has_dependencies, "minimap2 or samtools not available")

  temp_fastq <- create_test_fastq(10)
  on.exit(unlink(temp_fastq))

  result <- display_read_accuracy(temp_fastq, output_format = "both", table_format = "data.frame")

  expect_type(result, "list")
  expect_named(result, c("plot", "table"))
  expect_s3_class(result$plot, "ggplot")
  expect_s3_class(result$table, "data.frame")
})

# ============================================================================
# Theme Tests
# ============================================================================

test_that("display_read_accuracy works with all themes", {
  skip_if_not(has_dependencies, "minimap2 or samtools not available")

  temp_fastq <- create_test_fastq(10)
  on.exit(unlink(temp_fastq))

  themes <- c("sandstone", "simplex", "flatly")

  for (theme_name in themes) {
    result <- display_read_accuracy(temp_fastq, theme = theme_name, output_format = "plot")
    expect_s3_class(result, "ggplot")
  }
})

# ============================================================================
# File Output Tests
# ============================================================================

test_that("display_read_accuracy creates output directory if needed", {
  skip_if_not(has_dependencies, "minimap2 or samtools not available")

  temp_fastq <- create_test_fastq(10)
  temp_dir <- tempfile()
  on.exit({
    unlink(temp_fastq)
    unlink(temp_dir, recursive = TRUE)
  })

  expect_false(dir.exists(temp_dir))

  result <- display_read_accuracy(temp_fastq, output_dir = temp_dir, output_format = "plot")

  expect_true(dir.exists(temp_dir))
})

test_that("display_read_accuracy saves plot file when output_dir is specified", {
  skip_if_not(has_dependencies, "minimap2 or samtools not available")

  temp_fastq <- create_test_fastq(10)
  temp_dir <- tempfile()
  on.exit({
    unlink(temp_fastq)
    unlink(temp_dir, recursive = TRUE)
  })

  result <- display_read_accuracy(temp_fastq, output_dir = temp_dir, output_format = "plot")

  expect_true(file.exists(file.path(temp_dir, "accuracy_distribution.png")))
})

test_that("display_read_accuracy saves table file when output_dir is specified", {
  skip_if_not(has_dependencies, "minimap2 or samtools not available")

  temp_fastq <- create_test_fastq(10)
  temp_dir <- tempfile()
  on.exit({
    unlink(temp_fastq)
    unlink(temp_dir, recursive = TRUE)
  })

  result <- display_read_accuracy(temp_fastq, output_dir = temp_dir,
                                  output_format = "table", table_format = "data.frame")

  expect_true(file.exists(file.path(temp_dir, "accuracy_summary.csv")))
})

test_that("display_read_accuracy saves raw accuracies when output_dir is specified", {
  skip_if_not(has_dependencies, "minimap2 or samtools not available")

  temp_fastq <- create_test_fastq(10)
  temp_dir <- tempfile()
  on.exit({
    unlink(temp_fastq)
    unlink(temp_dir, recursive = TRUE)
  })

  result <- display_read_accuracy(temp_fastq, output_dir = temp_dir, output_format = "plot")

  expect_true(file.exists(file.path(temp_dir, "read_accuracies.txt")))

  # Check that file contains numeric values
  accuracies <- read.table(file.path(temp_dir, "read_accuracies.txt"))
  expect_true(all(accuracies[,1] >= 0 & accuracies[,1] <= 1))
})

test_that("display_read_accuracy handles output_dir with and without trailing slash", {
  skip_if_not(has_dependencies, "minimap2 or samtools not available")

  temp_fastq <- create_test_fastq(10)
  temp_dir1 <- tempfile()
  temp_dir2 <- paste0(tempfile(), "/")

  on.exit({
    unlink(temp_fastq)
    unlink(temp_dir1, recursive = TRUE)
    unlink(temp_dir2, recursive = TRUE)
  })

  result1 <- display_read_accuracy(temp_fastq, output_dir = temp_dir1, output_format = "plot")
  result2 <- display_read_accuracy(temp_fastq, output_dir = temp_dir2, output_format = "plot")

  expect_true(file.exists(file.path(temp_dir1, "accuracy_distribution.png")))
  expect_true(file.exists(file.path(temp_dir2, "accuracy_distribution.png")))
})

# ============================================================================
# Statistical Accuracy Tests
# ============================================================================

test_that("display_read_accuracy calculates reasonable accuracy values", {
  skip_if_not(has_dependencies, "minimap2 or samtools not available")

  temp_fastq <- create_test_fastq(20)
  on.exit(unlink(temp_fastq))

  result <- display_read_accuracy(temp_fastq, output_format = "table", table_format = "data.frame")

  mean_acc <- as.numeric(result$Value[result$Metric == "Mean Accuracy"])

  # For synthetic data with few errors, should be > 0.95
  expect_true(mean_acc > 0.95)
  expect_true(mean_acc <= 1.0)
})

test_that("display_read_accuracy table contains expected metrics", {
  skip_if_not(has_dependencies, "minimap2 or samtools not available")

  temp_fastq <- create_test_fastq(10)
  on.exit(unlink(temp_fastq))

  result <- display_read_accuracy(temp_fastq, output_format = "table", table_format = "data.frame")

  expected_metrics <- c(
    "Mean Accuracy", "Median Accuracy", "Standard Deviation",
    "Minimum", "Maximum", "5th Percentile", "25th Percentile",
    "75th Percentile", "95th Percentile", "Number of Reads"
  )

  expect_true(all(expected_metrics %in% result$Metric))
})

test_that("display_read_accuracy returns NULL for empty alignment", {
  skip_if_not(has_dependencies, "minimap2 or samtools not available")

  # Create a FASTQ with sequences that won't align
  temp_fastq <- tempfile(fileext = ".fastq")
  on.exit(unlink(temp_fastq))

  # Random sequences unlikely to align to DNA CS
  writeLines(c(
    "@read_1",
    "NNNNNNNNNNNNNNNNNNNN",
    "+",
    "IIIIIIIIIIIIIIIIIIII"
  ), temp_fastq)

  expect_output(
    result <- display_read_accuracy(temp_fastq),
    "No aligned reads found"
  )

  expect_null(result)
})

# ============================================================================
# Threading Tests
# ============================================================================

test_that("display_read_accuracy accepts different thread counts", {
  skip_if_not(has_dependencies, "minimap2 or samtools not available")

  temp_fastq <- create_test_fastq(10)
  on.exit(unlink(temp_fastq))

  result1 <- display_read_accuracy(temp_fastq, threads = 1, output_format = "plot")
  result2 <- display_read_accuracy(temp_fastq, threads = 4, output_format = "plot")

  expect_s3_class(result1, "ggplot")
  expect_s3_class(result2, "ggplot")
})

# ============================================================================
# Edge Cases
# ============================================================================

test_that("display_read_accuracy handles very small FASTQ files", {
  skip_if_not(has_dependencies, "minimap2 or samtools not available")

  temp_fastq <- create_test_fastq(1)
  on.exit(unlink(temp_fastq))

  result <- display_read_accuracy(temp_fastq, output_format = "plot")

  expect_s3_class(result, "ggplot")
})

test_that("display_read_accuracy plot has correct labels", {
  skip_if_not(has_dependencies, "minimap2 or samtools not available")

  temp_fastq <- create_test_fastq(10)
  on.exit(unlink(temp_fastq))

  result <- display_read_accuracy(temp_fastq, output_format = "plot")

  expect_equal(result$labels$title, "Distribution of Nanopore Read Accuracies")
  expect_equal(result$labels$x, "Read Accuracy")
  expect_equal(result$labels$y, "Frequency")
  expect_equal(result$labels$caption, "Analysis of DNA CS standard")
})

# ============================================================================
# Integration Test
# ============================================================================

test_that("display_read_accuracy full workflow with all outputs", {
  skip_if_not(has_dependencies, "minimap2 or samtools not available")

  temp_fastq <- create_test_fastq(15)
  temp_dir <- tempfile()

  on.exit({
    unlink(temp_fastq)
    unlink(temp_dir, recursive = TRUE)
  })

  result <- display_read_accuracy(
    temp_fastq,
    output_dir = temp_dir,
    threads = 2,
    output_format = "both",
    table_format = "data.frame",
    theme = "flatly"
  )

  # Check return structure
  expect_type(result, "list")
  expect_named(result, c("plot", "table"))
  expect_s3_class(result$plot, "ggplot")
  expect_s3_class(result$table, "data.frame")

  # Check files created
  expect_true(file.exists(file.path(temp_dir, "accuracy_distribution.png")))
  expect_true(file.exists(file.path(temp_dir, "accuracy_summary.csv")))
  expect_true(file.exists(file.path(temp_dir, "read_accuracies.txt")))

  # Check table content
  expect_equal(nrow(result$table), 10)
  expect_true(all(c("Metric", "Value") %in% names(result$table)))
})
