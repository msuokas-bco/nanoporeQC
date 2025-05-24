test_that("plot_read_accuracy handles missing files correctly", {
  expect_error(
    plot_read_accuracy("nonexistent.fastq"),
    "FASTQ file does not exist"
  )
})

test_that("plot_read_accuracy validates NULL input", {
  expect_error(
    plot_read_accuracy(NULL),
    "FASTQ file must be a single file path string"
  )
})
