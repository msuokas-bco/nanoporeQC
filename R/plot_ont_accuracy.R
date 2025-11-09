#' Plot Nanopore Read Accuracy (Deprecated)
#'
#' @description
#' \lifecycle{deprecated}
#'
#' This function has been deprecated and will be removed in a future version.
#' Please use \code{\link{display_read_accuracy}} instead, which better reflects
#' the expanded functionality (plots, tables, or both).
#'
#' @param fastq Character string. Path to the input FASTQ file.
#' @param output_dir Character string or NULL. Directory to save output files.
#' @param threads Integer. Number of threads to use for alignment.
#' @param output_format Character string. Format for output: "plot", "table", or "both".
#' @param table_format Character string. Format for table output: "kable" or "data.frame".
#' @param table_theme Character string. Theme for styling: "sandstone", "simplex", or "flatly".
#'
#' @return Same as \code{\link{display_read_accuracy}}
#'
#' @details
#' This function is a wrapper around \code{\link{display_read_accuracy}} and is
#' maintained for backward compatibility. It will be removed in version 1.0.0.
#'
#' Please update your code to use \code{display_read_accuracy()} instead:
#' \itemize{
#'   \item Old: \code{plot_read_accuracy("sample.fastq")}
#'   \item New: \code{display_read_accuracy("sample.fastq")}
#' }
#'
#' @examples
#' \dontrun{
#' # Deprecated - do not use
#' plot_read_accuracy("sample.fastq")
#'
#' # Use this instead
#' display_read_accuracy("sample.fastq")
#' }
#'
#' @seealso \code{\link{display_read_accuracy}} for the current function
#' @export

plot_read_accuracy <- function(
    fastq,
    output_dir = NULL,
    threads = 2,
    output_format = "plot",
    table_format = "kable",
    table_theme = "sandstone"
) {

  # Issue deprecation warning
  .Deprecated(
    new = "display_read_accuracy",
    package = "nanoporeQC",
    msg = paste(
      "plot_read_accuracy() is deprecated and will be removed in future versions.",
      "Please use display_read_accuracy() instead.",
      "The new name better reflects the function's ability to display plots, tables, or both.",
      sep = "\n  "
    )
  )

  # Call the new function
  display_read_accuracy(
    fastq = fastq,
    output_dir = output_dir,
    threads = threads,
    output_format = output_format,
    table_format = table_format,
    table_theme = table_theme
  )
}
