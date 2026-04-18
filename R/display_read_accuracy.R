#' Display Nanopore Read Accuracy
#'
#' Analyzes basecalled FASTQ files from Oxford Nanopore sequencers to calculate and visualize
#' read accuracy distributions by aligning reads to a DNA CS standard reference sequence.
#' Results can be displayed as plots, tables, or both with unified professional themes.
#'
#' @param fastq Character string. Path to the input FASTQ file.
#' @param output_dir Character string or NULL. Directory to save output files.
#'   If NULL (default), no files are saved to disk.
#' @param threads Integer. Number of threads to use for alignment.
#'   Default is 2.
#' @param output_format Character string. Format for output: "plot" (default), "table", or "both".
#'   "plot" returns a ggplot2 histogram, "table" returns a summary statistics table,
#'   "both" returns a list with both elements.
#' @param table_format Character string. Format for table output: "kable" (default) or "data.frame".
#'   Only used when output_format is "table" or "both".
#' @param theme Character string. Theme for plot and table styling. Available themes: "sandstone" (default),
#'   "simplex", "flatly", "journal", "lumen", "spacelab", "united", "default", "cerulean", "cosmo",
#'   "litera", "lux", "materia", "minty", "zephyr".
#'   The theme affects both the plot colors (histogram fill, lines, text) and table appearance.
#'
#' @return Depends on output_format:
#' \itemize{
#'   \item "plot": A ggplot2 object showing the distribution of read accuracies
#'   \item "table": An HTML table (browsable) or data.frame with summary statistics
#'   \item "both": A list with $plot and $table elements
#'   \item NULL if no aligned reads are found
#' }
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Creates a temporary reference sequence (DNA CS standard)
#'   \item Aligns FASTQ reads using minimap2
#'   \item Calculates per-read accuracy based on alignment length and edit distance
#'   \item Creates output based on specified format (histogram and/or table)
#' }
#'
#' Accuracy is calculated as: (alignment_length - edit_distance) / alignment_length
#' where alignment_length is derived from the CIGAR string and edit_distance is from
#' the NM tag. This method correctly accounts for substitutions, insertions, and deletions.
#'
#' The function requires external tools minimap2 and samtools to be installed
#' and available in the system PATH.
#'
#' When output_format includes "table", the following statistics are provided:
#' mean, median, standard deviation, min, max, and various percentiles (5th, 25th, 75th, 95th).
#'
#' The theme parameter creates visual consistency between plots and tables, offering 15
#' professional themes from the Quarto/Bootswatch collection. Each theme's colors are
#' algorithmically derived from a primary brand color to ensure visual harmony.
#'
#' @examples
#' \dontrun{
#' # Basic usage - returns ggplot object
#' p <- display_read_accuracy("sample.fastq")
#' print(p)
#'
#' # Get summary table only
#' table <- display_read_accuracy("sample.fastq", output_format = "table")
#' print(table)
#'
#' # Get both plot and table with matching flatly theme
#' results <- display_read_accuracy("sample.fastq",
#'                                  output_format = "both",
#'                                  theme = "flatly")
#' print(results$plot)
#' print(results$table)
#'
#' # Get data.frame for programmatic use
#' df <- display_read_accuracy("sample.fastq",
#'                             output_format = "table",
#'                             table_format = "data.frame")
#'
#' # Save results to directory with simplex theme
#' display_read_accuracy("sample.fastq",
#'                      output_dir = "./results/",
#'                      output_format = "both",
#'                      theme = "simplex")
#'
#' # Use in Quarto/R Markdown with matching theme
#' # In chunk with results='asis':
#' results <- display_read_accuracy("sample.fastq",
#'                                 output_format = "both",
#'                                 theme = "sandstone")
#' }
#'
#' @import ggplot2
#' @importFrom Rsamtools ScanBamParam scanBam
#' @importFrom stats median sd quantile
#' @importFrom utils write.table
#' @importFrom htmltools HTML browsable
#' @export

display_read_accuracy <- function(
    fastq,
    output_dir = NULL,
    threads = 2,
    output_format = "plot",
    table_format = "kable",
    theme = "sandstone"
) {

  # Validate output_format parameter
  output_format <- match.arg(output_format, c("plot", "table", "both"))
  table_format <- match.arg(table_format, c("kable", "data.frame"))
  theme <- match.arg(theme, names(THEME_SPECS))

  # Check for htmltools if needed
  if (output_format %in% c("table", "both") && table_format == "kable") {
    if (!requireNamespace("htmltools", quietly = TRUE)) {
      stop("Package 'htmltools' is required for HTML table output. Install it with: install.packages('htmltools')")
    }
  }

  # Check fastq file
  if (is.null(fastq) || !is.character(fastq) || length(fastq) != 1) {
    stop("FASTQ file must be a single file path string")
  }

  fastq <- path.expand(fastq)

  if (!file.exists(fastq)) {
    stop("FASTQ file does not exist")
  }

  # Set default output directory (only if specified)
  if (!is.null(output_dir)) {
    if (!endsWith(output_dir, "/")) {
      output_dir <- paste0(output_dir, "/")
    }
    # Create directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
  }

  # Step 1: Prepare temp reference sequence
  ref_fasta <- tempfile(fileext = ".fasta")
  reference_seq <- paste0(
    "GCCATCAGATTGTGTTTGTTAGTCGCTTTTTTTTTTTGGAATTTTTTTTTTGGAATTTTTTTTTTGCGCTAACAACCTCCTGCCGTTTTGCCCGTGCATATCGGTCACGAACAAATCTGATTACTAAACACAGTAGCCTGGATTTGTTCTATCAGTAATCGACCTTATTCCTAATTAAATAGAGCAAATCCCCTTATTGGGGGTAAGACATGAAGATGCCAGAAAAACATGACCTGTTGGCCGCCATTCTCGCGGCAAAGGAACAAGGCATCGGGGCAATCCTTGCGTTTGCAATGGCGTACCTTCGCGGCAGATATAATGGCGGTGCGTTTACAAAAACAGTAATCGACGCAACGATGTGCGCCATTATCGCCTAGTTCATTCGTGACCTTCTCGACTTCGCCGGACTAAGTAGCAATCTCGCTTATATAACGAGCGTGTTTATCGGCTACATCGGTACTGACTCGATTGGTTCGCTTATCAAACGCTTCGCTGCTAAAAAAGCCGGAGTAGAAGATGGTAGAAATCAATAATCAACGTAAGGCGTTCCTCGATATGCTGGCGTGGTCGGAGGGAACTGATAACGGACGTCAGAAAACCAGAAATCATGGTTATGACGTCATTGTAGGCGGAGAGCTATTTACTGATTACTCCGATCACCCTCGCAAACTTGTCACGCTAAACCCAAAACTCAAATCAACAGGCGCCGGACGCTACCAGCTTCTTTCCCGTTGGTGGGATGCCTACCGCAAGCAGCTTGGCCTGAAAGACTTCTCTCCGAAAAGTCAGGACGCTGTGGCATTGCAGCAGATTAAGGAGCGTGGCGCTTTACCTATGATTGATCGTGGTGATATCCGTCAGGCAATCGACCGTTGCAGCAATATCTGGGCTTCACTGCCGGGCGCTGGTTATGGTCAGTTCGAGCATAAGGCTGACAGCCTGATTGCAAAATTCAAAGAAGCGGGCGGAACGGTCAGAGAGATTGATGTATGAGCAGAGTCACCGCGATTATCTCCGCTCTGGTTATCTGCATCATCGTCTGCCTGTCATGGGCTGTTAATCATTACCGTGATAACGCCATTACCTACAAAGCCCAGCGCGACAAAAATGCCAGAGAACTGAAGCTGGCGAACGCGGCAATTACTGACATGCAGATGCGTCAGCGTGATGTTGCTGCGCTCGATGCAAAATACACGAAGGAGTTAGCTGATGCTAAAGCTGAAAATGATGCTCTGCGTGATGATGTTGCCGCTGGTCGTCGTCGGTTGCACATCAAAGCAGTCTGTCAGTCAGTGCGTGAAGCCACCACCGCCTCCGGCGTGGATAATGCAGCCTCCCCCCGACTGGCAGACACCGCTGAACGGGATTATTTCACCCTCAGAGAGAGGCTGATCACTATGCAAAAACAACTGGAAGGAACCCAGAAGTATATTAATGAGCAGTGCAGATAGAGTTGCCCATATCGATGGGCAACTCATGCAATTATTGTGAGCAATACACACGCGCTTCCAGCGGAGTATAAATGCCTAAAGTAATAAAACCGAGCAATCCATTTACGAATGTTTGCTGGGTTTCTGTTTTAACAACATTTTCTGCGCCGCCACAAATTTTGGCTGCATCGACAGTTTTCTTCTGCCCAATTCCAGAAACGAAGAAATGATGGGTGATGGTTTCCTTTGGTGCTACTGCTGCCGGTTTGTTTTGAACAGTAAACGTCTGTTGAGCACATCCTGTAATAAGCAGGGCCAGCGCAGTAGCGAGTAGCATTTTTTTCATGGTGTTATTCCCGATGCTTTTTGAAGTTCGCAGAATCGTATGTGTAGAAAATTAAACAAACCCTAAACAATGAGTTGAAATTTCATATTGTTAATATTTATTAATGTATGTCAGGTGCGATGAATCGTCATTGTATTCCCGGATTAACTATGTCCACAGCCCTGACGGGGAACTTCTCTGCGGGAGTGTCCGGGAATAATTAAAACGATGCACACAGGGTTTAGCGCGTACACGTATTGCATTATGCCAACGCCCCGGTGCTGACACGGAAGAAACCGGACGTTATGATTTAGCGTGGAAAGATTTGTGTAGTGTTCTGAATGCTCTCAGTAAATAGTAATGAATTATCAAAGGTATAGTAATATCTTTTATGTTCATGGATATTTGTAACCCATCGGAAAACTCCTGCTTTAGCAAGATTTTCCCTGTATTGCTGAAATGTGATTTCTCTTGATTTCAACCTATCATAGGACGTTTCTATAAGATGCGTGTTTCTTGAGAATTTAACATTTACAACCTTTTTAAGTCCTTTTATTAACACGGTGTTATCGTTTTCTAACACGATGTGAATATTATCTGTGGCTAGATAGTAAATATAATGTGAGACGTTGTGACGTTTTAGTTCAGAATAAAACAATTCACAGTCTAAATCTTTTCGCACTTGATCGAATATTTCTTTAAAAATGGCAACCTGAGCCATTGGTAAAACCTTCCATGTGATACGAGGGCGCGTAGTTTGCATTATCGTTTTTATCGTTTCAATCTGGTCTGACCTCCTTGTGTTTTGTTGATGATTTATGTCAAATATTAGGAATGTTTTCACTTAATAGTATTGGTTGCGTAACAAAGTGCGGTCCTGCTGGCATTCTGGAGGGAAATACAACCGACAGATGTATGTAAGGCCAACGTGCTCAAATCTTCATACAGAAAGATTTGAAGTAATATTTTAACCGCTAGATGAAGAGCAAGCGCATGGAGCGACAAAATGAATAAAGAACAATCTGCTGATGATCCCTCCGTGGATCTGATTCGTGTAAAAAATATGCTTAATAGCACCATTTCTATGAGTTACCCTGATGTTGTAATTGCATGTATAGAACATAAGGTGTCTCTGGAAGCATTCAGAGCAATTGAGGCAGCGTTGGTGAAGCACGATAATAATATGAAGGATTATTCCCTGGTGGTTGACTGATCACCATAACTGCTAATCATTCAAACTATTTAGTCTGTGACAGAGCCAACACGCAGTCTGTCACTGTCAGGAAAGTGGTAAAACTGCAACTCAATTACTGCAATGCCCTCGTAATTAAGTGAATTTACAATATCGTCCTGTTCGGAGGGAAGAACGCGGGATGTTCATTCTTCATCACTTTTAATTGATGTATATGCTCTCTTTTCTGACGTTAGTCTCCGACGGCAGGCTTCAATGACCCAGGCTGAGAAATTCCCGGACCCTTTTTGCTCAAGAGCGATGTTAATTTGTTCAATCATTTGGTTAGGAAAGCGGATGTTGCGGGTTGTTGTTCTGCGGGTTCTGTTCTTCGTTGACATGAGGTTGCCCCGTATTCAGTGTCGCTGATTTGTATTGTCTGAAGTTGTTTTTACGTTAAGTTGATGCAGATCAATTAATACGATACCTGCGTCATAATTGATTATTTGACGTGGTTTGATGGCCTCCACGCACGTTGTGATATGTAGATGATAATCATTATCACTTTACGGGTCCTTTCCGGTGAAAAAAAAGGTACCAAAAAAAACATCGTCGTGAGTAGTGAACCGTAAGC"
  )
  writeLines(c(">DNA_CS_standard", reference_seq), ref_fasta)

  # Step 2: Run alignment
  bam_out <- tempfile(fileext = ".bam")

  # Check if minimap2 and samtools are available
  if (system("which minimap2", ignore.stdout = TRUE) != 0) {
    stop("minimap2 is not installed or not in PATH")
  }
  if (system("which samtools", ignore.stdout = TRUE) != 0) {
    stop("samtools is not installed or not in PATH")
  }

  # Run alignment with better error handling
  cmd_align <- sprintf(
    "minimap2 -ax map-ont -t %d %s %s | samtools view -b -F 2308 | samtools sort -o %s",
    threads,
    shQuote(ref_fasta),
    shQuote(fastq),
    shQuote(bam_out)
  )

  result <- system(cmd_align, ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (result != 0) {
    stop("Alignment failed with exit code: ", result)
  }

  # Index the BAM file
  index_result <- system(sprintf("samtools index %s", shQuote(bam_out)), ignore.stdout = TRUE)
  if (index_result != 0) {
    stop("BAM indexing failed with exit code: ", index_result)
  }

  # Step 3: Read alignment with corrected parameters
  param <- ScanBamParam(
    what = c("qname", "flag", "mapq", "cigar", "seq", "qwidth"),
    tag = "NM"
  )

  alignments <- scanBam(bam_out, param = param)[[1]]

  # Function to calculate alignment length from CIGAR string
  get_alignment_length <- function(cigar_string) {
    if (is.na(cigar_string)) return(NA)

    # Extract operations and their lengths
    ops <- gregexpr("[0-9]+[MIDNSHP=X]", cigar_string)
    op_strings <- regmatches(cigar_string, ops)[[1]]

    aln_len <- 0

    for (op_str in op_strings) {
      len <- as.numeric(sub("[MIDNSHP=X]", "", op_str))
      op <- sub("[0-9]+", "", op_str)

      # Count M, =, X, and D operations (aligned to reference)
      if (op %in% c("M", "=", "X", "D")) {
        aln_len <- aln_len + len
      }
      # Skip I, S, H (not aligned to reference)
    }

    return(aln_len)
  }

  accuracies <- numeric()

  for (i in seq_along(alignments$qname)) {
    cigar <- alignments$cigar[i]
    edit_distance <- alignments$tag$NM[i]

    if (!is.na(edit_distance) && !is.na(cigar)) {
      alignment_length <- get_alignment_length(cigar)

      if (!is.na(alignment_length) && alignment_length > 0) {
        matches <- alignment_length - edit_distance

        if (matches >= 0) {
          accuracies <- c(accuracies, matches / alignment_length)
        }
      }
    }
  }

  # Clean up temporary files
  file.remove(ref_fasta)
  file.remove(bam_out)
  if (file.exists(paste0(bam_out, ".bai"))) {
    file.remove(paste0(bam_out, ".bai"))
  }

  if (length(accuracies) > 0) {
    # Calculate statistics
    median_accuracy <- median(accuracies)
    mean_accuracy <- mean(accuracies)
    std_dev <- ifelse(length(accuracies) > 1, sd(accuracies), 0)
    fifth <- quantile(accuracies, probs = 0.05)

    # Calculate additional statistics for table
    q25 <- quantile(accuracies, probs = 0.25)
    q75 <- quantile(accuracies, probs = 0.75)
    q95 <- quantile(accuracies, probs = 0.95)
    min_accuracy <- min(accuracies)
    max_accuracy <- max(accuracies)

    # Initialize results list
    results <- list()

    # Define theme configuration (used for both plot and table)
    theme_config <- get_theme_config(theme)

    # Create plot if requested
    if (output_format %in% c("plot", "both")) {
      accuracy_df <- data.frame(accuracy = accuracies)

      x_min <- 0.95
      x_max <- 1.00

      ont_accuracy <- ggplot(accuracy_df, aes(x = accuracy)) +
        geom_histogram(binwidth = 0.001,
                       fill = theme_config$fill_color,
                       color = "grey20",
                       alpha = 0.7) +
        geom_vline(xintercept = mean_accuracy,
                   color = theme_config$mean_line_color,
                   linetype = "dashed",
                   linewidth = 1) +
        geom_vline(xintercept = fifth,
                   color = theme_config$percentile_line_color,
                   linetype = "dashed",
                   linewidth = 1) +
        scale_x_continuous(limits = c(x_min, x_max)) +
        annotate("text", x = mean_accuracy, y = Inf,
                 label = sprintf("Mean: %.4f", mean_accuracy),
                 vjust = -0.5, hjust = 1.1,
                 color = theme_config$mean_line_color,
                 size = 3, angle = 90) +
        annotate("text", x = fifth, y = Inf,
                 label = sprintf("5th perc: %.4f", fifth),
                 vjust = -0.5, hjust = 1.1,
                 color = theme_config$percentile_line_color,
                 size = 3, angle = 90) +
        labs(title = "Distribution of Nanopore Read Accuracies",
             subtitle = sprintf("Mean: %.4f, Median: %.4f, 5th perc: %.4f, N: %d",
                                mean_accuracy, median_accuracy, fifth, length(accuracies)),
             x = "Read Accuracy",
             y = "Frequency",
             caption = "Analysis of DNA CS standard") +
        theme_minimal() +
        theme(
          plot.title = element_text(face = "bold",
                                    size = 14,
                                    color = theme_config$title_color),
          plot.subtitle = element_text(size = 10,
                                       color = theme_config$text_color),
          axis.title = element_text(size = 10,
                                    color = theme_config$text_color),
          axis.text = element_text(color = theme_config$text_color),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey90"),
          plot.background = element_rect(fill = theme_config$bg_color,
                                         color = NA),
          panel.background = element_rect(fill = "white",
                                          color = NA),
          legend.position = "none"
        )

      results$plot <- ont_accuracy

      if (!is.null(output_dir)) {
        ggsave(paste0(output_dir, "accuracy_distribution.png"), ont_accuracy,
               width = 10, height = 6, dpi = 300)
      }
    }

    # Create table if requested
    if (output_format %in% c("table", "both")) {
      summary_table <- data.frame(
        Metric = c("Mean Accuracy", "Median Accuracy", "Standard Deviation",
                   "Minimum", "Maximum", "5th Percentile", "25th Percentile",
                   "75th Percentile", "95th Percentile", "Number of Reads"),
        Value = c(
          sprintf("%.4f", mean_accuracy),
          sprintf("%.4f", median_accuracy),
          sprintf("%.4f", std_dev),
          sprintf("%.4f", min_accuracy),
          sprintf("%.4f", max_accuracy),
          sprintf("%.4f", fifth),
          sprintf("%.4f", q25),
          sprintf("%.4f", q75),
          sprintf("%.4f", q95),
          as.character(length(accuracies))
        ),
        stringsAsFactors = FALSE
      )

      if (table_format == "kable") {
        table_id <- paste0("nanopore-table-", format(Sys.time(), "%Y%m%d%H%M%S"),
                           sample(1000:9999, 1))

        table_html <- sprintf(
          '<style>
  #%s {
    font-family: %s;
    border-collapse: collapse;
    width: 100%%;
    max-width: 600px;
    margin: 20px 0;
    box-shadow: 0 2px 8px rgba(0,0,0,0.1);
    border-radius: 8px;
    overflow: hidden;
  }
  #%s caption {
    font-size: 1.3em;
    font-weight: bold;
    color: %s;
    text-align: left;
    padding: 12px 0;
    margin-bottom: 10px;
  }
  #%s thead th {
    background-color: %s;
    color: %s;
    font-weight: 600;
    text-align: left;
    padding: 12px 15px;
    font-size: 0.95em;
    letter-spacing: 0.5px;
  }
  #%s tbody td {
    padding: 10px 15px;
    border-bottom: 1px solid #E0E0E0;
  }
  #%s tbody tr:nth-child(even) {
    background-color: %s;
  }
  #%s tbody tr:hover {
    background-color: %s;
    transition: background-color 0.2s ease;
  }
  #%s tbody td:first-child {
    font-weight: 500;
    color: #333;
  }
  #%s tbody td:last-child {
    text-align: right;
    font-family: "Courier New", monospace;
    color: %s;
    font-weight: 600;
  }
</style>
<table id="%s">
  <caption>Sequence Accuracy Summary Statistics</caption>
  <thead>
    <tr>
      <th>Metric</th>
      <th>Value</th>
    </tr>
  </thead>
  <tbody>
%s
  </tbody>
</table>',
          table_id, theme_config$font_family,
          table_id, theme_config$title_color,
          table_id, theme_config$header_bg, theme_config$header_color,
          table_id,
          table_id, theme_config$stripe_color,
          table_id, theme_config$hover_color,
          table_id,
          table_id, theme_config$value_color,
          table_id,
          paste(sprintf("    <tr>\n      <td>%s</td>\n      <td>%s</td>\n    </tr>",
                        summary_table$Metric, summary_table$Value), collapse = "\n")
        )

        table_output <- htmltools::browsable(htmltools::HTML(table_html))

      } else {
        table_output <- summary_table
      }

      results$table <- table_output

      if (!is.null(output_dir)) {
        write.csv(summary_table,
                  paste0(output_dir, "accuracy_summary.csv"),
                  row.names = FALSE)
      }
    }

    # Save raw accuracies if output_dir is specified
    if (!is.null(output_dir)) {
      write.table(accuracies, paste0(output_dir, "read_accuracies.txt"),
                  row.names = FALSE, col.names = FALSE)
      cat("Output saved to:", output_dir, "\n")
    }

    # Return appropriate output based on format
    if (output_format == "plot") {
      return(results$plot)
    } else if (output_format == "table") {
      return(results$table)
    } else {
      return(results)
    }

  } else {
    cat("No aligned reads found.\n")
    return(NULL)
  }
}
