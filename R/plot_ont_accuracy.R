#' Plot Nanopore Read Accuracy
#'
#' Analyzes basecalled FASTQ files from Oxford Nanopore sequencera to calculate and visualize
#' read accuracy distributions by aligning reads to a DNA CS standard reference sequence.
#'
#' @param fastq Character string. Path to the input FASTQ file.
#' @param output_dir Character string or NULL. Directory to save output files.
#'   If NULL (default), no files are saved to disk.
#' @param threads Integer. Number of threads to use for alignment.
#'   Default is parallel::detectCores() - 2.
#'
#' @return A ggplot2 object showing the distribution of read accuracies,
#'   or NULL if no aligned reads are found.
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Creates a temporary reference sequence (DNA CS standard)
#'   \item Aligns FASTQ reads using minimap2
#'   \item Calculates per-read accuracy based on edit distance
#'   \item Creates a histogram showing accuracy distribution with statistical annotations
#' }
#'
#' The function requires external tools minimap2 and samtools to be installed
#' and available in the system PATH.
#'
#' @examples
#' \dontrun{
#' # Basic usage - returns ggplot object
#' p <- plot_read_accuracy("sample.fastq")
#' print(p)
#'
#' # Save results to directory
#' p <- plot_read_accuracy("sample.fastq", output_dir = "./results/")
#'
#' # Use in R Markdown
#' plot_read_accuracy("sample.fastq")
#'
#' # Customize further
#' p <- plot_read_accuracy("sample.fastq")
#' p + ggplot2::theme_minimal()
#' }
#'
#' @import ggplot2
#' @import ggthemes
#' @import Rsamtools
#' @importFrom parallel detectCores
#' @export
#'
plot_read_accuracy <- function(
    fastq,
    output_dir = NULL,
    threads = parallel::detectCores() - 2
) {
  # Load required libraries
  if (!require(Rsamtools, quietly = TRUE)) {
    stop("Rsamtools package is required but not installed")
  }
  if (!require(ggplot2, quietly = TRUE)) {
    stop("ggplot2 package is required but not installed")
  }
  if (!require(ggthemes, quietly = TRUE)) {
    stop("ggthemes package is required but not installed")
  }

  if (is.null(fastq) || !is.character(fastq) || length(fastq) != 1) {
    stop("FASTQ file must be a single file path string")
  }

  if (!file.exists(fastq)) {
    stop("FASTQ file does not exist")
  }

  # Set default output directory (only if specified)
  if (!is.null(output_dir)) {
    if (!endsWith(output_dir, "/")) {
      output_dir <- paste0(output_dir, "/")
    }
  }

  # Step 1: Prepare temp reference
  ref_fasta <- tempfile(fileext = ".fasta")
  reference_seq <- paste0(
    "GCCATCAGATTGTGTTTGTTAGTCGCTTTTTTTTTTTGGAATTTTTTTTTTGGAATTTTTTTTTTGCGCTAACAACCTCCTGCCGTTTTGCCCGTGCATATCGGTCACGAACAAATCTGATTACTAAACACAGTAGCCTGGATTTGTTCTATCAGTAATCGACCTTATTCCTAATTAAATAGAGCAAATCCCCTTATTGGGGGTAAGACATGAAGATGCCAGAAAAACATGACCTGTTGGCCGCCATTCTCGCGGCAAAGGAACAAGGCATCGGGGCAATCCTTGCGTTTGCAATGGCGTACCTTCGCGGCAGATATAATGGCGGTGCGTTTACAAAAACAGTAATCGACGCAACGATGTGCGCCATTATCGCCTAGTTCATTCGTGACCTTCTCGACTTCGCCGGACTAAGTAGCAATCTCGCTTATATAACGAGCGTGTTTATCGGCTACATCGGTACTGACTCGATTGGTTCGCTTATCAAACGCTTCGCTGCTAAAAAAGCCGGAGTAGAAGATGGTAGAAATCAATAATCAACGTAAGGCGTTCCTCGATATGCTGGCGTGGTCGGAGGGAACTGATAACGGACGTCAGAAAACCAGAAATCATGGTTATGACGTCATTGTAGGCGGAGAGCTATTTACTGATTACTCCGATCACCCTCGCAAACTTGTCACGCTAAACCCAAAACTCAAATCAACAGGCGCCGGACGCTACCAGCTTCTTTCCCGTTGGTGGGATGCCTACCGCAAGCAGCTTGGCCTGAAAGACTTCTCTCCGAAAAGTCAGGACGCTGTGGCATTGCAGCAGATTAAGGAGCGTGGCGCTTTACCTATGATTGATCGTGGTGATATCCGTCAGGCAATCGACCGTTGCAGCAATATCTGGGCTTCACTGCCGGGCGCTGGTTATGGTCAGTTCGAGCATAAGGCTGACAGCCTGATTGCAAAATTCAAAGAAGCGGGCGGAACGGTCAGAGAGATTGATGTATGAGCAGAGTCACCGCGATTATCTCCGCTCTGGTTATCTGCATCATCGTCTGCCTGTCATGGGCTGTTAATCATTACCGTGATAACGCCATTACCTACAAAGCCCAGCGCGACAAAAATGCCAGAGAACTGAAGCTGGCGAACGCGGCAATTACTGACATGCAGATGCGTCAGCGTGATGTTGCTGCGCTCGATGCAAAATACACGAAGGAGTTAGCTGATGCTAAAGCTGAAAATGATGCTCTGCGTGATGATGTTGCCGCTGGTCGTCGTCGGTTGCACATCAAAGCAGTCTGTCAGTCAGTGCGTGAAGCCACCACCGCCTCCGGCGTGGATAATGCAGCCTCCCCCCGACTGGCAGACACCGCTGAACGGGATTATTTCACCCTCAGAGAGAGGCTGATCACTATGCAAAAACAACTGGAAGGAACCCAGAAGTATATTAATGAGCAGTGCAGATAGAGTTGCCCATATCGATGGGCAACTCATGCAATTATTGTGAGCAATACACACGCGCTTCCAGCGGAGTATAAATGCCTAAAGTAATAAAACCGAGCAATCCATTTACGAATGTTTGCTGGGTTTCTGTTTTAACAACATTTTCTGCGCCGCCACAAATTTTGGCTGCATCGACAGTTTTCTTCTGCCCAATTCCAGAAACGAAGAAATGATGGGTGATGGTTTCCTTTGGTGCTACTGCTGCCGGTTTGTTTTGAACAGTAAACGTCTGTTGAGCACATCCTGTAATAAGCAGGGCCAGCGCAGTAGCGAGTAGCATTTTTTTCATGGTGTTATTCCCGATGCTTTTTGAAGTTCGCAGAATCGTATGTGTAGAAAATTAAACAAACCCTAAACAATGAGTTGAAATTTCATATTGTTAATATTTATTAATGTATGTCAGGTGCGATGAATCGTCATTGTATTCCCGGATTAACTATGTCCACAGCCCTGACGGGGAACTTCTCTGCGGGAGTGTCCGGGAATAATTAAAACGATGCACACAGGGTTTAGCGCGTACACGTATTGCATTATGCCAACGCCCCGGTGCTGACACGGAAGAAACCGGACGTTATGATTTAGCGTGGAAAGATTTGTGTAGTGTTCTGAATGCTCTCAGTAAATAGTAATGAATTATCAAAGGTATAGTAATATCTTTTATGTTCATGGATATTTGTAACCCATCGGAAAACTCCTGCTTTAGCAAGATTTTCCCTGTATTGCTGAAATGTGATTTCTCTTGATTTCAACCTATCATAGGACGTTTCTATAAGATGCGTGTTTCTTGAGAATTTAACATTTACAACCTTTTTAAGTCCTTTTATTAACACGGTGTTATCGTTTTCTAACACGATGTGAATATTATCTGTGGCTAGATAGTAAATATAATGTGAGACGTTGTGACGTTTTAGTTCAGAATAAAACAATTCACAGTCTAAATCTTTTCGCACTTGATCGAATATTTCTTTAAAAATGGCAACCTGAGCCATTGGTAAAACCTTCCATGTGATACGAGGGCGCGTAGTTTGCATTATCGTTTTTATCGTTTCAATCTGGTCTGACCTCCTTGTGTTTTGTTGATGATTTATGTCAAATATTAGGAATGTTTTCACTTAATAGTATTGGTTGCGTAACAAAGTGCGGTCCTGCTGGCATTCTGGAGGGAAATACAACCGACAGATGTATGTAAGGCCAACGTGCTCAAATCTTCATACAGAAAGATTTGAAGTAATATTTTAACCGCTAGATGAAGAGCAAGCGCATGGAGCGACAAAATGAATAAAGAACAATCTGCTGATGATCCCTCCGTGGATCTGATTCGTGTAAAAAATATGCTTAATAGCACCATTTCTATGAGTTACCCTGATGTTGTAATTGCATGTATAGAACATAAGGTGTCTCTGGAAGCATTCAGAGCAATTGAGGCAGCGTTGGTGAAGCACGATAATAATATGAAGGATTATTCCCTGGTGGTTGACTGATCACCATAACTGCTAATCATTCAAACTATTTAGTCTGTGACAGAGCCAACACGCAGTCTGTCACTGTCAGGAAAGTGGTAAAACTGCAACTCAATTACTGCAATGCCCTCGTAATTAAGTGAATTTACAATATCGTCCTGTTCGGAGGGAAGAACGCGGGATGTTCATTCTTCATCACTTTTAATTGATGTATATGCTCTCTTTTCTGACGTTAGTCTCCGACGGCAGGCTTCAATGACCCAGGCTGAGAAATTCCCGGACCCTTTTTGCTCAAGAGCGATGTTAATTTGTTCAATCATTTGGTTAGGAAAGCGGATGTTGCGGGTTGTTGTTCTGCGGGTTCTGTTCTTCGTTGACATGAGGTTGCCCCGTATTCAGTGTCGCTGATTTGTATTGTCTGAAGTTGTTTTTACGTTAAGTTGATGCAGATCAATTAATACGATACCTGCGTCATAATTGATTATTTGACGTGGTTTGATGGCCTCCACGCACGTTGTGATATGTAGATGATAATCATTATCACTTTACGGGTCCTTTCCGGTGAAAAAAAAGGTACCAAAAAAAACATCGTCGTGAGTAGTGAACCGTAAGC"
  )
  writeLines(c(">DNA_CS_standard", reference_seq), ref_fasta)

  # Get reference length from sequence
  reference_length <- nchar(reference_seq)

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
    "minimap2 -ax map-ont -t %d %s %s | samtools view -b -F 4 | samtools sort -o %s",
    threads,
    ref_fasta,
    fastq,
    bam_out
  )

  result <- system(cmd_align, ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (result != 0) {
    stop("Alignment failed with exit code: ", result)
  }

  # Index the BAM file
  index_result <- system(sprintf("samtools index %s", bam_out), ignore.stdout = TRUE)
  if (index_result != 0) {
    stop("BAM indexing failed with exit code: ", index_result)
  }

  # Step 3: Read alignment with corrected parameters
  param <- Rsamtools::ScanBamParam(
    what = c("qname", "flag", "mapq", "cigar", "seq", "qwidth"),
    tag = "NM"
  )

  alignments <- Rsamtools::scanBam(bam_out, param = param)[[1]]

  # Initialize variables
  total_matched_bases <- 0
  total_bases <- 0
  accuracies <- numeric()

  # Process alignments with better error handling
  for (i in seq_along(alignments$qname)) {
    sequence_length <- alignments$qwidth[i]
    edit_distance <- alignments$tag$NM[i]

    if (!is.na(edit_distance) && !is.na(sequence_length) && sequence_length > 0) {
      matches <- sequence_length - edit_distance

      # Ensure matches is not negative
      if (matches >= 0) {
        read_accuracy <- matches / sequence_length
        accuracies <- c(accuracies, read_accuracy)
        total_matched_bases <- total_matched_bases + matches
        total_bases <- total_bases + sequence_length
      }
    }
  }

  # Calculate overall accuracy
  overall_accuracy <- ifelse(total_bases > 0,
                             total_matched_bases / total_bases, 0)

  # Clean up temporary files
  file.remove(ref_fasta)
  file.remove(bam_out)
  if (file.exists(paste0(bam_out, ".bai"))) {
    file.remove(paste0(bam_out, ".bai"))
  }

  if (length(accuracies) > 0) {
    median_accuracy <- median(accuracies)
    mean_accuracy <- mean(accuracies)
    std_dev <- ifelse(length(accuracies) > 1, sd(accuracies), 0)
    two_std_dev <- 2 * std_dev
    lower_bound <- max(0, mean_accuracy - two_std_dev)
    upper_bound <- min(1, mean_accuracy + two_std_dev)

    accuracy_df <- data.frame(accuracy = accuracies)

    # Create plot with dynamic x-axis limits
    x_min <- 0.95
    x_max <- 1.00

    ont_accuracy <- ggplot(accuracy_df, aes(x = accuracy)) +
      geom_histogram(binwidth = 0.001, fill = "dodgerblue1", color = "grey20", alpha = 0.7) +
      geom_vline(xintercept = mean_accuracy, color = "red4", linetype = "dashed", linewidth = 1) +
      geom_vline(xintercept = lower_bound, color = "navyblue", linetype = "dashed", linewidth = 1) +
      geom_vline(xintercept = upper_bound, color = "navyblue", linetype = "dashed", linewidth = 1) +
      scale_x_continuous(limits = c(x_min, x_max)) +
      annotate("text", x = mean_accuracy, y = Inf,
               label = sprintf("Mean: %.4f", mean_accuracy),
               vjust = -0.5, hjust = 1.1, color = "red4", size = 3, angle = 90) +
      annotate("text", x = lower_bound, y = Inf,
               label = sprintf("-2σ: %.4f", lower_bound),
               vjust = -0.5, hjust = 1.1, color = "navyblue", size = 3, angle = 90) +
      annotate("text", x = upper_bound, y = Inf,
               label = sprintf("+2σ: %.4f", upper_bound),
               vjust = 1.5, hjust = 1.1, color = "navyblue", size = 3, angle = 90) +
      labs(title = "Distribution of Nanopore Read Accuracies",
           subtitle = sprintf("Mean: %.4f, Median: %.4f, SD: %.4f, N: %d",
                              mean_accuracy, median_accuracy, std_dev, length(accuracies)),
           x = "Read Accuracy",
           y = "Frequency",
           caption = "Analysis of DNA CS standard") +
      theme_fivethirtyeight() +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10),
        axis.title = element_text(size = 10),
        panel.grid.minor = element_blank(),
        legend.position = "none"
      )

    # Save output only if output_dir is specified
    if (!is.null(output_dir)) {
      ggsave(paste0(output_dir, "accuracy_distribution.png"), ont_accuracy,
             width = 10, height = 6, dpi = 300)
      write.table(accuracies, paste0(output_dir, "read_accuracies.txt"),
                  row.names = FALSE, col.names = FALSE)
      cat("Plots and data saved to:", output_dir, "\n")
    }

    # Return the ggplot object directly
    return(ont_accuracy)

  } else {
    cat("No aligned reads found.\n")
    return(NULL)
  }
}
