% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_ont_accuracy.R
\name{plot_read_accuracy}
\alias{plot_read_accuracy}
\title{Plot Nanopore Read Accuracy}
\usage{
plot_read_accuracy(
  fastq,
  output_dir = NULL,
  threads = parallel::detectCores() - 2
)
}
\arguments{
\item{fastq}{Character string. Path to the input FASTQ file.}

\item{output_dir}{Character string or NULL. Directory to save output files.
If NULL (default), no files are saved to disk.}

\item{threads}{Integer. Number of threads to use for alignment.
Default is parallel::detectCores() - 2.}
}
\value{
A ggplot2 object showing the distribution of read accuracies,
  or NULL if no aligned reads are found.
}
\description{
Analyzes basecalled FASTQ files from Oxford Nanopore sequencera to calculate and visualize
read accuracy distributions by aligning reads to a DNA CS standard reference sequence.
}
\details{
This function performs the following steps:
\enumerate{
  \item Creates a temporary reference sequence (DNA CS standard)
  \item Aligns FASTQ reads using minimap2
  \item Calculates per-read accuracy based on edit distance
  \item Creates a histogram showing accuracy distribution with statistical annotations
}

The function requires external tools minimap2 and samtools to be installed
and available in the system PATH.
}
\examples{
\dontrun{
# Basic usage - returns ggplot object
p <- plot_read_accuracy("sample.fastq")
print(p)

# Save results to directory
p <- plot_read_accuracy("sample.fastq", output_dir = "./results/")

# Use in R Markdown
plot_read_accuracy("sample.fastq")

# Customize further
p <- plot_read_accuracy("sample.fastq")
p + ggplot2::theme_minimal()
}

}
