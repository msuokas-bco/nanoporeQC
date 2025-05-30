---
editor_options: 
  markdown: 
    wrap: 72
---

# NanoporeQC package

Package contains a R function for analyzing and visualizing read
accuracy distribution from Oxford Nanopore sequencing data by aligning
reads to a DNA CS (Control Standard) reference sequence.

## Overview

The `plot_read_accuracy()` function processes basecalled FASTQ files
from Oxford Nanopore sequencers to calculate per-read accuracy metrics
and generate publication-ready visualizations. It aligns sequencing
reads to a standardized DNA control sequence and computes accuracy based
on edit distance calculations.

## Features

-   **Automated alignment** using minimap2 with optimized ONT parameters
-   Supports both compressed and non-compressed fastq
-   **Per-read accuracy calculation** based on edit distance from
    reference
-   **Statistical analysis** including mean, median, and standard
    deviation
-   **Publication-ready plots** with customizable theming
-   **Flexible output options** - return plot or save to disk
-   **Parallel processing** support for faster alignment

## Prerequisites

Sequencing library needs to have DNA CS control included. FastQ file
needs to contain these sequences.

### R Packages

``` r
install.packages(c("ggplot2", "ggthemes"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rsamtools")
```

### External Tools

The function requires these command-line tools to be installed and
available in your system PATH:

-   **minimap2** - Fast sequence alignment tool
-   **samtools** - SAM/BAM file processing

#### Installation on different systems:

**Ubuntu/Debian:**

``` bash
sudo apt-get update
sudo apt-get install minimap2 samtools
```

**macOS (with Homebrew):**

``` bash
brew install minimap2 samtools
```

**Conda (cross-platform):**

``` bash
conda install -c bioconda minimap2 samtools
```

## Usage

### Basic Usage

``` r
# Load the function (assuming it's in your environment)
source("path/to/plot_read_accuracy.R")

# Analyze a FASTQ file and display the plot
plot <- plot_read_accuracy("sample.fastq")
print(plot)
```

### Save Results to Directory

``` r
# Analyze and save results to specified directory
plot_read_accuracy("sample.fastq", output_dir = "./results/")
```

### Customize Threading

``` r
# Use specific number of threads for alignment
plot_read_accuracy("sample.fastq", threads = 8)
```

### Integration with R Markdown

``` r
# Direct plotting in R Markdown documents
plot_read_accuracy("sample.fastq")
```

### Further Customization

``` r
# Get the plot object and customize further
p <- plot_read_accuracy("sample.fastq")
p + ggplot2::theme_minimal() + 
    ggplot2::labs(title = "Custom Title")
```

## Function Arguments

| Parameter | Type | Default | Description |
|------------------|----------------|----------------|----------------------|
| `fastq` | character | required | Path to input FASTQ file |
| `output_dir` | character/NULL | NULL | Directory to save output files (optional) |
| `threads` | integer | `detectCores() - 2` | Number of threads for alignment |

## Output

### Return Value

-   **ggplot2 object**: Histogram showing read accuracy distribution
    with statistical annotations
-   **NULL**: If no aligned reads are found

### Saved Files (when `output_dir` specified)

-   `accuracy_distribution.png`: High-resolution plot (300 DPI, 10×6
    inches)
-   `read_accuracies.txt`: Raw accuracy values for each read

### Plot Features

-   **Histogram**: Distribution of per-read accuracies (binwidth =
    0.001)
-   **Statistical lines**: Mean (red dashed) and ±2σ bounds (blue
    dashed)
-   **Annotations**: Mean, median, standard deviation, and sample size
-   **Professional theme**: Using ggthemes::theme_fivethirtyeight()

## Technical Details

### Analysis Workflow

1.  **Reference preparation**: Creates temporary FASTA file with DNA CS
    standard sequence
2.  **Sequence alignment**: Uses minimap2 with ONT-optimized parameters
    (`-ax map-ont`)
3.  **Accuracy calculation**: Computes per-read accuracy as
    `(read_length - edit_distance) / read_length`
4.  **Statistical analysis**: Calculates descriptive statistics and
    confidence bounds
5.  **Visualization**: Generates histogram with statistical annotations
6.  **Cleanup**: Removes temporary files automatically

### DNA CS Standard

The function uses a predefined DNA Control Standard sequence (3560 bp)
that serves as the reference for accuracy calculations. This
standardized approach ensures consistent and comparable results across
different sequencing runs.

### Alignment Parameters

-   **Algorithm**: minimap2 with ONT preset (`-ax map-ont`)
-   **Filtering**: Excludes unmapped reads (`-F 4`)
-   **Sorting**: Coordinate-sorted BAM output for efficient processing

## Example Output

The function generates a histogram showing: - X-axis: Read accuracy
(typically 0.95-1.00 range) - Y-axis: Frequency of reads - Red dashed
line: Mean accuracy - Blue dashed lines: ±2 standard deviation bounds -
Subtitle: Summary statistics (Mean, Median, SD, N)

Typical results show: - Mean accuracy: \~0.98-0.99 for recent ONT
chemistry - Standard deviation: \~0.01-0.02 - Distribution: Usually
right-skewed with peak near 0.99

## Error Handling

The function includes comprehensive error checking for: - Missing or
invalid input files - Unavailable external tools (minimap2, samtools) -
Failed alignment or indexing operations - Missing required R packages -
Invalid parameter values

## Performance Considerations

-   **Memory usage**: Scales with FASTQ file size and number of reads
-   **Processing time**: Depends on read count and available CPU cores
-   **Temporary files**: Created in system temp directory and cleaned up
    automatically
-   **Threading**: Automatically uses available cores minus 2 for system
    stability

## Troubleshooting

### Common Issues

**"minimap2 is not installed or not in PATH"** - Install minimap2 and
ensure it's in your system PATH - Test with `which minimap2` in terminal

**"samtools is not installed or not in PATH"** - Install samtools and
ensure it's in your system PATH - Test with `which samtools` in terminal

**"No aligned reads found"** - Check FASTQ file quality and format -
Verify reads are ONT format and not severely degraded - Consider
different reference sequence if using non-standard samples

**Memory issues with large files** - Reduce thread count to free up
memory - Process files in smaller batches - Consider using a machine
with more RAM

## Citation

If you use this function in your research, please cite the relevant
tools: - minimap2: Li, H. (2018). Bioinformatics, 34(18), 3094-3100. -
samtools: Li, H., et al. (2009). Bioinformatics, 25(16), 2078-2079.

## License

This project is licensed under the GNU General Public License v3.0
(GPL-3) - see the [LICENSE](LICENSE) file for details.
