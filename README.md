# NanoporeQC package

A R package for analyzing and visualizing read accuracy distributions from Oxford Nanopore sequencing data. The package aligns basecalled FASTQ files to a DNA CS standard reference sequence and provides accuracy metrics through customizable plots and tables.

## Features

-   Automated Alignment: Uses minimap2 for read alignment to DNA CS standard
-   Comprehensive Statistics: Calculates mean, median, percentiles, and distribution metrics
-   Flexible Output Formats: Generate plots, tables, or both
-   Themes: Built-in themes (sandstone, simplex, flatly) with unified styling
-   Export Capabilities: Save results as PNG, CSV files and per-read accuracy values in a text file

## Installation

``` r
# Install from GitHub
# install.packages("devtools")
devtools::install_github("msuokas-bco/nanoporeQC")
```

### System Requirements

The package requires external tools to be installed and available in your environment:

-   minimap2: Read alignment tool
-   samtools: BAM file processing

Installation on Ubuntu/Debian (or alternatively use conda):

``` bash
sudo apt-get install minimap2 samtools
```

Installation on macOS (using Homebrew):

``` bash
brew install minimap2 samtools
```

Installation on Windows: Use WSL (Windows Subsystem for Linux) or install via conda.

## Quick Start

``` r
library(nanoporeQC)

# Basic usage - returns a ggplot2 histogram
plot <- display_read_accuracy("sample.fastq")
print(plot)

# Get summary statistics table
table <- display_read_accuracy("sample.fastq", output_format = "table")
print(table)

# Get both plot and table with matching theme
results <- display_read_accuracy(
  "sample.fastq",
  output_format = "both",
  theme = "flatly"
)
print(results$plot)
print(results$table)
```

## Usage

### Function Parameters

``` r
display_read_accuracy(
  fastq,
  output_dir = NULL,
  threads = 2,
  output_format = "plot",
  table_format = "kable",
  theme = "sandstone"
)
```

Parameters:

-   `fastq`: Path to the input FASTQ file (required)
-   `output_dir`: Directory to save output files (default: NULL, no files saved)
-   `threads`: Number of threads for alignment (default: 2)
-   `output_format`: Output type - "plot", "table", or "both" (default: "plot")
-   `table_format`: Table format - "kable" (HTML) or "data.frame" (default: "kable")
-   `theme`: Visual theme - "sandstone", "simplex", "flatly", "journal", "lumen", "spacelab" or "united" (default: "sandstone")

### Output Formats

#### Plot Output

Returns a ggplot2 histogram showing:

-   Distribution of read accuracies
-   Mean accuracy line (dashed)
-   5th percentile line (dashed)
-   Summary statistics in subtitle
-   Customizable theme colors

#### Table Output

Returns summary statistics including:

-   Mean, median, standard deviation
-   Min, max values
-   Percentiles (5th, 25th, 75th, 95th)
-   Total number of reads

When `table_format = "kable"`, returns an interactive HTML table with professional styling.

#### Both Output

Returns a list with `$plot` and `$table` elements.

### Examples

#### Save Results to Directory

``` r
display_read_accuracy(
  "sample.fastq",
  output_dir = "./results/",
  output_format = "both",
  threads = 4
)
```

This creates:

-   `accuracy_distribution.png` - High-resolution plot (300 DPI)
-   `accuracy_summary.csv` - Summary statistics table
-   `read_accuracies.txt` - Raw per-read accuracy values

#### Get Data Frame for Further Analysis

``` r
stats <- display_read_accuracy(
  "sample.fastq",
  output_format = "table",
  table_format = "data.frame"
)

# Access specific metrics
mean_acc <- as.numeric(stats$Value[stats$Metric == "Mean Accuracy"])
```

#### Use in Quarto/R Markdown

```` markdown
```{r, results='asis'}
results <- display_read_accuracy(
  "sample.fastq",
  output_format = "both",
  theme = "sandstone"
)

# Plot will render automatically
print(results$plot)

# Table will render as HTML
results$table
```
````

#### Compare Different Themes

``` r
library(patchwork)

p1 <- display_read_accuracy("sample.fastq", theme = "sandstone")
p2 <- display_read_accuracy("sample.fastq", theme = "simplex")
p3 <- display_read_accuracy("sample.fastq", theme = "flatly")
p4 <- display_read_accuracy("sample.fastq", theme = "journal")
p5 <- display_read_accuracy("sample.fastq", theme = "lumen")
p6 <- display_read_accuracy("sample.fastq", theme = "spacelab")
p7 <- display_read_accuracy("sample.fastq", theme = "united")


p1 / p2 / p3 / p4 / p5 / p6 / p7
```

## Themes

The package includes seven themes that are inspired by Bootstrap 5 themes also used by Quarto.

## Understanding the Analysis

### Workflow

1.  Reference Preparation: Creates temporary DNA CS standard reference sequence
2.  Read Alignment: Aligns FASTQ reads using minimap2 with ONT preset
3.  Accuracy Calculation: Computes per-read accuracy from edit distance (NM tag)
4.  Visualization: Generates customized plots and/or tables
5.  Export: Optionally saves results to disk

### Accuracy Calculation

Read accuracy is calculated as:

```         
accuracy = (alignment_length - edit_distance) / alignment_length
```

Where edit_distance includes mismatches, insertions, and deletions.

### Statistics Provided

-   Mean: Average accuracy across all reads
-   Median: Middle value of accuracy distribution
-   Standard Deviation: Measure of accuracy variability
-   Percentiles: 5th, 25th, 75th, 95th percentiles
-   Range: Minimum and maximum accuracies

## Troubleshooting

### "minimap2 is not installed or not in PATH"

Ensure minimap2 is installed and accessible:

``` bash
which minimap2  # Should return path to minimap2
```

### "samtools is not installed or not in PATH"

Ensure samtools is installed:

``` bash
which samtools  # Should return path to samtools
```

### "No aligned reads found"

Possible causes:

-   FASTQ file is empty or corrupted
-   Reads don't match the DNA CS standard reference
-   File format issues (ensure file is proper FASTQ format)

### Memory Issues with Large Files

Increase available threads for faster processing:

``` r
display_read_accuracy("large_file.fastq", threads = 8)
```

## Dependencies

### R Packages

-   ggplot2
-   ggthemes
-   Rsamtools
-   htmltools (for HTML table output)
-   stats
-   utils
-   parallel

### External Tools

-   minimap2 (≥2.17)
-   samtools (≥1.10)

## Citation

If you use this function in your research, please cite also the relevant tools: - minimap2: Li, H. (2018). Bioinformatics, 34(18), 3094-3100. - samtools: Li, H., et al. (2009). Bioinformatics, 25(16), 2078-2079.

## License

This project is licensed under the GNU General Public License v3.0 (GPL-3) - see the [LICENSE](LICENSE) file for details.
