
# SWGACNV

`SWGACNV` is an R package designed to facilitate Copy Number Variation (CNV) detection from BAM files 
generated through SWGA sequencing, specifically in malaria-related studies. 

## Installation

You can install the development version of SWGACNV like so:

``` r
install.packages("devtools")
devtools::install_github("Rocamadourr/SWGACNV")

```

## Function Overview

### `bam_to_csv()`

Converts a folder of BAM files issued from SWGA into CSV files formatted for CNV analysis.
It can also be done manually with [`SAMtools`](http://www.htslib.org/).

#### **Parameters:**

- `bam_folder`: Path to the folder containing the BAM files issued from SWGA.
- `gff_path`: Path to the reference GFF file (e.g., from [MalariaGEN](https://www.malariagen.net/data_package/open-dataset-plasmodium-falciparum-v70/) under the section "annotation").
- `output_folder`: *(optional)* Path to the output folder. Defaults to the working directory if not provided.

#### **Output:**

A `.csv` file for each BAM input, containing:
- Chromosome (`seqnames`)
- Position (`pos`)
- Number of reads (`count`)

These CSV files should then be cleaned using the `cds_cleaning()` function.

### `cds_cleaning()`

Cleans coverage CSV files by removing non-coding regions, keeping only coding sequences (CDS) of the genome.
This function is optional but we strongly recommend using it after the `bam_to_csv` function to eliminate any background noise from non-coding region of the genome.

#### **Parameters:**

- `csv_folder`: Path to the folder containing samples CSV files converted from BAM. Each file must include at least the following columns: `seqnames`, `pos`, and `count`.
- `gff_path`: Path to the reference GFF file (e.g., from [MalariaGEN](https://www.malariagen.net/data_package/open-dataset-plasmodium-falciparum-v70/) under the section "annotation").
- `output_folder`: *(optional)* Path to the output folder. Defaults to the working directory if not provided.

#### **Output:**

A `.csv` file per sample, containing the same columns as before but only for the CDS of the genome.
These CSV files should then be used as input for the `cnv_analysis()` and the `create_profile()` function.

### `create_profile()`

Generates SWGA profiles from multiple coverage CSV files for CNV analyze with the `cnv_analysis()` function.

#### **Parameters:**

- `profile_csv_folder`: Path to a folder containing the CSV files. Each file should have 3 columns: `seqnames`, `pos`, and `count`.
- `gene_position` *optional* Path to the file containing each gene to analyze with their start and end. It should have 3 column : `gene`, `start`, `end`.
- `chromosome` *optional* Vector of chromosome numbers to analyze (e.g., c(1, 2, 3)). Defaults to 1:14.
- `output_folder`: *optional* Path to the folder where the profile files will be saved. Defaults to the working directory if not provided.

#### **Output:**

Creates one profile file per chromosome (e.g., `profile01.csv`, `profile02.csv`, ..., `profile14.csv`), where each file contains 2 columns :
- Gene : gene name.
- Mean_profile : Mean normalized coverage per gene across all samples used in input.

These profiles are to be used with the `cnv_analysis()` function for CNV detection.


### `cnv_analysis()`

Performs CNV detection by comparing new SWGA samples to a reference profile created using `create_profile()`.

#### **Parameters:**

- `csv_folder`: Path to the folder containing sample CSV files with columns: `seqnames`, `pos`, and `count`.
- `region`: Name(s) of the reference region(s) to compare against. Either a string (e.g., `"BENIN"`) or a vector of strings (e.g., `c("BENIN", "TOGO")`).
- `chr`: ID of the chromosome to analyze (e.g., `"Pf3D7_01_v3"`, `"Pf3D7_02_v3"`...,`"Pf3D7_14_v3"`).
- `profile_folder` *(optionnal)* Path to the folder which contains a profile. This is useful if you want analyze your own sample against your own profile.
- `output_folder`: *(optional)* Path to the folder where the results will be saved. Defaults to the working directory if not provided.
#### **Output:**

- A CSV file with a ratio of the expression level of the new sample compared to the profile, and a z-score for each gene and sample.
  A higher ratio indicates greater gene expression relative to the profile while a lower ratio suggest the opposite.
  The z-score indicates how the gene is expressed compared to the profile. Positive z-scores reflect gene overexpression while negative values suggest underexpression relative to the profile.  
- One JPG plot per sample per region showing CNV highlights with gene-level annotations.

These results help identify amplifications or deletions in sample genomes based on significant z-scores (typically >2 or <-2).


## Exemple

```r
library(SWGACNV)

# Step 1: Convert BAM files to coverage CSV
bam_folder <- "data/bam_files"
gff_path <- "data/Pfalciparum_annotation.gff"
output_folder <- "results/csv_files"

bam_to_csv(bam_folder, gff_path, output_folder)

# Step 2: Create a reference profile from coverage data
profile_csv_folder <- "results/csv_files"
output_path <- "results/profiles"

create_profile(profile_csv_folder, output_path)

# Step 3: Detect CNVs in new samples
csv_folder <- "data/new_samples_csv"
region <- "BENIN"
chr <- "Pf3D7_01_v3"
output_folder <- "results/cnv_plots"

cnv_analysis(csv_folder, region, chr, output_folder)

```

## Dependencies

This package depends on:

- `Rsamtools`
- `GenomicRanges`
- `GenomicFeatures`
- `txdbmaker`
- `IRanges`
- `dplyr`

## Author

Created by Noé MATHIEUX and Romain COPPÉE.<br>
GitHub: [Rocamadourr²](https://github.com/Rocamadourr)
