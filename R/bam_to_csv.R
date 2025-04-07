
#' @title Converting BMA to csv for CNV identification.
#' @description This function convert BAM files issued from a Swga sequencing to a specific csv format (chromosome/position/number of reads), rendering it usable by the cnv_analysis function for CNV identification.
#'
#' @param bam_folder Path to the folder containing BAM files from SWGA.
#' @param gff_path The path to the gff file. It can be downloaded on MalariaGen website.
#' @param output_folder (optional) Path to the output folder. Defaults to the working directory if not provided.
#'
#' @return A csv file for each BAM, usable by the cnv_analysis function.
#'
#' @import dplyr
#' @importFrom Rsamtools BamFile ScanBamParam indexBam pileup PileupParam
#' @importFrom GenomicFeatures genes
#' @importFrom GenomicRanges seqnames
#' @importFrom IRanges IRanges IRangesList
#' @importFrom txdbmaker makeTxDbFromGFF
#'
#' @export

bam_to_csv <- function(bam_folder, gff_path, output_folder = NULL) {
  if (!dir.exists(bam_folder)) stop("bam_folder not found.")
  bam_files <- list.files(bam_folder, pattern = "\\.bam$", full.names = TRUE) # Loading BAM files

  if (!file.exists(gff_path)) stop("GFF file not found.")
  # Output path
  if (is.null(output_folder)) {
    output_folder <- file.path(getwd())  # Using working directory by default
  }else{
    if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  }


  txdb <- txdbmaker::makeTxDbFromGFF(gff_path, format = "gff")
  genes_ranges <- genes(txdb)

  iranges_list <- split(genes_ranges, seqnames(genes_ranges))
  iranges_list <- lapply(iranges_list, function(gr) {
    IRanges(start = start(gr), end = end(gr))
  })
  iranges_list <- IRangesList(iranges_list)

  which <- iranges_list
  scanbam_paramazer <- ScanBamParam(which = which)
  pileup_paramsazer <- PileupParam(
    min_nucleotide_depth = 0,
    distinguish_nucleotides = FALSE,
    distinguish_strands = FALSE
  )

  for (bam_file in bam_files) {
    cat("Treating :", bam_file, "\n")
    bam <- BamFile(bam_file)

    # Verify if the index file exist and create it if it does not.
    bai_file <- paste0(bam_file, ".bai")
    if (!file.exists(bai_file)) {
      cat("Indexing BAM file :", bam_file, "\n")
      indexBam(bam_file)
      if (!file.exists(bai_file)) {
        stop("ERROR : index file could not be properly generated for ",
             bam_file)
      }
      cat("Index generated for :", bam_file, "\n")
    } else {
      cat("Index already existing for :", bam_file, "\n")
    }

    rm(bam)
    gc()  # Cleaning the memory to enforce the reloading of the BAI file.

    # Verifying that the index is recognized.
    if (!file.exists(bai_file)) {
      stop("Index file could not be found even after creation for :",
           bam_file)
    }


    pileup_data <- pileup(
      bam_file,
      index = bai_file,
      scanBamParam = scanbam_paramazer,
      pileupParam = pileup_paramsazer
    )

    # Cleaning the columns.
    pileup_data <- pileup_data %>% select(-which_label)

    # Deleting duplicates.
    pileup_data <- pileup_data %>%
      group_by(seqnames, pos) %>%
      slice(1) %>%
      ungroup()


    output_coverage <- file.path(output_folder, paste0(basename(bam_file), ".csv"))
    write.csv(pileup_data, output_coverage, row.names = FALSE, quote = TRUE)
    cat("Filed saved at :", output_coverage, "\n")
  }

  cat("All BAM were converted succesfully.\n")

  # Clean the environment
  rm(
    bam_files,
    txdb,
    genes_ranges,
    iranges_list,
    which,
    scanbam_paramazer,
    pileup_paramsazer,
    bam_file,
    bai_file,
    pileup_data,
    output_coverage
  )
  gc()

}
