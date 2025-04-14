#' @title Converting BAM to CSV without GFF, for CNV identification.
#' @description This function converts BAM files from SWGA sequencing into a simple csv format (chromosome / position / number of reads), filtering out low coverage positions. It can be used as input for the cnv_analysis function.
#'
#' @param bam_folder Path to the folder containing BAM files from SWGA.
#' @param output_folder (optional) Path to the output folder. Defaults to the working directory if not provided.
#' @param min_coverage (optional) Minimum read depth to keep a position. Default is 5.
#'
#' @return A csv file for each BAM, containing only positions with at least `min_coverage` reads.
#'
#' @import dplyr
#' @importFrom Rsamtools BamFile ScanBamParam indexBam pileup PileupParam
#'
#' @export
bam_to_csv_nogff <- function(bam_folder, output_folder = NULL, min_coverage = 5) {
  if (!dir.exists(bam_folder)) stop("bam_folder not found.")
  bam_files <- list.files(bam_folder, pattern = "\\.bam$", full.names = TRUE)

  if (is.null(output_folder)) {
    output_folder <- getwd()
  } else {
    if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  }

  scan_param <- ScanBamParam()
  pileup_param <- PileupParam(
    min_nucleotide_depth = 0,
    distinguish_nucleotides = FALSE,
    distinguish_strands = FALSE
  )

  for (bam_file in bam_files) {
    cat("Treating :", bam_file, "\n")
    bam <- BamFile(bam_file)
    bai_file <- paste0(bam_file, ".bai")

    if (!file.exists(bai_file)) {
      cat("Indexing BAM file :", bam_file, "\n")
      indexBam(bam_file)
      if (!file.exists(bai_file)) {
        stop("ERROR: Could not create index for ", bam_file)
      }
      cat("Index generated for :", bam_file, "\n")
    } else {
      cat("Index already exists for :", bam_file, "\n")
    }

    rm(bam)
    gc()

    pileup_data <- pileup(
      bam_file,
      index = bai_file,
      scanBamParam = scan_param,
      pileupParam = pileup_param
    )

    # Filtrage dynamique par seuil de couverture
    pileup_data <- pileup_data %>%
      filter(count >= min_coverage)

    output_coverage <- file.path(output_folder, paste0(basename(bam_file), ".csv"))
    write.csv(pileup_data, output_coverage, row.names = FALSE, quote = TRUE)
    cat("File saved at :", output_coverage, "\n")
  }

  cat("All BAM files were successfully converted.\n")

  rm(bam_files, scan_param, pileup_param, bam_file, bai_file, pileup_data, output_coverage)
  gc()
}
