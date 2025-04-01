#' @title Removing the sample introns.
#' @description Remove the introns from the sample in CSV format. It should be used after the `bam_to_csv` function for a more accurate analysis.
#'
#' @param bam_folder Path to the folder containing samples CSV files converted from BAM.
#' @param gff_path The path to the gff file. It can be downloaded on MalariaGen website.
#' @param output_folder (optional) Path to the output folder. Defaults to the working directory if not provided.
#'
#' @return A csv file without introns, usable by the cnv_analysis function.
#'
#' @import dplyr
#' @import GenomicRanges
#' @import rtracklayer
#'
#' @export

cds_cleaning <- function(csv_folder, gff_path, output_folder) {
  # Vérification des dossiers
  if (!dir.exists(csv_folder)) stop("csv_folder does not exist.")
  if (!file.exists(gff_path)) stop("GFF file not found.")
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

  # Output path
  if (is.null(output_folder)) {
    output_folder <- file.path(getwd())  # Using working directory by default
  }

  # Load the gff file
  gff <- rtracklayer::import(gff_path)
  gff_df <- as.data.frame(gff)

  # Select lines where CDS appear
  cds_df <- dplyr::filter(gff_df, type == "CDS")

  # Create a Grange item for the CDS
  gr_cds <- GenomicRanges::GRanges(
    seqnames = cds_df$seqnames,
    ranges = IRanges::IRanges(start = cds_df$start, end = cds_df$end)
  )

  # Loop on each csv file in the folder
  csv_files <- list.files(csv_folder, pattern = "\\.csv$", full.names = TRUE)
  for (file in csv_files) {
    csv <- tryCatch(read.csv(file), error = function(e) NULL)
    if (is.null(csv)) {
      warning(paste("Erreur de lecture du fichier :", file))
      next
    }

    if (!all(c("seqnames", "pos") %in% names(csv))) {
      warning(paste("Fichier ignoré (colonnes manquantes) :", file))
      next
    }

    # Créer GRanges pour le CSV
    gr_csv <- GenomicRanges::GRanges(
      seqnames = csv$seqnames,
      ranges = IRanges::IRanges(start = csv$pos, end = csv$pos)
    )

    # Trouver les overlaps
    hits <- GenomicRanges::findOverlaps(gr_csv, gr_cds)

    # Filtrer le CSV
    filtered_csv <- csv[queryHits(hits), ]

    # Nom de sortie
    output_name <- basename(file)
    output_file <- file.path(output_folder, sub("\\.csv$", "_CDS.csv", output_name))


    write.csv(filtered_csv, output_file, row.names = FALSE)
    message("Clean file saved at : ", output_file)
  }
}
