
#' @title Profile creation for the cnv_analysis() function.
#' @description This function allow the user to create his own profile to use in the cnv_analysis function.
#'
#' @param profile_csv_folder Path to a folder containing the CSV files. Each file should have 3 columns : `seqnames`, `pos`, and `count`.
#' @param chromosomes (optional) Vector of chromosome numbers to analyze (e.g., 1:14 or c(1, 2, 3)). Defaults to 1:14.
#' @param gene_position (optional)Path to the file containing each gene to analyze with their start and end. It should have 3 column : `gene`, `start`, `end`.
#' @param output_folder (optional) Path to the folder where the profile files will be saved. Defaults to the working directory if not provided.
#'
#' @return A dataframe representing the SWGA profile of the chosen samples.
#'
#' @import dplyr
#' @importFrom utils install.packages read.csv write.csv
#' @importFrom purrr map_dfc
#' @importFrom pracma trapz
#' @importFrom tools file_path_sans_ext
#'
#' @export
create_profile <- function(profile_csv_folder, chromosomes = 1:14, gene_position = NULL, output_folder = NULL){

  # Checking if the path are correct
  if (!dir.exists(profile_csv_folder)) stop("profile_csv_folder not found.")
  # Chromosomes check
  if (!is.numeric(chromosomes) || any(chromosomes < 1 | chromosomes > 14 | chromosomes != as.integer(chromosomes))) {
    stop("'chromosomes' must be an integer vector between 1 and 14 (e.g., 1:14 or c(1, 3, 5)).")
  }



  # load file with genes start/end
  if (is.null(gene_position)) {
    # Load the csv file containing 3 columns : gene name, the positions at which they start and end.
    genes_file_path <- system.file("extdata", "genes_positions.csv", package = "SWGACNV")
    df_genes <- read.csv(genes_file_path)
  }else{
    if (!file.exists(gene_position)) stop("gene_position file not found.")
    df_genes <- read.csv(gene_position) ####### ATESTER #####################################################
  }

  # Load a list of the csv files path. Those files contain 3 columns : the chromosome, the position in the genome and the number of read.
  coverage_files <- list.files(profile_csv_folder, pattern = "\\.csv$", full.names = TRUE)# select all files ending by .csv



  # Loop to pre-load each files once.
  df_coverage_list <- list()
  for (cov_file in coverage_files) {
    cat("Pre-loading :", cov_file, "\n")
    cov_file_name <- file_path_sans_ext(basename(cov_file))  # Unique name.
    df_coverage_list[[cov_file_name]] <- read.csv(cov_file)
  }

  # Loop for each chromosome
  for (i in chromosomes){
    seqname_value <- paste0("Pf3D7_", sprintf("%02d", i), "_v3")
    chr_prefix <- sub("_v3$", "", seqname_value)
    chr_prefix <- toupper(chr_prefix)

    df_genes$gene <- toupper(df_genes$gene)
    df_genes_chr <- df_genes[grep(paste0("^", toupper(chr_prefix)), df_genes$gene), ]

    auc_profile <- data.frame(gene = df_genes_chr$gene)# append the Area Under the Curve for each samples in a column "AUC_samplename"

    cat ("Profiling :", seqname_value, "\n")

    # Calculate the AUC for each gene of each sample.
    for (cov_file in coverage_files) {
      auc_new_col <- sub("\\.sorted\\.bam_coverage$", "", file_path_sans_ext(basename(cov_file)))
      auc_profile[[auc_new_col]] <- rep(NA, nrow(auc_profile))

      cat ("Treating :", cov_file, "\n")
      cov_file_name <- file_path_sans_ext(basename(cov_file)) # Je suis pas sur de cette ligne a supp si ça bug
      df_coverage_temp <- df_coverage_list[[cov_file_name]] # Load the csv from the pre-made list.
      df_coverage_temp <- subset(df_coverage_temp, toupper(seqnames) == toupper(seqname_value))#Trop long

      for (y in 1:nrow(df_genes_chr)) {
        print(y)
        idx <- which(df_coverage_temp$pos >= df_genes_chr$start[y]
                     & df_coverage_temp$pos <= df_genes_chr$end[y]
        )
        auc_value <- trapz(df_coverage_temp$pos[idx], df_coverage_temp$count[idx])
        gene_length <- (df_genes_chr$end[y] - df_genes_chr$start[y] + 1)
        auc_profile[[auc_new_col]][y] <- auc_value / gene_length
      }
    }

    ################## calcule du cnv_analysis profile ##################

    # Correcting sample bias.
    GlobalMeans <- colMeans(auc_profile[, 2:ncol(auc_profile)])
    Correction_Factors <- mean(GlobalMeans) / GlobalMeans
    # Applying the correction to each column.
    auc_profile <- auc_profile %>%
      mutate(across(2:ncol(.), ~ . * Correction_Factors[cur_column()]))

    # Calculating the mean.
    auc_profile <- auc_profile %>%
      rowwise() %>%  # Apply the mean row by row
      mutate(
        Mean_Profile = mean(c_across(-1), na.rm = TRUE),  # Mean of every columns except the first one for every row.
      ) %>%
      ungroup()  #Get out of rowwise.
    auc_profile <- auc_profile %>% select(gene, Mean_Profile)


    # Saving the output file
    if (is.null(output_folder)) {
      output_file <- file.path(getwd(), paste0("profile", sprintf("%02d", i), ".csv"))  # Saving in the wd by default
    } else {
      if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
      output_file <- file.path(output_folder, paste0("profile", sprintf("%02d", i), ".csv"))
    }
    write.csv(auc_profile, output_file, row.names = FALSE)
    cat("Filed saved as :", output_file, "\n")
  }

  # Clean the environment
  rm(
    genes_file_path,
    df_genes,
    coverage_files,
    df_coverage_list,
    i,
    seqname_value,
    chr_prefix,
    df_genes_chr,
    auc_profile,
    cov_file,
    auc_new_col,
    cov_file_name,
    df_coverage_temp,
    y,
    idx,
    auc_value,
    gene_length,
    GlobalMeans,
    Correction_Factors,
    output_file
  )
  gc()

}
