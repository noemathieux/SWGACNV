
#' @title Profile creation for the zscoreProfile function.
#' @description This function allow the user to create his own profile to use in the zscoreProfile function.
#'
#' @param profile_csv_folder A folder containing the sample in csv format (seqnames/pos/count) used to create a new profile.
#' @param output_path (optionnal) Name and path of the output file.
#'
#' @return A dataframe representing the SWGA profile of the chosen samples.
#'
#' @importFrom utils install.packages read.csv write.csv
#' @importFrom dplyr select bind_cols %>%
#' @importFrom purrr map_dfc
#' @importFrom pracma trapz
#' @importFrom tools file_path_sans_ext
#'
#' @export
create_profile <- function(profile_csv_folder, output_path = NULL){

  # Load the csv file containing 3 columns : gene name, the positions at which they start and end.
  genes_file_path <- system.file("extdata", "genes_positions.csv", package = "SWGACNV")
  df_genes <- read.csv(genes_file_path)

  # Load a list of the csv files path. Thoses files contain 3 columns : the chromosome, the position in the genome and the number of read.
  coverage_files <- list.files(profile_csv_folder, pattern = "\\.csv$", full.names = TRUE)# select all files ending by .csv



  # Loop to pre-load each files once.
  df_coverage_list <- list()
  for (cov_file in coverage_files) {
    cat("Pre-loading :", cov_file, "\n")
    cov_file_name <- file_path_sans_ext(basename(cov_file))  # Unique name.
    df_coverage_list[[cov_file_name]] <- read.csv(cov_file)
  }

  # Loop for each chromosome
  for (i in 1:14){
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

    ################## calcule du zscore profile ##################

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
    if (is.null(output_path)) {
      output_file <- file.path(getwd(), paste0("profile", sprintf("%02d", i), ".csv"))  # Saving in the wd by default
    } else {
      output_file <- file.path(output_path, paste0("profile", sprintf("%02d", i), ".csv"))
    }
    write.csv(auc_profile, output_file, row.names = FALSE)
    cat("Filed saved as :", output_file, "\n")
  }
}
