
#' @title Profile creation for the cnv_analysis() function.
#' @description This function allow the user to create his own profile to use in the cnv_analysis function.
#'
#' @param profile_csv_folder Path to a folder containing the CSV files. Each file should have 3 columns : `seqnames`, `pos`, and `count`.
#' @param chromosomes (optional) Vector of chromosome numbers to analyze (e.g., 1:14 or c(1, 4, 6)). Defaults to 1:14.
#' @param gene_position (optional)Path to the file containing each gene to analyze with their start and end. It should have 3 column : `gene`, `start`, `end`. Default to the one provided with the package.
#' @param output_folder (optional) Path to the folder where the profile files will be saved. Defaults to the working directory if not provided.
#'
#' @return A dataframe representing the SWGA profile of the chosen samples.
#'
#' @importFrom dplyr mutate across
#' @importFrom magrittr %>%
#' @importFrom utils read.csv write.csv
#' @importFrom pracma trapz
#' @importFrom tools file_path_sans_ext
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth ggtitle xlab ylab theme_bw theme ggsave
#' @importFrom grDevices png dev.off
#' @importFrom graphics abline
#' @importFrom stats lm
#'
#' @export
create_profile <- function(profile_csv_folder, chromosomes = 1:14, gene_position = NULL, output_folder = NULL){

  # Checking if the path are correct
  if (!dir.exists(profile_csv_folder)) stop("'profile_csv_folder' not found.")
  if (!is.null(output_folder) && !dir.exists(output_folder)) stop("'output_folder' not found.")
  # Chromosomes check
  if (!is.numeric(chromosomes) || any(chromosomes < 1 | chromosomes > 14 | chromosomes != as.integer(chromosomes))) {
    stop("'chromosomes' must be an integer vector between 1 and 14 (e.g., 1:14 or c(1, 3, 5)).")
  }


  # Load the CSV containing the start and end of each genes.
  if (is.null(gene_position)) {
    # Load the csv file containing 3 columns : gene name, the positions at which they start and end.
    genes_file_path <- system.file("extdata", "genes_positions.csv", package = "SWGACNV")
    df_genes <- read.csv(genes_file_path)
  }else{
    if (!file.exists(gene_position)) stop("'gene_position' file not found.")
    df_genes <- read.csv(gene_position)
  }

  # Load a list of the csv files path. Those files contain 3 columns : the chromosome, the position in the genome and the number of read.
  coverage_files <- list.files(profile_csv_folder, pattern = "\\.csv$", full.names = TRUE)# select all files ending by .csv
  # If only one file is found, duplicate it virtually to allow computation
  if (length(coverage_files) == 1) {
    original_file <- coverage_files[1]
    duplicate_file <- tempfile(fileext = ".csv")
    file.copy(original_file, duplicate_file)
    coverage_files <- c(original_file, duplicate_file)
    message("Only one reference file detected, at least one more would be needed for a better analysis.")
  }



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
      df_coverage_temp <- subset(df_coverage_temp, toupper(seqnames) == toupper(seqname_value))

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
    ################## Profile calculation ##################
    # Correcting sample bias.
    GlobalMeans <- colMeans(auc_profile[, 2:ncol(auc_profile)]) # Mean of each column(sample)
    Correction_Factors <- mean(GlobalMeans) / GlobalMeans # Standardization
    # Applying the correction to each column.
    for (col in names(auc_profile)[2:ncol(auc_profile)]) {
      auc_profile[[col]] <- auc_profile[[col]] * Correction_Factors[[col]] # Standardization is needed because of the variable depth
    }

    # Stock the AUC means for each chr for later standardisation in cnv_analysis()
    if (!exists("mean_auc_all_chr")) mean_auc_all_chr <- list()
    mean_auc_all_chr[[paste0("chr", sprintf("%02d", i))]] <- mean(GlobalMeans)

    ### Comparing to profile.

    gene_names <- auc_profile$gene
    n_samples <- ncol(auc_profile) - 1 # Separate the genes names from the samples values.
    ratio_list <- list()

    # Loop on each sample column.
    for (l in 2:ncol(auc_profile)) {
      vector_sample <- auc_profile[[l]]

      # Ratios (y / x for each pair of genes in the sample).
      ratios_sample <- outer(vector_sample, vector_sample, FUN = function(x, y) y / x)

      # Add the names.
      rownames(ratios_sample) <- gene_names
      colnames(ratios_sample) <- gene_names

      # List of matrix.
      ratio_list[[l - 1]] <- ratios_sample
    }

    # Convert in 3D array : [genes, genes, n_samples]
    array_ratios <- simplify2array(ratio_list)

    ## The goal is to find the 10 more stable genes to use as a reference for cnv analysis.
    # Calculation of the variance for each pair of gene between all samples.
    ratio_var <- apply(array_ratios, c(1, 2), var, na.rm = TRUE)
    rownames(ratio_var) <- gene_names
    colnames(ratio_var) <- gene_names
    # Mean of the var for each row (using the rowe instead of col because we will divide by these gene later so we need stability when dividing).
    mean_var_per_gene_row <- rowMeans(ratio_var, na.rm = TRUE)
    # 10 genes with the lowest var.
    top10_genes_row <- names(sort(mean_var_per_gene_row))[1:10]

    # Mean of the AUC ratios, will be used for comparison against new samples.
    ratio_mean <- apply(array_ratios, c(1, 2), mean, na.rm = TRUE)
    # Keep only the 10 most stable genes for comparison.
    ratio_filtered <- ratio_mean[top10_genes_row, , drop = FALSE]


    # Saving the output file
    if (is.null(output_folder)) {
      output_file <- file.path(getwd(), paste0("profile_chr", sprintf("%02d", i), ".csv"))  # Saving in the wd by default
    } else {
      if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
      output_file <- file.path(output_folder, paste0("profile_chr", sprintf("%02d", i), ".csv"))
    }

    write.csv(ratio_filtered, output_file, row.names = TRUE)
    cat("Filed saved as :", output_file, "\n")
  }


  ### CORRELATION ###
  auc_values <- auc_profile[, -1]  # Withdraw the column "gene"
  sample_names <- colnames(auc_values)

  # Double loop for each pair of samples
  for (s1 in 1:(ncol(auc_values) - 1)) {
    for (s2 in (s1 + 1):ncol(auc_values)) {
      x <- auc_values[[s1]]
      y <- auc_values[[s2]]
      cor_val <- cor(x, y, method = "pearson", use = "complete.obs")

      df_plot <- data.frame(x = x, y = y)

      plot_file <- file.path(output_folder, paste0("correlation_", sample_names[s1], "_vs_", sample_names[s2], "_chr", sprintf("%02d", i), ".png"))

      p <- ggplot(df_plot, aes(x = x, y = y)) +
        geom_point(color = "darkgreen", size = 2) +
        geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
        ggtitle(paste0(sample_names[s1], " vs ", sample_names[s2],
                       "\nPearson correlation : ", round(cor_val, 3))) +
        xlab(sample_names[s1]) +
        ylab(sample_names[s2]) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))

      ggsave(plot_file, plot = p, width = 10, height = 6)
      cat("Plot enregistré :", plot_file, "\n")
    }
  }
  ###################

  #return(invisible(mean_auc_all_chr))
  return(mean_auc_all_chr)

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
