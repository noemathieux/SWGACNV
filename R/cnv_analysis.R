
#' @title CNV identification for SWGA sequencing.
#' @description Performs CNV detection by comparing new SWGA samples to a reference profile created using `create_profile()`.
#'
#' @param csv_folder Path to the folder containing sample CSV files with columns: `seqnames`, `pos`, and `count`.
#' @param chr ID of the chromosome to analyze (e.g., `"Pf3D7_01_v3"`, or `"Pf3D7_02_v3"`..., or `"Pf3D7_14_v3"`).
#' @param mean_profile Vector of integer containing the mean AUC for each chromosome. Returned by create_profile() and used for standardization.
#' @param gene_position (optional)Path to the file containing each gene to analyze with their start and end. It should have 3 column : `gene`, `start`, `end`. Default to the one provided with the package.
#' @param profile_folder (optional) Path to the folder which contains a profile. This is useful if you want analyze your own sample against your own profile.
#' @param output_folder (optional)  Path to the folder where the profile files will be saved. Defaults to the working directory if not provided.
#'
#' @return A dataframe and a plot containing the score for each gene of each sample.
#'
#' @importFrom dplyr mutate across
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot aes geom_point geom_hline labs theme_minimal theme element_text ggsave scale_x_discrete scale_color_identity element_rect
#' @importFrom ggrepel geom_text_repel
#' @importFrom utils read.csv write.csv
#' @importFrom purrr map_dfc
#' @importFrom pracma trapz
#' @importFrom tools file_path_sans_ext
#' @importFrom stats cor var
#'
#' @export

cnv_analysis <- function(csv_folder, chr, mean_profile, profile_folder = NULL, gene_position = NULL, output_folder = NULL) {
  # Checking if files and folder exists
  if (!dir.exists(csv_folder)) stop("csv_folder not found.")
  if (!is.null(profile_folder)) {
    if (!dir.exists(profile_folder)) stop("profile_folder not found.")
  }
  chr <- toupper(chr)
  # Load the csv containing the score of each region.
  chr_number <- sub("PF3D7_", "", chr)
  chr_number <- sub("_V3", "", chr_number)
  chr_number <- sprintf("%02d", as.numeric(chr_number)) # Extract the chromosome name to load the right profile.

  # Profile selection
  # Build the profile name.
  profile_filename <- paste0("profile_chr", chr_number, ".csv")

  if (is.null(profile_folder)) {
    profile_path <- system.file("extdata", profile_filename, package = "SWGACNV")
    # Load the corresponding profile.
    if (!file.exists(profile_path)) {
      stop(paste(
        "ERROR :",
        profile_filename,
        "does not exist. The error may come from",
        chr
      ))
    }
  } else {
    profile_path <- file.path(profile_folder, profile_filename)
  }

  # Load the profile and add the row name.
  profile <- read.csv(profile_path, check.names = FALSE)
  colnames(profile)[1] <- "gene"
  rownames(profile) <- profile$gene
  profile <- profile[, -1, drop = FALSE]
  top10_genes_row <- rownames(profile)  # Extract the gene name for later filtering.



  # Load the mean AUC for standardization
  mean_profile_chr <- mean_profile[[paste0("chr", chr_number)]]

  # Load the CSV containing the start and end of each genes.
    if (is.null(gene_position)) {
    # Load the csv file containing 3 columns : gene name, the positions at which they start and end.
    genes_file_path <- system.file("extdata", "genes_positions.csv", package = "SWGACNV")
    df_genes <- read.csv(genes_file_path)
  }else{
    if (!file.exists(gene_position)) stop("'gene_position' file not found.")
    df_genes <- read.csv(gene_position)
  }

  # Selecting the wanted chromosome.
  chr_prefix <- sub("_V3$", "", chr)
  df_genes <- df_genes[grep(paste0("^", chr_prefix), df_genes$gene), ]

  # Load a list of each sample CSV path.
  coverage_files <- list.files(csv_folder, pattern = "\\.csv$", full.names = TRUE)
  auc_results <- data.frame(gene = df_genes$gene)

  # Calculate the AUC for each gene of each sample. Loop for each file.
  for (cov_file in coverage_files) {
    auc_new_col <- sub("\\.sorted\\.bam_coverage$",
                       "",
                       file_path_sans_ext(basename(cov_file)))
    auc_results[[auc_new_col]] <- NA

    cat ("Treating :", cov_file, "\n")
    df_coverage_temp <- read.csv(cov_file)
    df_coverage_temp <- subset(df_coverage_temp, toupper(seqnames) == toupper(chr))
    for (i in 1:nrow(df_genes)) {
      print(i)
      idx <- which(
        df_coverage_temp$pos >= df_genes$start[i]
        & df_coverage_temp$pos <= df_genes$end[i]
      )
      auc_value <- trapz(df_coverage_temp$pos[idx], df_coverage_temp$count[idx])
      gene_length <- (df_genes$end[i] - df_genes$start[i] + 1)
      auc_results[[auc_new_col]][i] <- auc_value / gene_length
    }
  }
  ################## New sample score ##################

    # Correcting sample bias.
    mean_sample <- colMeans(as.matrix(auc_results[, 2:ncol(auc_results)]), na.rm = TRUE) # Mean of each column(sample)
    #mean_sample <- colMeans(auc_results[, 2:ncol(auc_results)]) # Mean of each column(sample)
    Correction_Factors <- mean_profile_chr / mean_sample # Standardization
    # Applying the correction to each column.
    for (col in names(auc_results)[2:ncol(auc_results)]) {
      auc_results[[col]] <- auc_results[[col]] * Correction_Factors[[col]]
    }

    gene_names <- auc_results$gene

    # Loop on each sample to create a plot
    for (sample in colnames(auc_results)[-1]) {

      values_sample <- auc_results[[sample]]

      ratios_sample <- outer(values_sample, values_sample, FUN = function(x, y) y / x)
      rownames(ratios_sample) <- gene_names
      colnames(ratios_sample) <- gene_names
      ratio_filtered_sample <- ratios_sample[top10_genes_row, , drop = FALSE]

      # Aline the columns and row in both of the matrix.
      ratio_filtered_sample <- ratio_filtered_sample[rownames(profile), colnames(profile)]

      ratio_filtered_sample <- as.matrix(ratio_filtered_sample)
      profile <- as.matrix(profile)

      mode(ratio_filtered_sample) <- "numeric"
      mode(profile) <- "numeric"
      matrix_diff <- ratio_filtered_sample / profile

      # Output path.
      if (is.null(output_folder)) {
        output_file <- file.path(getwd(), paste0("CNVresults_", sample, chr_number, ".csv"))  # Using working directory by default.
      } else {
        if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
        output_file <- file.path(output_folder, paste0("CNVresults_", sample, "_chr", chr_number, ".csv"))  # Using output path.
      }
      write.csv(matrix_diff, output_file, row.names = TRUE)

      #Keep only the best value to eliminate error
      mean_values <- apply(matrix_diff, 2, function(vec) {
        vec_sorted <- vec[order(abs(vec - 1))]  # Sort values by their distance to 1.
        mean(vec_sorted[1:8], na.rm = TRUE)    # Mean of the 8 best.
      })

      # DF for ggplot
      df_plot <- data.frame(
        gene = names(mean_values),
        mean_ratio = mean_values
      )

      # Color swap and display name for gene with 50% more or less variation
      df_plot$color <- ifelse(df_plot$mean_ratio > 1.5 | df_plot$mean_ratio < 0.5, "red", "blue")
      df_plot$label <- ifelse(df_plot$mean_ratio > 1.5 | df_plot$mean_ratio < 0.5, df_plot$gene, NA)
      # display only 1 in 10 gene in x
      df_plot$gene <- factor(df_plot$gene, levels = df_plot$gene)  # Keep the same order so label do not get mixed up with the wrong value
      x_labels <- levels(df_plot$gene)
      x_labels[!(seq_along(x_labels) %% 10 == 1)] <- ""  # Display only 1 out of 10 labels

      # Plot x=gene names y=ratio btw profile and sample
      p <- ggplot(df_plot, aes(x = gene, y = mean_ratio, color = color)) +
        geom_point() +
        geom_text_repel(
          aes(label = label),
          size = 3,
          box.padding = 0.3,
          point.padding = 0.2,
          max.overlaps = Inf,
          segment.color = "red",
          segment.size = 0.3
        ) +
        geom_hline(yintercept = 1, color = "blue", linetype = "dashed", linewidth = 0.6) +
        geom_hline(yintercept = 1.5, color = "darkred", linetype = "dotted", linewidth = 0.6) +
        geom_hline(yintercept = 0.5, color = "darkred", linetype = "dotted", linewidth = 0.6) +
        scale_color_identity() +
        scale_x_discrete(labels = x_labels) +
        labs(
          title = paste("Mean CNV Ratio per Gene -", sample),
          x = "Gene",
          y = "Mean CNV Ratio"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA)
        )

      # Path to saving the plots
      if (is.null(output_folder)) {
        plot_path <- file.path(getwd(), paste0("MeanRatioPlot_", sample, "_chr", chr_number, ".png"))  # Using working directory by default.
      } else {
        plot_path <- file.path(output_folder, paste0("MeanRatioPlot_", sample, "_chr", chr_number, ".png"))  # Using output path.
      }
      if (file.exists(plot_path)) {
        file.remove(plot_path)  # Delete the file if already existing
      }

      # Saving the plot, supress geom_text_repel label messages
      suppressMessages(suppressWarnings(ggsave(plot_path, plot = p, width = 10, height = 6, bg = "white")))


  # # Clean the environment
  # rm(
  #   chr,
  #   chr_number,
  #   profile_filename,
  #   profile_path,
  #   profile,
  #   selected_region,
  #   genes_file_path,
  #   df_genes,
  #   chr_prefix,
  #   coverage_files,
  #   output_file,
  #   cov_file,
  #   auc_new_col,
  #   df_coverage_temp,
  #   i,
  #   idx,
  #   auc_value,
  #   gene_length,
  #   Mean_Profil_Global,
  #   CNV_results,
  #   regions,
  #   results_list,
  #   moy_reg,
  #   transformed_results,
  #   sample,
  #   ratio_col,
  #   zscore_col,
  #   p,
  #   plot_path
  # )
  #
  # gc()

}
}
