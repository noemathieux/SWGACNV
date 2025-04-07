
#' @title CNV identification for SWGA sequencing.
#' @description Performs CNV detection by comparing new SWGA samples to a reference profile created using `create_profile()`.
#'
#' @param csv_folder Path to the folder containing sample CSV files with columns: `seqnames`, `pos`, and `count`.
#' @param region Name(s) of the reference region(s) to compare against. Either a string (e.g., `"BENIN"`) or a vector of strings (e.g., `c("BENIN", "TOGO")`).
#' @param chr ID of the chromosome to analyze (e.g., `"Pf3D7_01_v3"`, or `"Pf3D7_02_v3"`..., or `"Pf3D7_14_v3"`).
#' @param profile_folder (optional) Path to the folder which contains a profile. This is useful if you want analyze your own sample against your own profile.
#' @param output_folder (optional)  Path to the folder where the profile files will be saved. Defaults to the working directory if not provided.
#'
#' @return A dataframe and a plot containing the score for each gene of each sample.
#'
#' @import dplyr
#' @importFrom ggplot2 ggplot aes geom_point geom_hline scale_color_manual labs theme_minimal theme guides element_text ggsave scale_x_discrete geom_abline
#' @importFrom ggrepel geom_text_repel
#' @importFrom utils read.csv write.csv
#' @importFrom purrr map_dfc
#' @importFrom pracma trapz
#' @importFrom tools file_path_sans_ext
#' @importFrom stats cor
#'
#' @export

cnv_analysis <- function(csv_folder, region, chr, profile_folder = NULL, output_folder = NULL) {
  ## test commit
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

  ## Profile selection
  # Build the profile name.
  profile_filename <- paste0("profile", chr_number, ".csv")

  if (is.null(profile_folder)) {
    profile_path <- system.file("extdata", profile_filename, package = "SWGACNV")  # Modifier le dossier si besoin
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

  profiles <- read.csv(profile_path)
  # Verify if the region exist in the columns.
  if (!all(region %in% colnames(profiles))) {
    stop(
      paste(
        "Error : Some columns specified in 'region' are not avalaible as a profile. Look for synthax error."
      )
    )
  }
  # Select the requested region columns and the genes.
  selected_region <- data.frame(profiles$gene, profiles[, region, drop = FALSE])
  colnames(selected_region)[1] <- "gene" # Rename the column "profiles$gene" in "gene"






  # Load the CSV containing the start and end of each genes.
  genes_file_path <- system.file("extdata", "genes_positions.csv", package = "SWGACNV")
  df_genes <- read.csv(genes_file_path)

  # Selecting the wanted chromosome.
  chr_prefix <- sub("_V3$", "", chr)
  df_genes <- df_genes[grep(paste0("^", chr_prefix), df_genes$gene), ]
  selected_region <- selected_region[grep(paste0("^", chr_prefix), selected_region$gene), ]

  # Load a list of each sample CSV path.
  coverage_files <- list.files(csv_folder, pattern = "\\.csv$", full.names = TRUE)
  auc_results <- data.frame(gene = df_genes$gene)

  # Calculate the AUC for each gene of each sample.
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

  Mean_Profil_Global <- colMeans(selected_region[, -1, drop = FALSE], na.rm = TRUE)# Profile mean

  # Extract the "gene" column
  CNV_results <- auc_results %>% select(gene)

  # List the regions names
  regions <- names(Mean_Profil_Global)

  # Stock the results
  results_list <- list()

  # Loop on each names of Mean_Profil_Global
  for (moy_reg in regions) {
    transformed_results <- auc_results %>%
      select(-gene) %>%
      map_dfc( ~ {
        Mean_Sample_Global <- mean(.x, na.rm = TRUE)  # Sample mean
        Correction_Factor <- Mean_Profil_Global[moy_reg] / Mean_Sample_Global  # Use the current region value

        # New values calculation
        Corrige <- .x * Correction_Factor
        Ratio <- Corrige / selected_region[[moy_reg]]
        Mean_Ratio <- mean(Ratio, na.rm = TRUE)
        SD_Ratio <- sd(Ratio, na.rm = TRUE)
        Z_Ratio <- (Ratio - Mean_Ratio) / SD_Ratio

        # Retourner un dataframe avec les 3 colonnes
        return(data.frame(Corrige, Ratio, Z_Ratio))
      })

    # Add the sample's name in the column
    colnames(transformed_results) <- paste0(rep(colnames(auc_results)[-1], each = 3),
                                            #PB ????
                                            "_",
                                            moy_reg,
                                            "_",
                                            rep(c("Adjusted", "Ratio", "z-score"), times = length(colnames(auc_results)[-1])))

    # Stocker dans une liste
    results_list[[moy_reg]] <- transformed_results

    # Fuse the results with the genes
    CNV_results <- bind_cols(auc_results %>% select(gene), transformed_results)

    # Loop on each sample to create a plot
    for (sample in colnames(auc_results)[-1]) {
      # Build the columns for this sample
      ratio_col <- paste0(sample, "_", moy_reg, "_Ratio")
      zscore_col <- paste0(sample, "_", moy_reg, "_z-score")

      # Verify if the columns exist
      if (!(ratio_col %in% colnames(CNV_results)) |
          !(zscore_col %in% colnames(CNV_results))) {
        next  # Sauter l'echantillon s'il manque des colonnes
      }

      # Output path.
      if (is.null(output_folder)) {
        output_file <- file.path(getwd(), paste0("CNVresults_", sample, chr_number, ".csv"))  # Using working directory by default.
      } else {
        if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
        output_file <- file.path(output_folder, paste0("CNVresults_", sample, chr_number, ".csv"))  # Using output path.
      }

      # Path to saving the plots
      if (is.null(output_folder)) {
        plot_path <- file.path(getwd(), paste0(sample, chr_number, "_", moy_reg, ".jpg"))  # Using working directory by default.
        plot_path2 <- file.path(getwd(), paste0(sample, " vs ", moy_reg, chr_number, "_", ".jpg"))
      } else {
        plot_path <- file.path(output_folder, paste0(sample, chr_number, "_", moy_reg, ".jpg"))  # Using output path.
        plot_path2 <- file.path(output_folder, paste0(sample, " vs ", moy_reg, chr_number, "_", ".jpg"))

      }
      if (file.exists(plot_path)) {
        file.remove(plot_path)  # Delete the file if already existing
      }
      if (file.exists(plot_path2)) {
        file.remove(plot_path2)  # Delete the file if already existing
      }

      # Plot showing expression level
      p <- ggplot(CNV_results, aes(x = gene, y = .data[[ratio_col]])) +
        geom_point(aes(color = (.data[[zscore_col]] > 2 |
                                  .data[[zscore_col]] < -2)), size = 1.5) +  # Taille des points reduite
        geom_hline(
          yintercept = 1,
          linetype = "dashed",
          color = "gray50"
        ) +

        # Display the gene name and the score if it is significant.
        geom_text_repel(
          aes(label = ifelse((.data[[zscore_col]] > 2 |
                                .data[[zscore_col]] < -2), paste0(gene, "\nZ=", round(.data[[zscore_col]], 2)), ""
          )),
          color = "red",
          size = 4,
          max.overlaps = Inf,
          force = 5,
          nudge_y = 0.2
        ) +
        scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
        labs(
          title = paste("CNV -", sample, "-", moy_reg),
          x = "Genes",
          y = "Sequencing ratio (sample/profil)"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_x_discrete(breaks = CNV_results$gene[seq(1, nrow(CNV_results), by = 25)]) +
        guides(color = "none")

      # Plot showing if the sample and the profiles are correlated
      correlation_value <- cor(
        CNV_results[[2]],
        selected_region[[moy_reg]],
        method = "pearson"
      )


      p2 <- ggplot(CNV_results, aes(x = .data[[colnames(CNV_results)[2]]])) +
        geom_point(aes(y = selected_region[[moy_reg]]), color = "darkgreen", size = 2) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 1) +
        labs(
          title = paste(sample, " vs ", moy_reg, "Pearson correlation : ", round(correlation_value, 3)),
          x = colnames(CNV_results)[2],
          y = paste("Profil moyen - ", moy_reg)
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

      print(CNV_results[[2]])
      print(selected_region[[moy_reg]])

      # Saving the plots
      ggsave(
        plot_path,
        plot = p,
        width = 12,
        height = 6,
        dpi = 300
      )
      ggsave(
        plot_path2,
        plot = p2,
        width = 12,
        height = 6,
        dpi = 300
      )
      # Saving the output file
      write.csv(CNV_results, output_file, row.names = FALSE)
      cat("Filed saved as :", output_file, "\n")
    }
  }


  # Clean the environment
  rm(
    chr,
    chr_number,
    profile_filename,
    profile_path,
    profiles,
    selected_region,
    genes_file_path,
    df_genes,
    chr_prefix,
    coverage_files,
    auc_results,
    output_file,
    cov_file,
    auc_new_col,
    df_coverage_temp,
    i,
    idx,
    auc_value,
    gene_length,
    Mean_Profil_Global,
    CNV_results,
    regions,
    results_list,
    moy_reg,
    transformed_results,
    sample,
    ratio_col,
    zscore_col,
    p,
    plot_path
  )
  gc()

}
