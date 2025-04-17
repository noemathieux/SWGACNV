utils::globalVariables(c(".", "Mean_Profile", "which_label", "gene", "seqnames", "type", "mean_ratio", "color", "label"))

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("SWGACNV package correctly loaded.")
}

