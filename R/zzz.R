utils::globalVariables(c(".", "Mean_Profile", "which_label", "gene", "seqnames"))

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("SWGACNV package correctly loaded.")
}

