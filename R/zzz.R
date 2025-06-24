utils::globalVariables(c(".", "Mean_Profile", "which_label", "gene", "seqnames", "type", "mean_ratio", "color", "label", "count"))

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("SWGACNV package correctly loaded.")
}

if (.Platform$OS.type == "windows" && !pkgbuild::has_build_tools()) {
  stop("Rtools is required to install this package. Please install Rtools from https://cran.r-project.org/bin/windows/Rtools/")
}

