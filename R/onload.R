.onLoad <- function(libname, pkgname) {
    # Load the dataset into the global environment
    data("hg38Tohg19", package = pkgname, envir = .GlobalEnv)
    data("hg19Tohg38", package = pkgname, envir = .GlobalEnv)
}
