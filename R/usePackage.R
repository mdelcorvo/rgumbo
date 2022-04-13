#' usePackage
#'
#' Function to automatically check, install from both CRAN and Bioconductor repositories (if necessary) and load required packages
#'
#' \code{usePackage(pkgs)}.
#'
#' @param pkgs Numeric vector (or single value) with the name of the packages
#' @return a messagge explaining if all packages are installed.
#'

#' @examples
#' usePackage(pkgs = c('data.table','purrr'))
#' @name usePackage
NULL

#' @export
#' @rdname usePackage
"%btwn%" <- function(x, rng) {
  between(x, rng, inclusive = TRUE)
}

#' @export
#' @rdname usePackage
usePackage <- function(pkgs) {
    
    ## Install BiocManager package
    if (!requireNamespace("BiocManager", quietly = TRUE))
		install.packages("BiocManager",repos = "http://cran.rstudio.com/",lib = .libPaths()[1])
		
    ## Install packages from CRAN or Bioconductor repository                      
    isInstalled <- pkgs %in% installed.packages(lib.loc = .libPaths()[1])[, 1]
    BiocManager::install(pkgs[!isInstalled],lib = .libPaths()[1],update = FALSE,dependencies = TRUE)                                                                
    
    pkg.load <- lapply(pkgs, FUN = function(x) {
        x[!(x %in% installed.packages(.libPaths()[1])[, "Package"])]
    })
    
    if (length(unlist(pkg.load)) == 0) {
        cat("All required packages are installed \n")
    } else {
        cat(unlist(pkg.load), ": failed to install")
    }
    ## Load all packages
    lapply(pkgs, require, character.only = TRUE)
}
