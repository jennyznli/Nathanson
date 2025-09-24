# ============================================================
# Install and load all required CRAN and Bioconductor packages
# ============================================================

cat("Loading required packages...\n")

set.seed(123)
DATE <- format(Sys.Date(), "%Y%m%d")
TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")

packages <- list(
    cran = c(
        "dplyr", "ggplot2", "readr", "tidyr", "plotly", "readxl", "here", "stringr",
        "Rtsne", "pheatmap", "RColorBrewer", "ggpubr", "ggrepel",  "lubridate", "data.table",
         "tidyverse",  "pals",  "impute"
    ),
    bioc = c(
        "BiocParallel"
    )
)

load_packages <- function(pkg_list) {
    missing_cran <- pkg_list$cran[!sapply(pkg_list$cran, requireNamespace, quietly = TRUE)]
    missing_bioc <- pkg_list$bioc[!sapply(pkg_list$bioc, requireNamespace, quietly = TRUE)]

    if (length(missing_cran) > 0) {
        install.packages(missing_cran, quiet = TRUE)
        cat("Installed missing CRAN packages: ", paste(missing_cran, collapse = ", "), "\n")
    }

    if (length(missing_bioc) > 0) {
        BiocManager::install(missing_bioc, ask = FALSE, update = FALSE, quiet = TRUE)
        cat("Installed missing Bioconductor packages: ", paste(missing_bioc, collapse = ", "), "\n")
    }

    all_pkgs <- c(pkg_list$cran, pkg_list$bioc)
    failed_pkgs <- c()

    invisible(sapply(all_pkgs, function(x) {
        if (!requireNamespace(x, quietly = TRUE)) {
            failed_pkgs <<- c(failed_pkgs, x)
        } else {
            suppressPackageStartupMessages(library(x, character.only = TRUE))
        }
    }))

    if (length(failed_pkgs) > 0) {
        cat("Failed to load the following packages: ", paste(failed_pkgs, collapse = ", "), "\n")
    } else {
        message("All packages loaded successfully!")
    }
}

load_packages(packages)
