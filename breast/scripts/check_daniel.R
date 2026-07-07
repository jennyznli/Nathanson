# ========================
# Packages
# ========================
packages <- c("tidyr", "dplyr", "plotly", "readr", "readxl", "here",
              "stringr", "ggplot2",  "impute", "pals", "geneplotter", "plotly"
)
purrr::walk(packages, ~ require(.x, character.only = TRUE))
here()

DATE <- format(Sys.Date(), "%Y%m%d")

# ============================================================
# FUNCTIONS
# ============================================================

icd <- read.table(file.path("PMBB", "2.0", "Phenotype", "PMBB-Release-2020-2.3_phenotype_icd-10-matrix.txt"), header = FALSE)

icd10[1:10, 1:10]
