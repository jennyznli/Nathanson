# ========================
# GRAY - GERMLINE WT GENE CHECK
#
# QC-filters post-VEP germline calls for a set of genes of interest
# ========================
library(here)
setwd(here("gray"))
source(here("R", "config.R"))

# ========================
# CONFIG
# ========================
ss <- read_excel(here("gray", "ss", "pmbb_1099_breastcancercases_wildtypesamples_pmbb_20260623_pmbbid_jenny.xlsx"))
cov <- fread(file.path(here("PMBB", "4.0", "PMBB-Release-2026-4.0_phenotype_covariates.txt")), header = TRUE)

dim(ss)
dim(cov)

sum(ss$PMBB_ID %in% cov$person_id)
# all in we good

# ========================
# LOOP THROUGH + CASE
# ========================
export_filtered_for_dana <- function(filepath, out_dir, ids) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    files <- list.files(filepath, pattern = "_qc\\.csv$")

    for (f in files) {
        file <- read.csv(file.path(filepath, f), header = TRUE)
        file_filtered <- file %>% filter(Sample.ID %in% ids)
        out_name <- sub("\\.csv$", "_gray.csv", f)
        write.csv(file_filtered, file.path(out_dir, out_name), row.names = FALSE)
    }
}

export_filtered_for_dana(
    filepath = here("gray", "data", "qc"),
    out_dir = here("gray", "data", "dana"),
    ids = ss$PMBB_ID
)

add_flags <- function(df, filepath, id_col, pattern = "_qc\\.csv$") {
    files <- list.files(filepath, pattern = pattern)
    # if it doesn't have a VUS col already the df , add it and initialize as 0
    if (!"VUS" %in% names(df)) {
        df$VUS <- 0
    }
    # if it does already just  through each variant in each file and mark as 1 if the variant is VUS in the CLINVAR.SIG col
    # so that means looping through not just getting a vector of every ID in the file
    for (f in files) {
        file <- read.csv(file.path(filepath, f), header = TRUE)
        ids <- (unique(file$Sample.ID))
        gene <- unique(file$Gene)
        if (length(gene) != 1) {
            warning("Skipping ", f, ": expected exactly one gene, got ", length(gene))
            next
        }
        col_name <- paste0(gene, "_Carrier")
        df[[col_name]] <- ifelse(is.na(df[[id_col]]), NA, ifelse(df[[id_col]] %in% ids, 1, 0))
        print(paste0(gene, " Cases: ", sum(df[[col_name]], na.rm = TRUE)))

        vus_ids <- character(0)
        for (i in seq_len(nrow(file))) {
            if (grepl("Uncertain_significance", file$AutoGVP[i], ignore.case = TRUE)) {
                vus_ids <- c(vus_ids, file$Sample.ID[i])
            }
        }
        vus_ids <- unique(vus_ids)
        df$VUS <- ifelse(!is.na(df[[id_col]]) & df[[id_col]] %in% vus_ids, 1, df$VUS)
    }

    return(df)
}

ss2 <- add_flags(df = ss, filepath = here("gray", "data", "dana"), id_col = "PMBB_ID", pattern = "_gray\\.csv$")
# [1] "ATM Cases: 11"
# [1] "BRCA1 Cases: 28"
# [1] "BRCA2 Cases: 43"
# [1] "CHEK2 Cases: 24"
# [1] "PALB2 Cases: 7"

# Mismatch flag: cross-check our pipeline's carrier calls (BRCA1_Carrier /
# BRCA2_Carrier) against Dana's independent calls (Dana_BRCA12, BRCA1_Case,
# BRCA2_Case). %in% is used for the 1/0 checks below because it's NA-safe
# (NA %in% 1 is FALSE, not NA) -- the BRCA1/2_Case comparisons are guarded
# with explicit !is.na() since != would otherwise propagate NA and leave
# Mismatch as NA instead of a clean 0/1 for any row missing either value.
ss2 <- ss2 %>%
    mutate(
        Mismatch = as.integer(
            (Dana_BRCA12 %in% 1 & !(BRCA1_Carrier %in% 1) & !(BRCA2_Carrier %in% 1)) |
                (Dana_BRCA12 %in% 0 & (BRCA1_Carrier %in% 1 | BRCA2_Carrier %in% 1)) |
                (!is.na(BRCA1_Carrier) & !is.na(BRCA1_Case) & BRCA1_Carrier != BRCA1_Case) |
                (!is.na(BRCA2_Carrier) & !is.na(BRCA2_Case) & BRCA2_Carrier != BRCA2_Case)
        )
    )
print(paste0("Mismatches: ", sum(ss2$Mismatch, na.rm = TRUE)))

cov2 <- add_flags(df = cov, filepath = here("gray", "data", "qc"), id_col = "person_id")
# [1] "ATM Cases: 638"
# [1] "BARD1 Cases: 151"
# [1] "BRCA1 Cases: 955"
# [1] "BRCA2 Cases: 1189"
# [1] "CHEK2 Cases: 1173"
# [1] "PALB2 Cases: 204"
# [1] "RAD51C Cases: 112"

write_xlsx(ss2, here("gray", "ss", "pmbb_1099_breastcancercases_wildtypesamples_pmbb_20260828_pmbbid_jenny.xlsx"))
write.table(cov2, here("gray", "ss", "PMBB-Release-2026-4.0_phenotype_covariates_wt_check.txt"))

sum(ss$BRCA1_BRCA2_Case, na.rm = TRUE)
# 0
sum(ss$BRCA1_Case, na.rm = TRUE)
# 28
sum(ss$BRCA2_Case, na.rm = TRUE)
# 44
sum(ss$BRCA12_Case, na.rm = TRUE)
# 72


