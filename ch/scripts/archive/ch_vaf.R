# ============================================================
# VAF & clonal characteristics — carrier vs non-carrier
# (restricted to CHIP-positive individuals)
# ============================================================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))

library(coin)

# ============================================================
# READ DATA
# ============================================================
cov <- read_excel(file.path(OUT_DIR, "pmbb_brca12_cov_chip_df.xlsx")) %>% filter(Strata == 1)
dim(cov)

vars <- read_excel(file.path(OUT_DIR, "ch_seq_wl_art_minad4_vars.xlsx"))

cov <- cov %>%
    mutate(
        CHIP_Binary = person_id %in% vars$Sample.ID,
        CHIP_Count  = sapply(person_id, function(id) sum(vars$Sample.ID == id))
    )

cat("Overall cohort: N =", nrow(cov),
    "| CHIP+ =", sum(cov$CHIP_Binary),
    sprintf("(%.1f%%)\n", 100 * mean(cov$CHIP_Binary)))
# Overall cohort: N = 3004 | CHIP+ = 193 (6.4%)

VAF_COL <- "Sample.AltFrac"   # <-- update if needed

# ============================================================
# RESTRICT TO CHIP+ INDIVIDUALS
# ============================================================
chip_pos <- cov %>% filter(CHIP_Binary == TRUE)
cat("CHIP+ individuals:", nrow(chip_pos), "\n")
cat("  Carriers:     ", sum(chip_pos$BRCA12_Case == 1), "\n")
cat("  Non-carriers: ", sum(chip_pos$BRCA12_Case == 0), "\n")

# Join variant data to CHIP+ individuals
vars_chip <- vars %>%
    inner_join(cov %>% dplyr::select(person_id, BRCA12_Case, BRCA1_Case,
                                     BRCA2_Case, Carrier),
               by = c("Sample.ID" = "person_id"))

# Per-person summary: max VAF and total mutation count
person_summary <- vars_chip %>%
    group_by(Sample.ID, BRCA12_Case, BRCA1_Case, BRCA2_Case, Carrier, Sample_age) %>%
    summarise(
        max_vaf  = max(.data[[VAF_COL]], na.rm = TRUE),
        n_muts   = n(),
        .groups  = "drop"
    )

