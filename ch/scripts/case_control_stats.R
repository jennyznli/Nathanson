# ============================================================
# Basic stats for matched cohort
# ============================================================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))

# ============================================================
# READ DATA
# ============================================================
cov <- read_excel(file.path("ch", "data", "pmbb_brca12_cov_chip_df.xlsx")) %>% filter(Sequenced_gender == "Female", Strata == 1)
dim(cov)
# 1924

vars <- read_excel(file.path("ch", "data", "ch_seq_wl_art_minad4_vars.xlsx"))
dim(vars)
# 217

# ============================================================
# BASIC COUNTS
# ============================================================
table(cov$BRCA12_Case)
# 0    1
# 1519  405

table(cov$Carrier)
# BRCA1 BRCA1+BRCA2       BRCA2 Non-carrier
# 204           1         200        1519

cov <- cov %>%
    mutate(
        CHIP_Binary = person_id %in% vars$Sample.ID,
        CHIP_Count  = sapply(person_id, function(id) sum(vars$Sample.ID == id))
    )

table(cov$CHIP_Binary)
# FALSE  TRUE
# 1811   113

# ========================
# QUICK STATS
# ========================
sink(file.path("ch", "data", "ch_female_summary_stats.log"), split = TRUE)
cat("Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# --- helper: simple count + % table for any variable ---
basic_crosstab <- function(df, var) {
    tab  <- table(df[[var]], useNA = "ifany")
    prop <- prop.table(tab) * 100
    out  <- data.frame(
        Level = names(tab),
        N     = as.integer(tab),
        Pct   = sprintf("%.1f%%", prop)
    )
    cat("\n", var, ":\n")
    print(out, row.names = FALSE)
    invisible(out)
}

chip_crosstab <- function(df, strat_var) {
    tab_n    <- table(df[[strat_var]], df$CHIP_Binary)
    tab_prop <- prop.table(tab_n, margin = 1) * 100
    combined <- matrix(
        sprintf("%d (%.1f%%)", tab_n, tab_prop),
        nrow = nrow(tab_n),
        dimnames = list(rownames(tab_n), c("No CHIP", "CHIP"))
    )
    total_n    <- colSums(tab_n)
    total_prop <- total_n / sum(total_n) * 100
    total_row  <- matrix(
        sprintf("%d (%.1f%%)", total_n, total_prop),
        nrow = 1,
        dimnames = list("Total", c("No CHIP", "CHIP"))
    )
    out <- rbind(combined, total_row)
    cat("\nCHIP by", strat_var, ":\n")
    print(out, quote = FALSE, right = FALSE)
    invisible(out)
}

# --- basic demographic counts ---
cat("\n--- Basic counts ---\n")
for (sv in c("Sequenced_gender", "Carrier", "Batch")) {
    basic_crosstab(cov, sv)
}

# --- CHIP prevalence ---
cat("\n --- Overall CHIP prevalence --- \n\n")
ov_n    <- table(cov$CHIP_Binary)
ov_prop <- prop.table(ov_n) * 100
cat(sprintf("  No CHIP : %d (%.1f%%)\n", ov_n["FALSE"], ov_prop["FALSE"]))
cat(sprintf("  CHIP    : %d (%.1f%%)\n", ov_n["TRUE"],  ov_prop["TRUE"]))

for (sv in c("BRCA12_Case", "Carrier", "Batch")) {
    chip_crosstab(cov, sv)
}

sink()


