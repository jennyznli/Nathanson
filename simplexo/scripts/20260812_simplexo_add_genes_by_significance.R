# ========================
# SIMPLEXO
# M1 (PATHOGENIC) ADD-ONLY GENE RESULTS, ORDERED BY SIGNIFICANCE
# ========================
library(here)
library(data.table)

setwd(here())

reports_dir <- here("simplexo", "reports", "reports2")

in_files <- list.files(
    reports_dir,
    pattern = "^SIMPLEXO4_all_array2_.*_gene_based_results_ADDonly_20260811\\.csv$",
    full.names = TRUE
)
in_files

for (in_file in in_files) {
    stem <- basename(in_file) |>
        sub("^SIMPLEXO4_all_(.*)_gene_based_results_ADDonly_20260811\\.csv$", "\\1", x = _)

    res <- fread(in_file)
    res2 <- res |>
        filter(!is.na(LOG10P_ADJ), MASK %in% c("M1", "M3")) |>
        arrange(desc(LOG10P_ADJ))

    out_file <- file.path(reports_dir, paste0("SIMPLEXO4_all_", stem, "_M1M3_ADDonly_by_significance.xlsx"))
    write_xlsx(res2, out_file)
}

in_file <- file.path(reports_dir, "SIMPLEXO4_all_array2_gene_based_results_ADDonly_20260811.csv")

res <- fread(in_file)
res2 <- res |> filter(!is.na(LOG10P_ADJ), MASK %in% c("M1", "M3")) |>
    arrange(desc(LOG10P_ADJ))

write_xlsx(res2, file.path(reports_dir, "SIMPLEXO4_all_array2_M1M3_ADDonly_by_significance.xlsx"))



