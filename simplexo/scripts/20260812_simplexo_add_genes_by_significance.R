# ========================
# SIMPLEXO
# M1 (PATHOGENIC) ADD-ONLY GENE RESULTS, ORDERED BY SIGNIFICANCE
# ========================
library(here)
library(data.table)

setwd(here())

reports_dir <- here("simplexo", "reports", "reports")

in_file <- file.path(reports_dir, "SIMPLEXO4_all_array1_gene_based_results_ADDonly_20260811.csv")

res <- fread(in_file)
res2 <- res |> filter(!is.na(LOG10P_ADJ), MASK %in% c("M1", "M3")) |>
    arrange(desc(LOG10P_ADJ))

write_xlsx(res2, file.path(reports_dir, "SIMPLEXO4_all_array1_M1M3_ADDonly_by_significance.xlsx"))
