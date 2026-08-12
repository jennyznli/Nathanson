# ========================
# GREY - INITIAL QC
# ========================
library(here)
setwd(here("grey"))
source(here("R", "load_packages.R"))

# ========================
# LOAD ALL DATA FILES
# ========================
ss <- read_excel(here("grey", "ss", "20260810_grey_master.xlsx"))

# ========================
# ADD MEAN TARGET COVERAGE
# ========================
# metrics1 <- read_csv(here("grey", "data", "029-Gray_PreCancer.Target-HYB-PreCancer.20240129.metrics_summary.csv"))
# metrics2 <- read_excel(here("grey", "data", "029-Gray_PreCancer.Target-HYB-PreCancer.20240923.metrics_summary.xlsx"),
#                        sheet = "Sheet1")
#
# metrics1 <- metrics1 |> select(SAMPLE, MEAN_TARGET_COVERAGE)
# metrics2 <- metrics2 |> filter(!is.na(SAMPLE)) |> select(SAMPLE, MEAN_TARGET_COVERAGE)
#
# ss <- ss |>
#     left_join(metrics1, by = c("Sample_ID" = "SAMPLE")) |>
#     left_join(metrics2, by = c("Sample_ID" = "SAMPLE"), suffix = c("_1", "_2")) |>
#     mutate(Mean_Target_Coverage = coalesce(MEAN_TARGET_COVERAGE_1, MEAN_TARGET_COVERAGE_2)) |>
#     select(-MEAN_TARGET_COVERAGE_1, -MEAN_TARGET_COVERAGE_2)
#
# # sanity check: every sample should have a coverage value, and no sample should
# # match both metrics files (they cover disjoint sequencing batches)
# stopifnot(all(!is.na(ss$Mean_Target_Coverage)))
#
# write_xlsx(ss, here("grey", "ss", "20260810_grey_master_clean.xlsx"))
# length(unique(ss$Donor_ID))
# # 116 unique donors
# ========================
# ADD MEAN TARGET COVERAGE
# ========================
table(ss$Case)
# 0  1
# 53 66

table(ss$Gene_Mutation)
# ? BRCA1 BRCA2  TP53
# 2    26    37     1

summary(ss$Mean_Target_Coverage)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1018    2254    2462    2521    2754    4245







