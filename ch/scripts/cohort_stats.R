# ============================================================
# Basic stats for matched cohort
# ============================================================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
library(moments)
library(tableone)

# ============================================================
# READ DATA
# ============================================================
cov <- read.csv(file.path("ch", "data", "pmbb_brca12_cov_df.csv"), row.names = 1) %>%
    filter(Strata %in% c(1, 2))

# cov <- cov %>%
#     mutate(
#         CHIP_Binary = person_id %in% vars$Sample.ID,
#         CHIP_Count  = sapply(person_id, function(id) sum(vars$Sample.ID == id))
#     )

cov_matched <- read_excel(file.path("ch", "data", "pmbb_brca12_cov_chip_df.xlsx"))
dim(cov_matched)
# 3004 2

# ============================================================
# TABLE 1
# ============================================================
sink(file.path("ch", "data", "ch_unaff_summary_stats.log"), split = TRUE)
cat("Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# --- Define variables for Table 1 ---
vars <- c(
    "Age",              # replace with your actual age column name
    "Sequenced_gender",
    "Smoke_History",
    "Batch",
    "Strata",
    "Class",         # replace with your actual ancestry column name
    "CHIP_Binary"
)

cat_vars <- c(
    "Sequenced_gender",
    "Smoke_History",
    "Batch",
    "Strata",
    "Class",
    "CHIP_Binary"
)

# --- Generate stratified Table 1 ---
tab1 <- CreateTableOne(
    vars       = vars,
    strata     = "BRCA12_Case",    # carrier vs non-carrier
    data       = cov_matched,
    factorVars = cat_vars,
    addOverall = TRUE               # adds an "Overall" column too
)

print(tab1,
      showAllLevels = TRUE,
      quote         = FALSE,
      noSpaces      = TRUE,
      printToggle   = TRUE,
      formatOptions = list(big.mark = ","))

# --- Export to CSV ---
tab1_out <- print(tab1,
                  showAllLevels = TRUE,
                  quote         = FALSE,
                  noSpaces      = TRUE,
                  printToggle   = FALSE)

write.csv(tab1_out,
          file = file.path("ch", "data", "table1_by_carrier.csv"))

sink()

# ============================================================
# Wilcoxon - age in case vs. controls
# ============================================================
wilcox.test(Sample_age ~ BRCA12_Case, data = cov)
# Wilcoxon rank sum test with continuity correction
#
# data:  Sample_age by BRCA12_Case
# W = 15295432, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0
cov %>%
    group_by(BRCA12_Case) %>%
    summarise(
        median_age = median(Sample_age),
        q25 = quantile(Sample_age, 0.25),
        q75 = quantile(Sample_age, 0.75)
    )
# BRCA12_Case median_age   q25   q75
# <int>      <dbl> <dbl> <dbl>
# 1           0       52.7  35.9  64.4
# 2           1       40.8  32.5  51.8

wilcox.test(Sample_age ~ BRCA12_Case, data = cov_matched)
# Wilcoxon rank sum test with continuity correction
#
# data:  Sample_age by BRCA12_Case
# W = 784444, p-value = 0.8702
# alternative hypothesis: true location shift is not equal to 0

cov_matched %>%
    group_by(BRCA12_Case) %>%
    summarise(
        median_age = median(Sample_age),
        q25 = quantile(Sample_age, 0.25),
        q75 = quantile(Sample_age, 0.75)
    )
# # A tibble: 2 × 4
# BRCA12_Case median_age   q25   q75
# <dbl>      <dbl> <dbl> <dbl>
#     1           0       40.8  31.6  53.0
# 2           1       40.8  32.3  52

# ============================================================
# SHAPIRO-WILK - age carriers
# ============================================================
# Normality tests by group
shapiro.test(cov$Sample_age[cov$BRCA12_Case == 1])   # carriers
# shapiro.test(cov$Sample_age[cov$BRCA12_Case == 0])   # controls
# Shapiro-Wilk normality test
#
# data:  cov$Sample_age[cov$BRCA12_Case == 1]
# W = 0.96908, p-value = 6.957e-12

# Skewness
skewness(cov$Sample_age[cov$BRCA12_Case == 1])

# Visual confirmation
ggplot(cov, aes(x = Sample_age, fill = factor(BRCA12_Case))) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c("0" = "grey60", "1" = "#B71C1C"),
                      labels = c("Controls", "gBRCA1/2 carriers")) +
    labs(x = "Age", y = "Density", fill = NULL) +
    theme_minimal()

# ============================================================
# CHI - SQUARED for freeze
# ============================================================
chisq.test(table(cov$BRCA12_Case, cov$Batch))
# Pearson's Chi-squared test with Yates' continuity correction
#
# data:  freeze_table
# X-squared = 1793.1, df = 1, p-value < 2.2e-16

chisq.test(table(cov_matched$BRCA12_Case, cov_matched$Batch))
# Pearson's Chi-squared test with Yates' continuity correction
#
# data:  table(cov_matched$BRCA12_Case, cov_matched$Batch)
# X-squared = 3.5579, df = 1, p-value = 0.05926

# ============================================================
# CHI - SQUARED for sex
# ============================================================
chisq.test(table(cov$BRCA12_Case, cov$Sequenced_gender))

# Pearson's Chi-squared test with Yates' continuity correction
#
# data:  table(cov$BRCA12_Case, cov$Sequenced_gender)
# X-squared = 250.4, df = 1, p-value < 2.2e-16

chisq.test(table(cov_matched$BRCA12_Case, cov_matched$Sequenced_gender))
# Pearson's Chi-squared test with Yates' continuity correction
#
# data:  table(cov_matched$BRCA12_Case, cov_matched$Sequenced_gender)
# X-squared = 3.527, df = 1, p-value = 0.06038

table(cov_matched$Sequenced_gender, cov_matched$BRCA12_Case)

# ============================================================
# CHIP PREVALENCE
# ============================================================
chisq.test(table(cov$BRCA12_Case, cov$Sequenced_gender))

# Pearson's Chi-squared test with Yates' continuity correction
#
# data:  table(cov$BRCA12_Case, cov$Sequenced_gender)
# X-squared = 250.4, df = 1, p-value < 2.2e-16

chisq.test(table(cov_matched$BRCA12_Case, cov_matched$Sequenced_gender))
# Pearson's Chi-squared test with Yates' continuity correction
#
# data:  table(cov_matched$BRCA12_Case, cov_matched$Sequenced_gender)
# X-squared = 3.527, df = 1, p-value = 0.06038

table(cov_matched$Sequenced_gender, cov_matched$BRCA12_Case)



