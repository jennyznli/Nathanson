# ========================
# CONTROL MATCHING
# ========================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
source(here("ch", "scripts", "brca_qc2.R"))

# ========================
# READ DATA
# ========================
b1_case <- read.csv(file.path("ch", "data", "pmbb_brca1_case.csv"))$x
b2_case <- read.csv(file.path("ch", "data", "pmbb_brca2_case.csv"))$x

COV <- here("PMBB", "3.0", "PMBB-Release-2024-3.0_covariates.txt")
PER <- here("PMBB", "3.0", "PMBB-Release-2024-3.0_phenotype_person.txt")
cov <- fread(COV, header = TRUE)
person <- fread(PER, header = TRUE)

flags <- read_excel(file.path("PMBB", "pmbb_flag_tables.xlsx"))
consent <- flags %>% filter(PMBB_Consent == 1)
consent_ids <- consent$PMBB_ID
length(consent_ids)
# 53062

up <- read.csv(here("simplexo", "data", "simplexo_up_map.csv"))
up$SampNum <- as.numeric(up$SampNum)

# ========================
# ADD CASE/CONTROL FLAGS
# ========================
# Add flag columns to cov df indicating BRCA1/BRCA2 case status
cov <- cov %>%
    mutate(
        BRCA1_Case = ifelse(person_id %in% b1_case$PMBB_ID, 1, 0),
        BRCA2_Case = ifelse(person_id %in% b2_case$PMBB_ID, 1, 0),
        BRCA1_Control = ifelse(!person_id %in% b1_case$PMBB_ID, 1, 0),
        BRCA2_Control = ifelse(!person_id %in% b2_case$PMBB_ID, 1, 0)
    )

# Check case counts
table(cov$BRCA1_Case)
table(cov$BRCA2_Case)

# ========================
# ADD CASE/CONTROL FLAGS
# ========================
# Overall cohort characteristics
table(carrier_status, CH_status)
summary(age ~ carrier_status)
summary(sex ~ carrier_status)

# Stratify by cancer status
table(carrier_status, cancer_diagnosis)
table(carrier_status, treatment_status)

library(MatchIt)

# Match each carrier to 1-4 non-carriers
matched <- matchit(carrier_status ~ Sample_age + Sequneced_gender + PC1 + PC2 + PC3
                   + Batch + Smoke_History,
                   data = cohort,
                   method = "nearest",
                   distance = "mahalanobis",
                   caliper = c(age = 2),  # +/- 2 years
                   ratio = 4)  # 1:4 matching

matched_data <- match.data(matched)

# Check balance
summary(matched)

# Bin age into decades
cohort$age_decade <- cut(cohort$age, breaks = seq(20, 90, 10))

# Calculate distribution in cases
case_dist <- cohort %>%
    filter(carrier_status == 1) %>%
    count(age_decade, sex, batch) %>%
    mutate(prop = n/sum(n))

# Sample controls proportionally
controls_matched <- cohort %>%
    filter(carrier_status == 0) %>%
    group_by(age_decade, sex, batch) %>%
    sample_n(size = case_count_per_bin * 4, replace = FALSE) %>%
    ungroup()

final_cohort <- bind_rows(cases, controls_matched)

# check after
# Compare distributions after matching
library(tableone)

vars <- c("age", "sex", "PC1", "PC2", "PC3", "batch")
tab <- CreateTableOne(vars = vars,
                      strata = "carrier_status",
                      data = matched_data,
                      test = TRUE)
print(tab, smd = TRUE)  # SMD < 0.1 is well-balanced

# Overall CH prevalence
matched_data %>%
    group_by(carrier_status) %>%
    summarise(
        n = n(),
        CH_n = sum(CH_status),
        CH_prev = mean(CH_status),
        median_age = median(age)
    )

# Stratified by age decade and sex
prev_table <- matched_data %>%
    group_by(carrier_status, age_decade, sex) %>%
    summarise(
        n = n(),
        CH_prev = mean(CH_status),
        .groups = "drop"
    )

# Visualize
ggplot(prev_table, aes(x = age_decade, y = CH_prev,
                       color = factor(carrier_status), group = carrier_status)) +
    geom_line() + geom_point() +
    facet_wrap(~sex) +
    labs(y = "CH Prevalence", x = "Age Decade")


# No cancer subgroup
no_cancer <- matched_data %>% filter(cancer_diagnosis == 0)

# With cancer subgroup
cancer <- matched_data %>% filter(cancer_diagnosis == 1)

# Further stratify by treatment
treated <- cancer %>% filter(treatment_status == 1)
untreated <- cancer %>% filter(treatment_status == 0)

# Run prevalence for each subgroup


### MULTIVARIATE
library(broom)

model1 <- glm(CH_status ~ carrier_status + age + sex +
                  PC1 + PC2 + PC3 + batch,
              data = matched_data,
              family = binomial())

# Get OR and 95% CI
results <- tidy(model1, exponentiate = TRUE, conf.int = TRUE)
results %>% filter(term == "carrier_status")

## stratified
# No cancer
model_no_cancer <- glm(CH_status ~ carrier_status + age + sex + PC1 + PC2 + PC3,
                       data = no_cancer,
                       family = binomial())

# Cancer, no treatment
model_untreated <- glm(CH_status ~ carrier_status + age + sex + PC1 + PC2 + PC3,
                       data = untreated,
                       family = binomial())

# Cancer, treated
model_treated <- glm(CH_status ~ carrier_status + age + sex + PC1 + PC2 + PC3,
                     data = treated,
                     family = binomial())

# Extract all results
all_models <- bind_rows(
    tidy(model1, conf.int = TRUE, exponentiate = TRUE) %>% mutate(model = "Overall"),
    tidy(model_no_cancer, conf.int = TRUE, exponentiate = TRUE) %>% mutate(model = "No cancer"),
    tidy(model_untreated, conf.int = TRUE, exponentiate = TRUE) %>% mutate(model = "Untreated"),
    tidy(model_treated, conf.int = TRUE, exponentiate = TRUE) %>% mutate(model = "Treated")
) %>%
    filter(term == "carrier_status")


CH_positive <- matched_data %>% filter(CH_status == 1)

# Number of clones
wilcox.test(n_clones ~ carrier_status, data = CH_positive)

# Median VAF
wilcox.test(median_VAF ~ carrier_status, data = CH_positive)

# Box plots
ggplot(CH_positive, aes(x = factor(carrier_status), y = median_VAF)) +
    geom_boxplot() + geom_jitter(width = 0.1, alpha = 0.3)

# Negative binomial for clone count
library(MASS)
nb_model <- glm.nb(n_clones ~ carrier_status + age + sex, data = CH_positive)

#### gene distribution analysis

# Long format with one row per mutation
mutations_long <- CH_positive %>%
    unnest(driver_genes)  # assuming you have this info

# Top genes in each group
gene_dist <- mutations_long %>%
    group_by(carrier_status, gene) %>%
    summarise(n = n()) %>%
    arrange(carrier_status, desc(n))

# Fisher's exact for specific genes (e.g., DNMT3A, TET2)
table(mutations_long$gene == "DNMT3A", mutations_long$carrier_status)
fisher.test(table(mutations_long$gene == "DNMT3A", mutations_long$carrier_status))

# Heatmap of gene frequencies
gene_matrix <- mutations_long %>%
    count(carrier_status, gene) %>%
    pivot_wider(names_from = carrier_status, values_from = n, values_fill = 0)

# Exclude prevalent cancers at time of blood draw (if timing data available)
# Different VAF thresholds for CH definition
# Include smoking status, BMI if available
# Test interaction term: carrier_status * cancer_diagnosis

