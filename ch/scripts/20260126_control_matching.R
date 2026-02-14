# ========================
# PROPENSITY SCORE MATCHING
# BRCA1/2 Carriers vs Non-carriers for CHIP Analysis
# ========================
library(here)
library(MatchIt)
library(tableone)
library(cobalt)

setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
source(here("ch", "scripts", "brca_qc2.R"))

# ========================
# READ DATA
# ========================
b1_case <- read.csv(file.path("ch", "data", "pmbb_brca1_case.csv"))$x
b2_case <- read.csv(file.path("ch", "data", "pmbb_brca2_case.csv"))$x
b12_case <- unique(c(b1_case, b2_case))

COV <- here("PMBB", "3.0", "PMBB-Release-2024-3.0_covariates.txt")
PER <- here("PMBB", "3.0", "PMBB-Release-2024-3.0_phenotype_person.txt")

cov <- fread(COV, header = TRUE)
person <- fread(PER, header = TRUE)

# ========================
# SMOKING HISTORY
# ========================
### PMBB ###
substance <- fread(here("PMBB", "3.0", "PMBB-Release-2024-3.1_phenotype_substance_use.txt"),
                   header = TRUE)

# turn "Not Asked" into NA values
substance <- substance %>%
    mutate(TOBACCO_USER = ifelse(TOBACCO_USER == "Not Asked", NA, TOBACCO_USER))

unique(substance$TOBACCO_USER)
# "Never"   "Quit"    "Yes"     NA        "Passive"

# collapse down values in this column
merge_duplicates <- function(df, group_col) {
    df %>%
        group_by(across(all_of(group_col))) %>%
        summarise(
            across(-any_of(group_col), ~ {
                non_na_vals <- unique(na.omit(as.character(.x)))
                if(length(non_na_vals) == 0) NA_character_ else paste(non_na_vals, collapse = ";")
            }),
            .groups = "drop"
        )
}

substance_merge <- merge_duplicates(substance, "person_id")
# unique(substance_merge$TOBACCO_USER)

# non smokers will only have Never any no other values
pmbb_non_smoker <- substance_merge %>% filter(TOBACCO_USER == "Never")
pmbb_non_smoker_ids <- sort(unique(pmbb_non_smoker$person_id))
length(pmbb_non_smoker_ids)
# 28369

# smokers will have anything other than Never and can't be empty
pmbb_smoker <- substance_merge %>% filter(!is.na(TOBACCO_USER) & TOBACCO_USER != "Never")
pmbb_smoker_ids <- sort(unique(pmbb_smoker$person_id))
length(pmbb_smoker_ids)
# 25638

### PROGENY ###
progeny <- read_excel(here("ch", "ss", "brca_carriers_ch_freq_w_seen_in_crep_20251020.xlsx"),
                      sheet = "Data_from_master_table")
progeny <- progeny %>% filter(DNA == "D", !is.na(SampNum))
up <- read.csv(here("simplexo", "data", "simplexo_up_map.csv"))
progeny$SampNum <- as.numeric(progeny$SampNum)

# join with those only in PMBB key
progeny_pmbb <- progeny %>%
    left_join(up, by = "SampNum") %>%
    filter(!is.na(PMBB_ID))
dim(progeny_pmbb)
# 1278

table(progeny_pmbb$SmokingEver)
# Never      No Unknown     Yes
# 775      62      22     367

# get rid of unknowns
progeny_pmbb_smoke <- progeny_pmbb %>%
    filter(SmokingEver != "Unknown", !is.na(SmokingEver))

# non smokers have Never, No
progeny_non_smoker_ids <- progeny_pmbb_smoke %>%
    filter(SmokingEver %in% c("Never", "No")) %>%
    pull(PMBB_ID) %>%
    unique()
length(progeny_non_smoker_ids)
# 837

# smokers have yes
progeny_smoker_ids <- progeny_pmbb_smoke %>%
    filter(SmokingEver == "Yes") %>%
    pull(PMBB_ID) %>%
    unique()
length(progeny_smoker_ids)
# 367

### MERGED SMOKERS ###
cov_smoke <- cov %>% mutate(
    PMBB_SmokeHistory = case_when(
        person_id %in% pmbb_smoker_ids ~ 1,
        person_id %in% pmbb_non_smoker_ids ~ 0,
        TRUE ~ NA_real_
    ),
    Progeny_SmokeHistory = case_when(
        person_id %in% progeny_smoker_ids ~ 1,
        person_id %in% progeny_non_smoker_ids ~ 0,
        TRUE ~ NA_real_
    )
)
table(cov_smoke$PMBB_SmokeHistory, useNA = "ifany")
# 0     1  <NA>
#     28369 25638  3163
table(cov_smoke$Progeny_SmokeHistory, useNA = "ifany")
# 0     1  <NA>
# 837   367 55966

table(PMBB = cov_smoke$PMBB_SmokeHistory,
      Progeny = cov_smoke$Progeny_SmokeHistory,
      useNA = "ifany")
# Progeny
# PMBB       0     1  <NA>
#     0      507    27 27835
# 1       15   207 25416
# <NA>   315   133  2715

# further examine
# x <- cov_smoke %>% filter(Progeny_SmokeHistory == 1, PMBB_SmokeHistory == 0)
# y <- cov_smoke %>% filter(Progeny_SmokeHistory == 0, PMBB_SmokeHistory == 1)

# COMBINE
# if any have smoking history, mark as smoker
# if they are non smokers in both, mark as non smoker
all_smoker_ids <- sort(union(pmbb_smoker_ids, progeny_smoker_ids))
all_non_smoker_ids <- setdiff(union(pmbb_non_smoker_ids, progeny_non_smoker_ids), all_smoker_ids)

length(all_smoker_ids)
# 25798
length(all_non_smoker_ids)
# 28657

# ========================
# CANCER FLAGS
# ========================
cancer_df <- read.csv(file.path("ch", "data", "pmbb_cancer_age.csv"),
                      row.names = NULL)[2:4]
dim(cancer_df)
# 25362     3

# ========================
# PREPARE MATCHING DATASET
# ========================
cov1 <- cov %>%
    mutate(
        BRCA12_Case = ifelse(person_id %in% b12_case, 1, 0),
        BRCA1_Case = ifelse(person_id %in% b1_case, 1, 0),
        BRCA2_Case = ifelse(person_id %in% b2_case, 1, 0),
        Carrier = case_when(
            person_id %in% b1_case & person_id %in% b2_case ~ "BRCA1+BRCA2",
            person_id %in% b1_case ~ "BRCA1",
            person_id %in% b2_case ~ "BRCA2",
            TRUE ~ "Non-carrier"
        ),
        Smoke_History = case_when(
            person_id %in% all_smoker_ids ~ 1,
            person_id %in% all_non_smoker_ids ~ 0,
            TRUE ~ NA_real_
        ),
        Cancer = ifelse(person_id %in% cancer_df$PMBB_ID, 1, 0),
    ) %>% left_join(
        cancer_df[1:2], by = c("person_id" = "PMBB_ID")
    ) %>% mutate(
        Class = ifelse(Class %in% c("UNKNOWN1", "UNKNOWN2", "UNKNOWN3"), "UNKNOWN", Class)
    )


table(cov1$Smoke_History, useNA = "ifany")
# 0     1  <NA>
# 28657 25798  2715
table(cov1$BRCA12_Case)
table(cov1$BRCA1_Case)
table(cov1$BRCA2_Case)

# ========================
# CREATE CANCER STRATA
# ========================
cov2 <- cov1 %>%
    mutate(
        Strata = case_when(
            # Strata 1: No cancer
            Cancer == 0 ~ 1,
            # Strata 2: Cancer developed >1 year AFTER sample collection
            Cancer == 1 & !is.na(Dx_Age) & Dx_Age > (Sample_age + 1) ~ 2,
            # Strata 3: Cancer diagnosed BEFORE/AT sample collection
            Cancer == 1 & !is.na(Dx_Age) & Dx_Age <= Sample_age ~ 3,
            # Strata 4: Cancer within 1 year of sample OR missing dx age
            Cancer == 1 ~ 4,
            # NA - missing DX age
            TRUE ~ NA_real_
        ),
        Strata_Label = case_when(
            Strata == 1 ~ "No cancer",
            Strata == 2 ~ "Post-sample cancer (>1yr)",
            Strata == 3 ~ "Pre-sample cancer",
            Strata == 4 ~ "Uncertain",
            TRUE ~ NA_character_
        )
    )

# ========================
# CHECK STRATA DISTRIBUTION
# ========================
table(cov2$Strata_Label, useNA = "ifany")
# No cancer Post-sample cancer (>1yr)         Pre-sample cancer                 Uncertain
# 31808                      2419                     20889                      2054

table(cov2$Class, useNA = "ifany")
# AFR      AMR      EAS      EUR      SAS UNKNOWN0 UNKNOWN1 UNKNOWN2
# 12113      761      907    42133      743      448       59        6

table(cov2$Smoke_History, useNA = "ifany")
# 0     1  <NA>
# 28657 25798  2715

t1 <- table(Carrier = cov2$Carrier, Strata = cov2$Strata_Label, useNA = "ifany")
print(t1)
# Strata
# Carrier       No cancer Post-sample cancer (>1yr) Pre-sample cancer Uncertain
# BRCA1             351                        61               426        23
# BRCA1+BRCA2         2                         0                 7         0
# BRCA2             370                        58               433        32
# Non-carrier     31085                      2300             20023      1999

write.csv(t1, file.path("ch", "data", "strata_by_germline_table.csv"))

t2 <- table(Carrier = cov2$Carrier, Smoking = cov2$Smoke_History, useNA = "ifany")
print(t2)
# Smoking
# Carrier           0     1  <NA>
#     BRCA1         509   274    78
# BRCA1+BRCA2     6     3     0
# BRCA2         557   278    58
# Non-carrier 27585 25243  2579

write.csv(t2, file.path("ch", "data", "smoking_by_germline_table.csv"))
# Strata
# Carrier           1     2
# BRCA1         150   711
# BRCA1+BRCA2     1     8
# BRCA2         259   634
# Non-carrier 43503 11904

t3 <- table(Carrier = cov2$Carrier, Strata = cov2$Batch, useNA = "ifany")
print(t3)

cov3 <- cov2 %>% filter(!is.na(Smoke_History))
dim(cov3)
# 54455

# ========================
# UNAFFECTED ANALYSIS
# ========================
unaff <- cov3 %>%
    filter(Strata %in% c(1, 2))

dim(unaff)
# 32000    21

table(unaff$BRCA12_Case)
# 0     1
# 31252   748

table(unaff$BRCA1_Case)
# 0     1
# 31637   363

table(unaff$BRCA2_Case)
# 0     1
# 31613   387

table(unaff$Sequenced_gender)
# Female   Male
# 16552  15448

table(unaff$Batch)
# 1     2
# 26076  5924

table(unaff$Smoke_History)
# 0     1
# 17667 14333

range(unaff$Sample_age)
# 18 - 90

# ========================
# CHECK INITIAL IMBALANCE
# ========================
match_formula <- as.formula(paste(
    "BRCA12_Case ~ Sample_age + Sequenced_gender + Smoke_History + Batch +",
    paste(paste0("PC", 1:6), collapse = " + ")
))

# no matching - prematch object
m.out0 <- matchit(
    match_formula,
    data = unaff,
    method = NULL,
    distance = "glm" # generalized linear model (logistic by default), try ?distance for others
)

# check balance
summary(m.out0)

# Summary of Balance for All Data:
#     Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean eCDF Max
# distance                      0.1247        0.0210          1.1376     4.4916    0.3577   0.6197
# Sample_age                   43.3591       51.2122         -0.5861     0.6182    0.1110   0.2597
# Sequenced_genderFemale        0.7727        0.5111          0.6242          .    0.2616   0.2616
# Sequenced_genderMale          0.2273        0.4889         -0.6242          .    0.2616   0.2616
# Smoke_History                 0.2901        0.4517         -0.3560          .    0.1616   0.1616
# Batch                         1.7647        1.1713          1.3981     1.2694    0.2967   0.5935
# PC1                           0.0108        0.0015          0.8674     0.2186    0.0644   0.1994
# PC2                           0.0443        0.0334          0.8341     0.2785    0.0621   0.2251
# PC3                          -0.0174       -0.0149         -0.2585     0.4982    0.0487   0.1525
# PC4                          -0.0095       -0.0035         -0.3929     1.1482    0.1628   0.2363
# PC5                           0.0058        0.0113         -0.4902     0.9954    0.1476   0.2123
# PC6                          -0.0004        0.0002         -0.0776     0.7755    0.0153   0.0351

# ========================
# PROPENSITY SCORE MATCHING
# ========================
# try out different ratios
m.out1 <- matchit(
    match_formula,
    data = unaff,
    method = "nearest",
    distance = "glm",
    ratio = 1,
    caliper = c(0.2, # use 0.2 STD of logit(PS)
                Sample_age = 5), # +- 5 years for sample age
    std.caliper = c(TRUE, FALSE), # goes with caliper field
    exact = ~ Sequenced_gender + Batch + Smoke_History,
    replace = FALSE
)
m.out2 <- matchit(
    match_formula,
    data = unaff,
    method = "nearest",
    distance = "glm",
    ratio = 2,
    caliper = c(0.2, # use 0.2 STD of logit(PS)
                Sample_age = 5), # +- 5 years for sample age
    std.caliper = c(TRUE, FALSE), # goes with caliper field
    exact = ~ Sequenced_gender + Batch + Smoke_History,
    replace = FALSE
)
m.out3 <- matchit(
    match_formula,
    data = unaff,
    method = "nearest",
    distance = "glm",
    ratio = 3,
    caliper = c(0.2, # use 0.2 STD of logit(PS)
                Sample_age = 5), # +- 5 years for sample age
    std.caliper = c(TRUE, FALSE), # goes with caliper field
    exact = ~ Sequenced_gender + Batch + Smoke_History,
    replace = FALSE
)
m.out4 <- matchit(
    match_formula,
    data = unaff,
    method = "nearest",
    distance = "glm",
    ratio = 4,
    caliper = c(0.2, # use 0.2 STD of logit(PS)
                Sample_age = 5), # +- 5 years for sample age
    std.caliper = c(TRUE, FALSE), # goes with caliper field
    exact = ~ Sequenced_gender + Batch + Smoke_History,
    replace = FALSE
)
m.out5 <- matchit(
    match_formula,
    data = unaff,
    method = "nearest",
    distance = "glm",
    ratio = 5,
    caliper = c(0.2, # use 0.2 STD of logit(PS)
                Sample_age = 5), # +- 5 years for sample age
    std.caliper = c(TRUE, FALSE), # goes with caliper field
    exact = ~ Sequenced_gender + Batch + Smoke_History,
    replace = FALSE
)
# not all treated will get 5 matches, so we're running out ...
summary(m.out1)
summary(m.out2)
summary(m.out3)
summary(m.out4)
summary(m.out5)

# ========================
# PLOTS
# ========================
### jitter
pdf(file.path("ch", "figures", "matching_jitter_plot.pdf"), width = 8, height = 5)
plot(m.out5, type = "jitter", interactive = FALSE)
dev.off()

### density plots
pdf(file.path("ch", "figures", "matching_density_plots.pdf"), width = 8, height = 5)
plot(m.out5, type = "density", interactive = FALSE,
     which.xs = ~Sample_age + Batch + Sequenced_gender)
dev.off()

### love plot - SMD before and after matching
love_plot <- love.plot(
    m.out5,
    stats = "mean.diffs",
    threshold = c(m = 0.1),  # 0.1 threshold line
    var.order = "unadjusted",
    abs = TRUE,
    stars = "std", # those with mean diff that have been standardized
    title = "Covariate Balance: Before vs After Matching",
    colors = c("#B2182B", "#2166AC")
)
ggsave(file.path("ch", "figures", "matching_love_plot.pdf"),
       love_plot, width = 10, height = 6)

### Balance table
bal_tab <- bal.tab(
    m.out5,
    m.threshold = 0.1,  # Flag variables with SMD > 0.1
    un = TRUE,           # Show unmatched balance
    disp.v.ratio = TRUE  # Show variance ratios
)
print(bal_tab)
# Balance Measures
# Type Diff.Un V.Ratio.Un Diff.Adj    M.Threshold V.Ratio.Adj
# distance              Distance  1.1376     4.4916   0.0078 Balanced, <0.1      1.0149
# Sample_age             Contin. -0.5861     0.6182   0.0214 Balanced, <0.1      0.9495
# Sequenced_gender_Male   Binary -0.2616          .   0.0000 Balanced, <0.1           .
# Smoke_History           Binary -0.1616          .   0.0000 Balanced, <0.1           .
# Batch_2                 Binary  0.5935          .  -0.0000 Balanced, <0.1           .
# PC1                    Contin.  0.8674     0.2186  -0.0007 Balanced, <0.1      0.8466
# PC2                    Contin.  0.8341     0.2785   0.0119 Balanced, <0.1      0.8035
# PC3                    Contin. -0.2585     0.4982  -0.0107 Balanced, <0.1      0.6296
# PC4                    Contin. -0.3929     1.1482  -0.0085 Balanced, <0.1      1.8775
# PC5                    Contin. -0.4902     0.9954  -0.0109 Balanced, <0.1      1.4326
# PC6                    Contin. -0.0776     0.7755  -0.0017 Balanced, <0.1      0.8602
#
# Balance tally for mean differences
# count
# Balanced, <0.1        11
# Not Balanced, >0.1     0
#
# Variable with the greatest mean difference
# Variable Diff.Adj    M.Threshold
# Sample_age   0.0214 Balanced, <0.1
#
# Sample sizes
# Control Treated
# All                  31252.      748
# Matched (ESS)         2387.1     732
# Matched (Unweighted)  2988.      732
# Unmatched            28264.       16


matched_data <- match.data(m.out5)
dim(matched_data)
# 3720

write.csv(matched_data,
          file.path("ch", "data", "psm_matched_data.csv"),
          row.names = FALSE)

# split into case/controls
case_ids <- matched_data %>%
    filter(BRCA12_Case == 1) %>%
    pull(person_id) %>%
    sort()

control_ids <- matched_data %>%
    filter(BRCA12_Case == 0) %>%
    pull(person_id) %>%
    sort()

write.csv(case_ids,
          file.path("ch", "data", "matched_case_ids.csv"),
          row.names = FALSE)
write.csv(control_ids,
          file.path("ch", "data", "matched_control_ids.csv"),
          row.names = FALSE)

# ========================
# BALANCE PLOTS
# ========================
# Helper function to extract data for plotting
extract_plot_data <- function(m.out, var, data) {
    treat_var <- as.character(m.out$formula[[2]])
    matched <- match.data(m.out)

    # Handle propensity score
    if (var == "distance") {
        data$distance <- m.out$distance[rownames(data)]
        matched$distance <- m.out$distance[rownames(matched)]
    }

    # Before matching
    before <- data.frame(
        value = data[[var]],
        group = ifelse(data[[treat_var]] == 1, "BRCA1/2 Carriers", "Non-carriers"),
        period = "Before"
    )

    # After matching
    after <- data.frame(
        value = matched[[var]],
        group = ifelse(matched[[treat_var]] == 1, "BRCA1/2 Carriers", "Non-carriers"),
        period = "After"
    )

    rbind(before, after) %>%
        mutate(period = factor(period, levels = c("Before", "After")))
}

# Plot continuous variables (density plots)
plot_continuous <- function(m.out, var, var_label, data, save_dir) {
    df <- extract_plot_data(m.out, var, data)

    p <- ggplot(df, aes(x = value, fill = group, color = group)) +
        geom_density(alpha = 0.3, linewidth = 0.8) +
        facet_wrap(~ period, nrow = 1) +
        scale_fill_manual(values = c("BRCA1/2 Carriers" = "#C2185B", "Non-carriers" = "#1565C0")) +
        scale_color_manual(values = c("BRCA1/2 Carriers" = "#C2185B", "Non-carriers" = "#1565C0")) +
        labs(title = var_label, x = var_label, y = "Density") +
        theme_minimal(base_size = 12) +
        theme(
            legend.title = element_blank(),
            legend.position = "bottom",
            strip.text = element_text(face = "bold"),
            plot.title = element_text(face = "bold", hjust = 0.5)
        )

    filename <- paste0("psm_", gsub("[^A-Za-z0-9]", "_", var), ".pdf")
    ggsave(file.path(save_dir, filename), p, width = 8, height = 5, dpi = 300)
    message("Saved: ", filename)
}

# Plot discrete variables (bar plots)
plot_discrete <- function(m.out, var, var_label, data, save_dir) {
    df <- extract_plot_data(m.out, var, data)

    # Calculate proportions
    df_prop <- df %>%
        group_by(period, group, value) %>%
        summarise(n = n(), .groups = "drop") %>%
        group_by(period, group) %>%
        mutate(prop = n / sum(n))

    p <- ggplot(df_prop, aes(x = factor(value), y = prop, fill = group)) +
        geom_col(position = "dodge", alpha = 0.8) +
        facet_wrap(~ period, nrow = 1) +
        scale_fill_manual(values = c("BRCA1/2 Carriers" = "#C2185B", "Non-carriers" = "#1565C0")) +
        labs(title = var_label, x = var_label, y = "Proportion") +
        theme_minimal(base_size = 12) +
        theme(
            legend.title = element_blank(),
            legend.position = "bottom",
            strip.text = element_text(face = "bold"),
            plot.title = element_text(face = "bold", hjust = 0.5)
        )

    filename <- paste0("psm_", gsub("[^A-Za-z0-9]", "_", var), ".pdf")
    ggsave(file.path(save_dir, filename), p, width = 10, height = 5, dpi = 300)
    message("Saved: ", filename)
}

# Main function to plot all variables
plot_psm_balance <- function(m.out, data,
                             continuous_vars = NULL,
                             discrete_vars = NULL,
                             var_labels = NULL,
                             save_dir = "ch/figures") {

    if (is.null(var_labels)) {
        var_labels <- c(continuous_vars, discrete_vars)
        names(var_labels) <- var_labels
    }

    # Plot continuous variables
    if (!is.null(continuous_vars)) {
        for (var in continuous_vars) {
            plot_continuous(m.out, var, var_labels[var], data, save_dir)
        }
    }

    # Plot discrete variables
    if (!is.null(discrete_vars)) {
        for (var in discrete_vars) {
            plot_discrete(m.out, var, var_labels[var], data, save_dir)
        }
    }
}

continuous_vars <- c("Sample_age", "distance", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6")
discrete_vars <- c("Sequenced_gender", "Batch", "Smoke_History", "Class")

var_labels <- c(
    "Sample_age" = "Age at Sample Collection",
    "distance" = "Propensity Score",
    "PC1" = "Principal Component 1",
    "PC2" = "Principal Component 2",
    "PC3" = "Principal Component 3",
    "PC4" = "Principal Component 4",
    "PC5" = "Principal Component 5",
    "PC6" = "Principal Component 6",
    "Sequenced_gender" = "Sex",
    "Batch" = "Sequencing Batch",
    "Smoke_History" = "Smoking History",
    "Class" = "Ancestry"
)

plot_psm_balance(
    m.out5,
    data = unaff,
    continuous_vars = continuous_vars,
    discrete_vars = discrete_vars,
    var_labels = var_labels,
    save_dir = file.path("ch", "figures")
)

