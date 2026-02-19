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
# case lists
b1_case <- read.csv(file.path("ch", "data", "pmbb_brca1_case.csv"))$x
b2_case <- read.csv(file.path("ch", "data", "pmbb_brca2_case.csv"))$x
b12_case <- unique(c(b1_case, b2_case))

# progeny - BRCA12
progeny <- read_excel(here("ch", "ss", "brca_carriers_ch_freq_w_seen_in_crep_20251020.xlsx"), sheet = "Data_from_master_table")
progeny <- progeny %>% filter(DNA == "D", !is.na(SampNum))
progeny$SampNum <- as.numeric(progeny$SampNum)

up <- read.csv(here("simplexo", "data", "simplexo_up_map.csv"))

# PMBB
COV <- here("PMBB", "3.0", "PMBB-Release-2024-3.0_covariates.txt")
PER <- here("PMBB", "3.0", "PMBB-Release-2024-3.0_phenotype_person.txt")
cov <- fread(COV, header = TRUE)
person <- fread(PER, header = TRUE)

# ========================
# SMOKING HISTORY
# ========================
### PMBB ###
substance <- fread(here("PMBB", "3.0", "PMBB-Release-2024-3.1_phenotype_substance_use.txt"), header = TRUE)

# turn "Not Asked" into NA values
substance <- substance %>%
    mutate(TOBACCO_USER = ifelse(TOBACCO_USER == "Not Asked", NA, TOBACCO_USER))

unique(substance$TOBACCO_USER)
# "Never"   "Quit"    "Yes"     NA        "Passive"

# collapse down values in this column
merge_duplicates2 <- function(df, group_col) {
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

substance_merge <- merge_duplicates2(substance, "person_id")

# non smokers will only have Never and no other values
pmbb_non_smoker <- substance_merge %>% filter(TOBACCO_USER == "Never")
pmbb_non_smoker_ids <- sort(unique(pmbb_non_smoker$person_id))
length(pmbb_non_smoker_ids)
# 28369

# smokers will have anything other than Never and can't be empty
pmbb_smoker <- substance_merge %>% filter(!is.na(TOBACCO_USER) & TOBACCO_USER != "Never")
pmbb_smoker_ids <- sort(unique(pmbb_smoker$person_id))
length(pmbb_smoker_ids)
# 25638

### PROGENY - BRCA12 ###
progeny_merged <- merge_duplicates(progeny, "SampNum")

# filter down to PMBB only
progeny_pmbb <- progeny_merged %>%
    left_join(up, by = "SampNum") %>%
    filter(!is.na(PMBB_ID))
length(unique(progeny_pmbb$PMBB_ID))
# 1278

table(progeny_pmbb$SmokingEver, useNA = "ifany")
# Never      No Unknown     Yes    <NA>
#     775      62      22     367      52

# get rid of unknowns or NA
progeny_pmbb <- progeny_pmbb %>%
    filter(SmokingEver != "Unknown", !is.na(SmokingEver))

# non smokers have Never, No
progeny_non_smoker <- progeny_pmbb %>% filter(SmokingEver %in% c("Never", "No"))
progeny_non_smoker_ids <- sort(unique(progeny_non_smoker$PMBB_ID))
length(progeny_non_smoker_ids)
# 837

# smokers have yes
progeny_smoker <- progeny_pmbb %>% filter(SmokingEver == "Yes")
progeny_smoker_ids <- sort(unique(progeny_smoker$PMBB_ID))
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

# investigate concordance
table(PMBB = cov_smoke$PMBB_SmokeHistory,
      Progeny = cov_smoke$Progeny_SmokeHistory,
      useNA = "ifany")
# Progeny
# PMBB       0     1  <NA>
#     0      507    27 27835
# 1       15   207 25416
# <NA>   315   133  2715

### MERGE SOURCES ###
# if any person has smoking history -> smoker
# if they are non smokers in both -> non smoker
all_smoker_ids <- sort(union(pmbb_smoker_ids, progeny_smoker_ids))
all_non_smoker_ids <- setdiff(union(pmbb_non_smoker_ids, progeny_non_smoker_ids), all_smoker_ids)

length(all_smoker_ids)
# 25798
length(all_non_smoker_ids)
# 28657

# ========================
# CANCER FLAGS
# ========================
cancer_df <- read.csv(file.path("ch", "data", "pmbb_cancer_age.csv"), row.names = NULL)[2:4]
dim(cancer_df)
# 29590

cancer_ids <- read.csv(here("ch", "data", "pmbb_ch_cancer_ids.csv"))$x
length(cancer_ids)
# 29606

# ========================
# PREPARE MATCHING DATASET
# ========================
cov1 <- cov %>%
    mutate(
        BRCA12_Case = ifelse(person_id %in% b12_case, 1, 0),
        BRCA1_Case = ifelse(person_id %in% b1_case, 1, 0),
        BRCA2_Case = ifelse(person_id %in% b2_case, 1, 0),
        BRCA1_BRCA2_Case = ifelse((person_id %in% b1_case & person_id %in% b2_case), 1, 0),
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
        Cancer = ifelse(person_id %in% cancer_ids, 1, 0),
    ) %>% left_join(
        cancer_df[1:2], by = c("person_id" = "PMBB_ID")
    ) %>% mutate(
        Class = ifelse(Class %in% c("UNKNOWN1", "UNKNOWN2", "UNKNOWN3", "UNKNOWN0"), "UNKNOWN", Class)
    )


table(cov1$Smoke_History, useNA = "ifany")
# 0     1  <NA>
# 28657 25798  2715

# ========================
# CREATE CANCER STRATA
# ========================
cov2 <- cov1 %>%
    mutate(
        Strata = case_when(
            # Strata 1: No cancer
            Cancer == 0 ~ 1,
            # Strata 2: Cancer developed >1 year AFTER sample collection
            Cancer == 1 & !is.na(CaDxAge) & CaDxAge > (Sample_age + 1) ~ 2,
            # Strata 3: Cancer diagnosed BEFORE/AT sample collection
            Cancer == 1 & !is.na(CaDxAge) & CaDxAge <= Sample_age ~ 3,
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
# No cancer Post-sample cancer (>1yr)         Pre-sample cancer
# 27564                      3399                     23596
# Uncertain
# 2611

table(cov2$Class, useNA = "ifany")
# AFR     AMR     EAS     EUR     SAS UNKNOWN
# 12113     761     907   42133     743     513

table(cov2$Smoke_History, useNA = "ifany")
# 0     1  <NA>
# 28657 25798  2715

t1 <- table(Carrier = cov2$Carrier, Strata = cov2$Strata_Label, useNA = "ifany")
print(t1)
# Strata
# Carrier       No cancer Post-sample cancer (>1yr) Pre-sample cancer Uncertain
# BRCA1             320                        69               439        33
# BRCA1+BRCA2         2                         0                 7         0
# BRCA2             322                        77               449        45
# Non-carrier     26920                      3253             22701      2533
write.csv(t1, file.path("ch", "data", "strata_by_germline_table.csv"))

t2 <- table(Carrier = cov2$Carrier, Smoking = cov2$Smoke_History, useNA = "ifany")
print(t2)
# Smoking
# Carrier        0     1  <NA>
#     BRCA1      509   274    78
# BRCA1+BRCA2     6     3     0
# BRCA2         557   278    58
# Non-carrier 27585 25243  2579

write.csv(t2, file.path("ch", "data", "smoking_by_germline_table.csv"))

t3 <- table(Carrier = cov2$Carrier, Smoking = cov2$Batch, useNA = "ifany")
print(t3)
# Strata
# Carrier           1     2
# BRCA1         150   711
# BRCA1+BRCA2     1     8
# BRCA2         259   634
# Non-carrier 43503 11904

dim(cov2)
cov3 <- cov2 %>% filter(!is.na(Smoke_History), !is.na(Sample_age))
dim(cov3)
# 54455

# ========================
# UNAFFECTED ANALYSIS
# ========================
unaff <- cov3 %>% filter(Strata %in% c(1, 2))

dim(unaff)
# 28822    21

table(unaff$BRCA12_Case)
# 0     1
# 28125   697

table(unaff$BRCA1_Case)
# 0     1
# 28482   340

table(unaff$BRCA2_Case)
# 0     1
# 28463   359

table(unaff$BRCA1_BRCA2_Case)
# 0     1
# 28820     2

table(unaff$Batch)
# 1     2
# 23384  5438

table(unaff$Sequenced_gender)
# Female   Male
# 14684  14138

table(unaff$Smoke_History)
# 0     1
# 16000 12822

range(unaff$Sample_age)
# 18 - 90

# ========================
# CHECK INITIAL IMBALANCE
# ========================
unaff$Sequenced_gender <- as.factor(unaff$Sequenced_gender)
unaff$Smoke_History <- as.factor(unaff$Smoke_History)
unaff$Batch <- as.factor(unaff$Batch)
unaff$PC1 <- as.numeric(unaff$PC1)
unaff$PC2 <- as.numeric(unaff$PC2)
unaff$PC3 <- as.numeric(unaff$PC3)
unaff$PC4 <- as.numeric(unaff$PC4)
unaff$PC5 <- as.numeric(unaff$PC5)
unaff$PC6 <- as.numeric(unaff$PC6)

match_formula <- as.formula(paste(
    "BRCA12_Case ~ Sample_age + Sequenced_gender + Smoke_History + Batch +",
    paste(paste0("PC", 1:6), collapse = " + ")
))

# prematch object
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
# distance                      0.1316        0.0215          1.1750     4.3665    0.3606   0.6333
# Sample_age                   42.9803       50.7815         -0.5848     0.6062    0.1105   0.2600
# Sequenced_genderFemale        0.7676        0.5031          0.6262          .    0.2645   0.2645
# Sequenced_genderMale          0.2324        0.4969         -0.6262          .    0.2645   0.2645
# Smoke_History0                0.7102        0.5513          0.3502          .    0.1589   0.1589
# Smoke_History1                0.2898        0.4487         -0.3502          .    0.1589   0.1589
# Batch1                        0.2181        0.8260         -1.4722          .    0.6079   0.6079
# Batch2                        0.7819        0.1740          1.4722          .    0.6079   0.6079
# PC1                           0.0106        0.0016          0.8262     0.2280    0.0651   0.1955
# PC2                           0.0443        0.0334          0.8327     0.2758    0.0610   0.2236
# PC3                          -0.0174       -0.0148         -0.2571     0.4991    0.0487   0.1542
# PC4                          -0.0096       -0.0035         -0.4206     1.0340    0.1631   0.2375
# PC5                           0.0059        0.0113         -0.4913     0.9661    0.1474   0.2160
# PC6                          -0.0005        0.0002         -0.0930     0.7457    0.0177   0.0375

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
summary(m.out2, un = FALSE)
# Summary of Balance for Matched Data:
#     Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean eCDF Max Std. Pair Dist.
# distance                      0.1277        0.1272          0.0054     1.0135    0.0005   0.0162          0.0116
# Sample_age                   42.9196       42.6119          0.0231     0.9825    0.0059   0.0310          0.1834
# Sequenced_genderFemale        0.7622        0.7622         -0.0000          .    0.0000   0.0000          0.0000
# Sequenced_genderMale          0.2378        0.2378         -0.0000          .    0.0000   0.0000          0.0000
# Smoke_History0                0.7120        0.7120         -0.0000          .    0.0000   0.0000          0.0000
# Smoke_History1                0.2880        0.2880         -0.0000          .    0.0000   0.0000          0.0000
# Batch1                        0.2245        0.2245         -0.0000          .    0.0000   0.0000          0.0000
# Batch2                        0.7755        0.7755         -0.0000          .    0.0000   0.0000          0.0000
# PC1                           0.0106        0.0105          0.0079     0.8576    0.0166   0.0443          0.2627
# PC2                           0.0443        0.0439          0.0245     0.7961    0.0333   0.0746          0.3288
# PC3                          -0.0175       -0.0168         -0.0672     0.5401    0.0184   0.0465          0.7989
# PC4                          -0.0094       -0.0094          0.0018     1.7935    0.0164   0.0406          0.4780
# PC5                           0.0063        0.0064         -0.0116     1.4383    0.0144   0.0340          0.5009
# PC6                          -0.0004       -0.0005          0.0101     0.8797    0.0080   0.0266          1.0881
#
# Sample Sizes:
#     Control Treated
# All           28125.       697
# Matched (ESS)  1252.27     677
# Matched        1299.       677
# Unmatched     26826.        20
# Discarded         0.         0

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
summary(m.out3, un = FALSE)
# Summary of Balance for Matched Data:
#     Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean eCDF Max Std. Pair Dist.
# distance                      0.1277        0.1271          0.0067     1.0161    0.0006   0.0155          0.0139
# Sample_age                   42.9196       42.6037          0.0237     0.9640    0.0067   0.0357          0.1861
# Sequenced_genderFemale        0.7622        0.7622          0.0000          .    0.0000   0.0000          0.0000
# Sequenced_genderMale          0.2378        0.2378         -0.0000          .    0.0000   0.0000          0.0000
# Smoke_History0                0.7120        0.7120         -0.0000          .    0.0000   0.0000          0.0000
# Smoke_History1                0.2880        0.2880         -0.0000          .    0.0000   0.0000          0.0000
# Batch1                        0.2245        0.2245         -0.0000          .    0.0000   0.0000          0.0000
# Batch2                        0.7755        0.7755         -0.0000          .    0.0000   0.0000          0.0000
# PC1                           0.0106        0.0108         -0.0130     0.9031    0.0193   0.0588          0.2710
# PC2                           0.0443        0.0441          0.0091     0.8237    0.0348   0.0886          0.3279
# PC3                          -0.0175       -0.0171         -0.0408     0.5843    0.0213   0.0515          0.7971
# PC4                          -0.0094       -0.0092         -0.0118     1.4590    0.0153   0.0384          0.4864
# PC5                           0.0063        0.0063         -0.0018     1.3522    0.0143   0.0318          0.5028
# PC6                          -0.0004       -0.0004         -0.0090     0.8722    0.0055   0.0226          1.0932
#
# Sample Sizes:
#     Control Treated
# All           28125.       697
# Matched (ESS)  1675.79     677
# Matched        1854.       677
# Unmatched     26271.        20
# Discarded         0.         0

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
summary(m.out4, un = FALSE)
# Summary of Balance for Matched Data:
#     Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean eCDF Max Std. Pair Dist.
# distance                      0.1277        0.1269          0.0082     1.0169    0.0007   0.0149          0.0164
# Sample_age                   42.9196       42.5391          0.0285     0.9560    0.0075   0.0404          0.1872
# Sequenced_genderFemale        0.7622        0.7622          0.0000          .    0.0000   0.0000          0.0000
# Sequenced_genderMale          0.2378        0.2378          0.0000          .    0.0000   0.0000          0.0000
# Smoke_History0                0.7120        0.7120          0.0000          .    0.0000   0.0000          0.0000
# Smoke_History1                0.2880        0.2880          0.0000          .    0.0000   0.0000          0.0000
# Batch1                        0.2245        0.2245          0.0000          .    0.0000   0.0000          0.0000
# Batch2                        0.7755        0.7755          0.0000          .    0.0000   0.0000          0.0000
# PC1                           0.0106        0.0109         -0.0214     0.9197    0.0218   0.0627          0.2903
# PC2                           0.0443        0.0441          0.0102     0.8095    0.0360   0.0907          0.3393
# PC3                          -0.0175       -0.0173         -0.0216     0.6315    0.0223   0.0545          0.7795
# PC4                          -0.0094       -0.0091         -0.0169     1.4157    0.0164   0.0430          0.5048
# PC5                           0.0063        0.0063          0.0003     1.3500    0.0145   0.0345          0.5187
# PC6                          -0.0004       -0.0003         -0.0105     0.8438    0.0070   0.0267          1.1066
#
# Sample Sizes:
#     Control Treated
# All           28125.       697
# Matched (ESS)  1957.97     677
# Matched        2327.       677
# Unmatched     25798.        20
# Discarded         0.         0

# Distribution of matches per carrier
match_counts <- matched_data %>%
    group_by(subclass) %>%
    summarise(
        n_total = n(),
        n_controls = sum(BRCA12_Case == 0),
        n_treated = sum(BRCA12_Case == 1)
    ) %>%
    pull(n_controls)

table(match_counts)
# 1   2   3   4
# 55  67  82 473

# 1:4 PROBIT
m.out42 <- matchit(
    match_formula,
    data = unaff,
    method = "nearest",
    distance = "glm",
    ratio = 4,
    caliper = c(0.2, # use 0.2 STD of logit(PS)
                Sample_age = 5), # +- 5 years for sample age
    std.caliper = c(TRUE, FALSE), # goes with caliper field
    exact = ~ Sequenced_gender + Batch + Smoke_History,
    replace = FALSE,
    link = "probit"
)
summary(m.out42, un = FALSE)
# Summary of Balance for Matched Data:
#     Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean eCDF Max Std. Pair Dist.
# distance                      0.1179        0.1171          0.0102     1.0185    0.0007   0.0158          0.0176
# Sample_age                   42.9561       42.5563          0.0300     0.9536    0.0080   0.0415          0.1916
# Sequenced_genderFemale        0.7635        0.7635         -0.0000          .    0.0000   0.0000          0.0000
# Sequenced_genderMale          0.2365        0.2365         -0.0000          .    0.0000   0.0000          0.0000
# Smoke_History0                0.7124        0.7124         -0.0000          .    0.0000   0.0000          0.0000
# Smoke_History1                0.2876        0.2876         -0.0000          .    0.0000   0.0000          0.0000
# Batch1                        0.2219        0.2219         -0.0000          .    0.0000   0.0000          0.0000
# Batch2                        0.7781        0.7781         -0.0000          .    0.0000   0.0000          0.0000
# PC1                           0.0106        0.0101          0.0534     0.7350    0.0151   0.0468          0.3830
# PC2                           0.0443        0.0437          0.0454     0.7305    0.0380   0.0842          0.4051
# PC3                          -0.0175       -0.0173         -0.0193     0.7119    0.0177   0.0518          0.7602
# PC4                          -0.0095       -0.0090         -0.0349     1.7819    0.0267   0.0562          0.5341
# PC5                           0.0061        0.0062         -0.0093     1.3686    0.0091   0.0270          0.4935
# PC6                          -0.0004       -0.0002         -0.0305     0.8443    0.0101   0.0293          1.1027
#
# Sample Sizes:
#     Control Treated
# All           28125.      697
# Matched (ESS)  1974.3     685
# Matched        2353.      685
# Unmatched     25772.       12
# Discarded         0.        0

x <- summary(m.out42)$sum.matched
write.csv(x, file.path("ch", "data", "ch_psm_matched4_balance_table2.csv"), row.names = FALSE)

### 1:5
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
summary(m.out5, un = FALSE)
# Call:
#     matchit(formula = match_formula, data = unaff, method = "nearest",
#             distance = "glm", link = "probit", exact = ~Sequenced_gender +
#                 Batch + Smoke_History, replace = FALSE, caliper = c(0.2,
#                                                                     Sample_age = 5), std.caliper = c(TRUE, FALSE), ratio = 5)
#
# Summary of Balance for Matched Data:
#     Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean eCDF Max Std. Pair Dist.
# distance                      0.1179        0.1171          0.0111     1.0186    0.0008   0.0172          0.0187
# Sample_age                   42.9561       42.5415          0.0311     0.9501    0.0081   0.0415          0.1915
# Sequenced_genderFemale        0.7635        0.7635          0.0000          .    0.0000   0.0000          0.0000
# Sequenced_genderMale          0.2365        0.2365          0.0000          .    0.0000   0.0000          0.0000
# Smoke_History0                0.7124        0.7124          0.0000          .    0.0000   0.0000          0.0000
# Smoke_History1                0.2876        0.2876          0.0000          .    0.0000   0.0000          0.0000
# Batch1                        0.2219        0.2219          0.0000          .    0.0000   0.0000          0.0000
# Batch2                        0.7781        0.7781          0.0000          .    0.0000   0.0000          0.0000
# PC1                           0.0106        0.0101          0.0498     0.7408    0.0151   0.0473          0.3995
# PC2                           0.0443        0.0437          0.0453     0.7311    0.0372   0.0829          0.4147
# PC3                          -0.0175       -0.0173         -0.0198     0.6877    0.0193   0.0533          0.7693
# PC4                          -0.0095       -0.0091         -0.0291     1.8541    0.0260   0.0570          0.5340
# PC5                           0.0061        0.0063         -0.0147     1.3956    0.0091   0.0252          0.4903
# PC6                          -0.0004       -0.0002         -0.0329     0.8436    0.0099   0.0299          1.1096
#
# Sample Sizes:
#     Control Treated
# All           28125.       697
# Matched (ESS)  2145.85     685
# Matched        2733.       685
# Unmatched     25392.        12
# Discarded         0.         0

# ========================
# EXTRACT MATCHED DATA
# ========================
matched_data <- match_data(m.out4)
head(matched_data)
dim(matched_data)
# 3004

write.csv(matched_data, file.path("ch", "data", "ch_psm_matched4_data.csv"), row.names = FALSE)

# split into case/controls
case_ids <- matched_data %>%
    filter(BRCA12_Case == 1) %>%
    pull(person_id) %>%
    sort()

control_ids <- matched_data %>%
    filter(BRCA12_Case == 0) %>%
    pull(person_id) %>%
    sort()

length(case_ids)
length(control_ids)
intersect(control_ids, case_ids)
all_ids <- sort(c(control_ids, case_ids))

write.csv(case_ids, file.path("ch", "data", "ch_psm_matched4_case_ids.csv"), row.names = FALSE)
write.csv(control_ids, file.path("ch", "data", "ch_psm_matched4_control_ids.csv"), row.names = FALSE)
write.csv(all_ids, file.path("ch", "data", "ch_psm_matched4_case_control_ids.csv"), row.names = FALSE)

# ========================
# CHECK AGE
# ========================
md <- match.data(m.out4)

# for each treated unit, max |age difference| among its matched controls
age_diff_check <- md %>%
    group_by(subclass) %>%
    summarize(age_range = max(Sample_age) - min(Sample_age),
              n = n(),
              .groups="drop") %>%
    summarize(max_age_range = max(age_range), p95 = quantile(age_range, .95))

age_diff_check
# max_age_range   p95
# <dbl> <dbl>
#     1           9.9  8.82

age_spread <- md %>%
    group_by(subclass) %>%
    summarise(
        age_range = max(Sample_age) - min(Sample_age),
        .groups = "drop"
    )

summary(age_spread$age_range)
# psm
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.100   3.900   5.000   5.298   6.900   9.900

quantile(age_spread$age_range, probs = c(0, .25, .5, .75, .9, .95, 1))
# psm
# 0%  25%  50%  75%  90%  95% 100%
# 0.10 3.90 5.00 6.90 8.30 8.82 9.90

# ========================
# PLOTS
# ========================
### love plot - SMD before and after matching
love_plot <- love.plot(
    m.out4,
    stats = "mean.diffs",
    threshold = c(m = 0.1),  # 0.1 threshold line
    var.order = "unadjusted",
    abs = TRUE,
    stars = "std", # those with mean diff that have been standardized
    colors = c("#B2182B", "#2166AC")
)
ggsave(file.path("ch", "figures", "matching_love_plot.pdf"),
       love_plot, width = 6, height = 4)

### Balance table
bal_tab <- bal.tab(
    m.out4,
    m.threshold = 0.1,  # Flag variables with SMD > 0.1
    un = TRUE,           # Show unmatched balance
    disp.v.ratio = TRUE  # Show variance ratios
)
print(bal_tab)
# Balance Measures
# Type Diff.Un V.Ratio.Un Diff.Adj    M.Threshold V.Ratio.Adj
# distance              Distance  1.1750     4.3665   0.0082 Balanced, <0.1      1.0169
# Sample_age             Contin. -0.5848     0.6062   0.0285 Balanced, <0.1      0.9560
# Sequenced_gender_Male   Binary -0.2645          .   0.0000 Balanced, <0.1           .
# Smoke_History           Binary -0.1589          .   0.0000 Balanced, <0.1           .
# Batch_2                 Binary  0.6079          .   0.0000 Balanced, <0.1           .
# PC1                    Contin.  0.8262     0.2280  -0.0214 Balanced, <0.1      0.9197
# PC2                    Contin.  0.8327     0.2758   0.0102 Balanced, <0.1      0.8095
# PC3                    Contin. -0.2571     0.4991  -0.0216 Balanced, <0.1      0.6315
# PC4                    Contin. -0.4206     1.0340  -0.0169 Balanced, <0.1      1.4157
# PC5                    Contin. -0.4913     0.9661   0.0003 Balanced, <0.1      1.3500
# PC6                    Contin. -0.0930     0.7457  -0.0105 Balanced, <0.1      0.8438
#
# Balance tally for mean differences
# count
# Balanced, <0.1        11
# Not Balanced, >0.1     0
#
# Variable with the greatest mean difference
# Variable Diff.Adj    M.Threshold
# Sample_age   0.0285 Balanced, <0.1
#
# Sample sizes
# Control Treated
# All                  28125.       697
# Matched (ESS)         1957.97     677
# Matched (Unweighted)  2327.       677
# Unmatched            25798.        20

write.csv(bal_tab$Balance, file.path("ch", "data", "ch_psm_matched4_balance_table.csv"), row.names = FALSE)

# ========================
# INDIVIDUAL BALANCE PLOTS
# ========================
extract_plot_data <- function(m.out, var, data) {
    treat_var <- as.character(m.out$formula[[2]])
    matched <- match.data(m.out)

    if (var == "distance") {
        data$distance <- m.out$distance[rownames(data)]
        matched$distance <- m.out$distance[rownames(matched)]
    }

    before <- data.frame(
        value = data[[var]],
        group = ifelse(data[[treat_var]] == 1, "BRCA1/2 Carriers", "Non-carriers"),
        period = "Before"
    )

    after <- data.frame(
        value = matched[[var]],
        group = ifelse(matched[[treat_var]] == 1, "BRCA1/2 Carriers", "Non-carriers"),
        period = "After"
    )

    rbind(before, after) %>%
        mutate(period = factor(period, levels = c("Before", "After")))
}

plot_continuous <- function(m.out, var, var_label, data, save_dir, name) {
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

    filename <- paste0(name, "_", gsub("[^A-Za-z0-9]", "_", var), ".pdf")
    ggsave(file.path(save_dir, filename), p, width = 8, height = 5, dpi = 300)
    message("Saved: ", filename)

    return(p)
}

plot_discrete <- function(m.out, var, var_label, data, save_dir, name) {
    df <- extract_plot_data(m.out, var, data)

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

    filename <- paste0(name, "_", gsub("[^A-Za-z0-9]", "_", var), ".pdf")
    ggsave(file.path(save_dir, filename), p, width = 10, height = 5, dpi = 300)
    message("Saved: ", filename)

    return(p)
}

# ========================
# COMPOSITE PLOT FUNCTION
# ========================
create_composite_balance_plot <- function(m.out, data,
                                          continuous_vars = NULL,
                                          discrete_vars = NULL,
                                          var_labels = NULL,
                                          save_dir = "ch/figures",
                                          ncol = 3,
                                          name = "psm") {

    library(patchwork)

    if (is.null(var_labels)) {
        var_labels <- c(continuous_vars, discrete_vars)
        names(var_labels) <- var_labels
    }

    plot_list <- list()
    all_vars <- names(var_labels)

    for (var in all_vars) {
        if (var %in% continuous_vars) {
            df <- extract_plot_data(m.out, var, data)

            p <- ggplot(df, aes(x = value, fill = group, color = group)) +
                geom_density(alpha = 0.3, linewidth = 0.6) +
                facet_wrap(~ period, nrow = 1) +
                scale_fill_manual(values = c("BRCA1/2 Carriers" = "#C2185B", "Non-carriers" = "#1565C0")) +
                scale_color_manual(values = c("BRCA1/2 Carriers" = "#C2185B", "Non-carriers" = "#1565C0")) +
                labs(title = var_labels[var], x = "", y = "Density") +
                theme_minimal(base_size = 10) +
                theme(
                    legend.position = "none",
                    strip.text = element_text(face = "bold", size = 8),
                    plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
                    axis.title = element_text(size = 8),
                    axis.text = element_text(size = 7)
                )

            plot_list[[var]] <- p

        } else if (var %in% discrete_vars) {
            df <- extract_plot_data(m.out, var, data)

            df_prop <- df %>%
                group_by(period, group, value) %>%
                summarise(n = n(), .groups = "drop") %>%
                group_by(period, group) %>%
                mutate(prop = n / sum(n))

            p <- ggplot(df_prop, aes(x = factor(value), y = prop, fill = group)) +
                geom_col(position = "dodge", alpha = 0.8) +
                facet_wrap(~ period, nrow = 1) +
                scale_fill_manual(values = c("BRCA1/2 Carriers" = "#C2185B", "Non-carriers" = "#1565C0")) +
                labs(title = var_labels[var], x = "", y = "Proportion") +
                theme_minimal(base_size = 10) +
                theme(
                    legend.position = "none",
                    strip.text = element_text(face = "bold", size = 8),
                    plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
                    axis.title = element_text(size = 8),
                    axis.text = element_text(size = 7),
                    axis.text.x = element_text(angle = 45, hjust = 1)
                )

            plot_list[[var]] <- p
        }
    }

    if (length(plot_list) > 0) {
        legend_plot <- ggplot(data.frame(x = 1, y = 1, group = c("BRCA1/2 Carriers", "Non-carriers")),
                              aes(x = x, y = y, fill = group)) +
            geom_point(size = 4, alpha = 0.8) +
            scale_fill_manual(values = c("BRCA1/2 Carriers" = "#C2185B", "Non-carriers" = "#1565C0")) +
            theme_void() +
            theme(
                legend.position = "bottom",
                legend.title = element_blank(),
                legend.text = element_text(size = 12)
            ) +
            guides(fill = guide_legend(override.aes = list(shape = 22, size = 5)))

        legend <- ggplotGrob(legend_plot)$grobs[[which(sapply(ggplotGrob(legend_plot)$grobs, function(x) x$name) == "guide-box")]]

        composite <- wrap_plots(plot_list, ncol = ncol) +
            plot_annotation(
                title = "Covariate Balance: Before vs After Propensity Score Matching",
                theme = theme(
                    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
                    plot.margin = margin(20, 20, 60, 20)
                )
            )

        final_plot <- composite / wrap_elements(legend)
        final_plot <- final_plot + plot_layout(heights = c(1, 0.05))

        filename <- paste0(name, "_composite_balance_plot.pdf")
        ggsave(file.path(save_dir, filename),
               final_plot,
               width = 4 * ncol,
               height = 3 * ceiling(length(plot_list) / ncol) + 1,
               dpi = 300)

        message("Saved composite plot: ", filename)

        return(final_plot)
    }
}

# ========================
# MAIN WRAPPER FUNCTION
# ========================
plot_psm_balance <- function(m.out, data,
                             continuous_vars = NULL,
                             discrete_vars = NULL,
                             var_labels = NULL,
                             save_dir = "ch/figures",
                             create_composite = TRUE,
                             composite_ncol = 3,
                             name = "psm") {

    if (is.null(var_labels)) {
        var_labels <- c(continuous_vars, discrete_vars)
        names(var_labels) <- var_labels
    }

    if (!is.null(continuous_vars)) {
        for (var in continuous_vars) {
            plot_continuous(m.out, var, var_labels[var], data, save_dir, name)
        }
    }

    if (!is.null(discrete_vars)) {
        for (var in discrete_vars) {
            plot_discrete(m.out, var, var_labels[var], data, save_dir, name)
        }
    }

    if (create_composite) {
        composite_plot <- create_composite_balance_plot(
            m.out, data, continuous_vars, discrete_vars,
            var_labels, save_dir, composite_ncol, name
        )
        return(composite_plot)
    }
}

# ========================
# PLOT
# ========================
continuous_vars <- c("Sample_age", "distance", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6")
discrete_vars <- c("Sequenced_gender", "Batch", "Smoke_History", "Class")

var_labels <- c(
    "distance" = "Propensity Score",
    "Sample_age" = "Age at Sample Collection",
    "Class" = "Ancestry",
    "Sequenced_gender" = "Sex",
    "Batch" = "Sequencing Batch",
    "Smoke_History" = "Smoking History",
    "PC1" = "PC1",
    "PC2" = "PC2",
    "PC3" = "PC3",
    "PC4" = "PC4",
    "PC5" = "PC5",
    "PC6" = "PC6"
)

composite_plot <- plot_psm_balance(
    m.out4,
    data = unaff,
    continuous_vars = continuous_vars,
    discrete_vars = discrete_vars,
    var_labels = var_labels,
    save_dir = file.path("ch", "figures"),
    create_composite = TRUE,
    composite_ncol = 3,
    name = "psm4"
)

# ========================
# IDENTIFY 2 PANEL OF NORMALS
# ========================
aml <- read_excel(here("ch", "ss", "20251110_aml_ngs.xlsx"))
cov3 %>% filter(person_id %in% aml$PMBB_ID) %>% dplyr::select(Batch) %>% unique()
# pilot study has all in freeze 2, so can't really tell...

### FREEZE 2 ###
# get the youngest 50 people in freeze 2
youngest_f2 <- cov3 %>% filter(Strata == 1, Smoke_History == 0, Batch == 1) %>%
    arrange(Sample_age) %>%
    slice_head(n = 50)
youngest_f2_ids <- sort(youngest_f2$person_id)

# check ancestry distribution
prop.table(table(youngest_f2$Class))
# AFR     AMR     EAS     EUR     SAS UNKNOWN
# 0.38    0.02    0.02    0.52    0.02    0.04

prop.table(table((unaff %>% filter(Batch == 1))$Class))
# AFR        AMR        EAS        EUR        SAS    UNKNOWN
# 0.27373418 0.01492473 0.01821758 0.66412932 0.01761888 0.01137530

# aml_controls <- c("PMBB3278234229783", "PMBB9930975965330")
# f2_pon <- sort(c(youngest_f2_ids[1:78], aml_controls))

write.csv(youngest_f2_ids, file.path("ch", "data", "pon_f2_pon_ids.csv"), row.names = FALSE)
# write.csv(f2_pon, file.path("ch", "data", "pon_f2_aml_pon_ids.csv"), row.names = FALSE)

### FREEZE 3 ###
youngest_f3 <- cov3 %>% filter(Strata == 1, Smoke_History == 0, Batch == 2) %>%
    arrange(Sample_age) %>%
    slice_head(n = 50)
youngest_f3_ids <- sort(youngest_f3$person_id)

# check ancestry distribution
prop.table(table(youngest_f3$Class))
# AFR     AMR     EAS     EUR     SAS UNKNOWN
# 0.06    0.04    0.02    0.82    0.04    0.02

prop.table(table((unaff %>% filter(Batch == 2))$Class))
# AFR        AMR        EAS        EUR        SAS    UNKNOWN
# 0.06914307 0.01930857 0.02427363 0.85582935 0.01838911 0.01305627

write.csv(youngest_f3_ids, file.path("ch", "data", "pon_f3_pon_ids.csv"), row.names = FALSE)

### save end df
pon_ids <- c(youngest_f2_ids, youngest_f3_ids)

pon_df <- data.frame(
    PMBB_ID = pon_ids,
    Batch   = ifelse(pon_ids %in% youngest_f2_ids, "2", "3"),
    stringsAsFactors = TRUE
)
write.csv(pon_df, file.path("ch", "data", "pon_f23_pon_df.csv"), row.names = FALSE)

# ========================
# IPW
# ========================
matched_data <- match_data(m.out4)
head(matched_data)
dim(matched_data)
# 3004

write.csv(matched_data, file.path("ch", "data", "psm_matched_data.csv"), row.names = FALSE)

# split into case/controls
case_ids <- matched_data %>%
    filter(BRCA12_Case == 1) %>%
    pull(person_id) %>%
    sort()

control_ids <- matched_data %>%
    filter(BRCA12_Case == 0) %>%
    pull(person_id) %>%
    sort()

write.csv(case_ids, file.path("ch", "data", "ch_matched_case_ids.csv"), row.names = FALSE)
write.csv(control_ids, file.path("ch", "data", "ch_matched_control_ids.csv"), row.names = FALSE)

