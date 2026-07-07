# ========================
# CONFIGURATION
# ========================
N_WORKERS <- 20

# BASE_DIR <- "/home/lijz/thyroid"
BASE_DIR <- "~/Documents/HPC_share/thyroid"
R_DIR <- file.path(BASE_DIR, "R")
SS_DIR <- file.path(BASE_DIR, "ss")
DATA_DIR <- file.path(BASE_DIR, "data")
IDAT_DIR <- file.path(BASE_DIR, "IDAT")

source(file.path(R_DIR, "functions.R"))
source(file.path(R_DIR, "color_keys.R"))
source(file.path(R_DIR, "load_packages.R"))

# ============================================================
# CLEAN THCA DATA
# ============================================================
key <- read.table(file.path(DATA_DIR, "key.txt"), header = FALSE)
master <- read_xlsx(file.path(SS_DIR, "adult_master.xlsx"))
# sum(!(key$V2 %in% master$Patient_ID))
key2 <- key[grepl("(01A$)", key$V3),]
master2 <- master[grepl("(01$)", master$Sample_ID),]

# ========================
# LOAD DATA
# ========================
ped_pvals <- as.data.frame(readRDS(file.path(DATA_DIR, "ped_pvals.rds")))
ped_sdf <- as.data.frame(readRDS(file.path(DATA_DIR, "ped_sigdf_mask.rds")))
ped_betas <- as.data.frame(readRDS(file.path(DATA_DIR, "ped_betas_mask.rds")))

ss <- as.data.frame(read_excel(file.path(SS_DIR, "20231102_thyroid_master.xlsx")))
ss_new <- read_excel(file.path(SS_DIR, "pediatric_master.xlsx")) %>%
    dplyr::filter(is.na(Lymph_Node))

# ========================
#  SAMPLE QC STATS - SHOULD THIS BE MASKED OR NOT MASKED PROBES ??
# ========================
### PVALS ###
mean_pvals <- apply(ped_pvals, 2, mean, na.rm = TRUE)
write.csv(mean_pvals, file.path(DATA_DIR, "ped_pvals_mean.csv"))

# higher pvals...may want to remove?
mean_pvals[!(mean_pvals < 0.05)]
# 207219750052_R01C01 207686140043_R06C01 209392270039_R04C01 209392270039_R05C01 209392270039_R06C01
# 0.13334055          0.06265042          0.05187620          0.08870988          0.11362932
# 209392270039_R08C01 209443890176_R03C01 209443890176_R06C01 209443890176_R07C01
# 0.06982606          0.05454046          0.07019812          0.05532481

### PROBE SUCCESS RATE ###
rate <- probeSuccessRate(sdf, mask=FALSE)

### MISSING PROBES ###
missing <- apply(ped_betas, 2, sum(is.na))
na_counts <- apply(ped_betas, 2, function(x) sum(is.na(x)))
na_fractions <- na_counts / nrow(ped_betas)

# ========================
# APPLY CURRENT CLASSIFIER
# ========================
ped_betas <- as.data.frame(readRDS(file.path(DATA_DIR, "ped_betas_imp.rds")))
MODEL_DIR <- file.path(DATA_DIR, "final_driver_model")
imp <- readRDS(file.path(MODEL_DIR, "fold_importance.rds"))

sel_probes3k <- rownames(imp %>% head(3000))
ped_betas3k  <- ped_betas[sel_probes3k, ]

pred <- predict(model, t(ped_betas3k))
prob <- predict(model, t(ped_betas3k), type = "prob")

pred_df <- data.frame(
    Sample_ID = colnames(ped_betas3k),
    Predicted_Class = pred
)
prob_df <- as.data.frame(prob)
prob_df$Sample_ID <- rownames(prob_df)

results_df <- merge(pred_df, prob_df, by = "Sample_ID")

write.csv(results_df, file.path(DATA_DIR, "new_batch_predictions.csv"), row.names = FALSE)

# ========================
# CLUSTERING ANALYSIS
# ========================
ped_betas <- as.data.frame(readRDS(file.path(DATA_DIR, "ped_betas_qc.rds")))


# ========================
# CLUSTERING ANALYSIS
# ========================
ped_betas <- as.data.frame(readRDS(file.path(DATA_DIR, "ped_betas_qc.rds")))




