# ============================================================
# CNA ANALYSIS - BATCH PROCESSING
# ============================================================
source("config.R")
library(survival)
library(survminer)

# ========================
# LOAD DATA & CONVERT
# ========================
master <- read_xlsx(here("tgct", "ss", "20260812_tgct_master.xlsx"))
ss <- read.csv(here("tgct", "ss", "prstcstageirelapse-stageipatients_data_2026-04-14_1446.csv"))

# ========================
# PREP DATA
# ========================
sum(is.na(ss$testicular_tumor_histology))
# 0
sum(is.na(ss$relapse_at_stage_1))
# 8
sum(is.na(ss$adjuvant_chemotherapy_s1))
# 342

ss$relapse_at_stage_1 <- as.factor(ss$relapse_at_stage_1)
ss$adjuvant_chemotherapy_s1 <- as.factor(ss$adjuvant_chemotherapy_s1)

ss <- ss |> mutate(
    date1_of_relapse = as.Date(date1_of_relapse, format = "%m/%d/%y"),
    date2_of_relapse = as.Date(date2_of_relapse, format = "%m/%d/%y"),
    date_of_diagnosis = as.Date(date_of_diagnosis, format = "%m/%d/%y"),
    date_last_contact = as.Date(date_last_contact, format = "%m/%d/%y")
) |>
    filter(!is.na(testicular_tumor_histology), !is.na(relapse_at_stage_1))

# time-to-event: time to relapse for events, time to last contact for censored
# (previously used date1_of_relapse for everyone, which is NA for the ~393
# patients who never relapsed -- that silently dropped all censored patients
# from the survival fit instead of treating them as censored)
ss$time_to_event <- as.numeric(difftime(
    if_else(as.character(ss$relapse_at_stage_1) == "1", ss$date1_of_relapse, ss$date_last_contact),
    ss$date_of_diagnosis,
    units = "days"
)) / 365.25
sum(is.na(ss$time_to_event))

# ========================
# KAPLAN SURVIVAL CURVES
# ========================
# join = deferred, TODO later: left join the PRS from master (1116 rows)

table(ss$relapse_at_stage_1)
# 0   1
# 393  69

# Create survival object
surv_obj <- Surv(
    time = ss$time_to_event,
    event = as.numeric(as.character(ss$relapse_at_stage_1))
)

# Overall relapse-free survival
fit <- survfit(surv_obj ~ 1, data = ss)

ggsurvplot(
    fit,
    data = ss,
    risk.table = TRUE,
    conf.int = TRUE,
    xlab = "Time",
    ylab = "Relapse-free survival probability"
)


# ========================
# STRATIFIED BY
# ========================





