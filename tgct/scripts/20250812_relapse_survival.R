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
    filter(!is.na(testicular_tumor_histology), !is.na(relapse_at_stage_1), !is.na(prs_2026))

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
# 336  56

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
# STRATIFIED BY PRS
# ========================
ss$prs_2026 <- as.numeric(ss$prs_2026 )
# create tertiles
prs_tertile_breaks <- quantile(ss$prs_2026, probs = c(1 / 3, 2 / 3), na.rm = TRUE)
ss$prs_tertile <- cut(
    ss$prs_2026,
    breaks = c(-Inf, prs_tertile_breaks, Inf),
    labels = c("Low", "Mid", "High")
)
table(ss$prs_tertile)
# Low  Mid High
# 131  130  131

fit_tertile <- survfit(surv_obj ~ prs_tertile, data = ss)

p_tertile <- ggsurvplot(
    fit_tertile,
    data = ss,
    risk.table = TRUE,
    conf.int = TRUE,
    pval = TRUE,
    xlab = "Time",
    ylab = "Relapse-free survival probability",
    legend.title = "PRS tertile",
    legend.labs = c("Low", "Mid", "High")
)

png(here("tgct", "figures", "tgct_relapse_survival_prs_tertile.png"), width = 7, height = 6.5, units = "in", res = 300)
print(p_tertile)
dev.off()

# also compare tertiles as a two-level variable: top tertile vs bottom two combined
ss$prs_group2 <- if_else(ss$prs_tertile == "High", "High", "Low/Mid")
table(ss$prs_group2)

fit_group2 <- survfit(surv_obj ~ prs_group2, data = ss)

p_group2 <- ggsurvplot(
    fit_group2,
    data = ss,
    risk.table = TRUE,
    conf.int = TRUE,
    pval = TRUE,
    xlab = "Time",
    ylab = "Relapse-free survival probability",
    legend.title = "PRS group",
    legend.labs = c("High", "Low/Mid")
)

png(here("tgct", "figures", "tgct_relapse_survival_prs_2level.png"), width = 7, height = 6.5, units = "in", res = 300)
print(p_group2)
dev.off()

# ========================
# SUBANALYSES - SEMINOMA
# ========================
# 1 = seminoma, 2/3 = nonseminoma
ss$histology_group <- factor(
    case_when(
        ss$testicular_tumor_histology == 1 ~ "Seminoma",
        ss$testicular_tumor_histology %in% c(2, 3) ~ "Nonseminoma",
        TRUE ~ NA_character_
    ),
    levels = c("Seminoma", "Nonseminoma")
)
table(ss$histology_group, useNA = "ifany")

fit_histology <- survfit(surv_obj ~ histology_group, data = ss)

p_histology <- ggsurvplot(
    fit_histology,
    data = ss,
    risk.table = TRUE,
    conf.int = TRUE,
    pval = TRUE,
    xlab = "Time",
    ylab = "Relapse-free survival probability",
    legend.title = "Histology",
    legend.labs = c("Seminoma", "Nonseminoma")
)

png(here("tgct", "figures", "tgct_relapse_survival_seminoma.png"), width = 7, height = 6.5, units = "in", res = 300)
print(p_histology)
dev.off()


