# ============================================================
# CNA ANALYSIS - BATCH PROCESSING
# ============================================================
library(here)
source("R/config.R")
library(survival)
library(survminer)

# ========================
# LOAD DATA & CONVERT
# ========================
master <- read_xlsx(here("tgct", "ss", "20260813_tgct_master.xlsx"))
ss_s1 <- read_xlsx(here("tgct", "ss", "prstcstageirelapse-stageipatients_data_2026-04-14_1446.xlsx"))
ss <- read_xlsx(here("tgct", "ss", "PRSTCStageIRelapse-AllPatientsAndData_DATA_2026-08-17_0843.xlsx"))
dim(ss)
# 1186

# need to join the PRS somehow...

# ========================
# PREP DATA
# ========================
not_kt <- master |> filter(KTID == "TGCT_Not_KT")
kt <- master |> filter(KTID != "TGCT_Not_KT")
dim(kt)
# 984
dim(not_kt)
# 108

ss <- ss |> left_join(kt, by = c("ktid" = "KTID")) |>
filter(ktid != "KT0002101")
# this guys has 20 year relasep... KT0002101

ss$relapse_at_stage_1 <- as.factor(ss$relapse_at_stage_1)
ss$adjuvant_chemotherapy_s1 <- as.factor(ss$adjuvant_chemotherapy_s1)
ss$PRS_Z <- as.numeric(ss$PRS_Z)
ss$PRS_Effect <- as.numeric(ss$PRS_Effect)

ss <- ss |> mutate(
    date1_of_relapse = as.Date(date1_of_relapse, format = "%m/%d/%y"),
    date2_of_relapse = as.Date(date2_of_relapse, format = "%m/%d/%y"),
    date_of_diagnosis = as.Date(date_of_diagnosis, format = "%m/%d/%y"),
    date_last_contact = as.Date(date_last_contact, format = "%m/%d/%y")
)

# time-to-event: time to relapse for events, time to last contact for censored
# (previously used date1_of_relapse for everyone, which is NA for the ~393
# patients who never relapsed -- that silently dropped all censored patients
# from the survival fit instead of treating them as censored)
ss$time_to_event <- as.numeric(difftime(
    if_else(as.character(ss$relapse_at_stage_1) == "1", ss$date1_of_relapse, ss$date_last_contact),
    ss$date_of_diagnosis,
    units = "days"
)) / 365.25

ss <- ss |> filter(!is.na(testicular_tumor_histology), !is.na(relapse_at_stage_1), !is.na(PRS_Effect),
                   !is.na(time_to_event))
dim(ss)
# 401

summary(ss$PRS_Effect)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.7992  3.0610  3.8335  3.8234  4.6220  6.5941

summary(ss$PRS_Z)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#   25.48  126.50  163.00  162.15  198.50  288.98      15

# lowkey what ...
sum(is.na(ss$adjuvant_chemotherapy_s1))
# 342

# ========================
# SURVIVAL - overall
# ========================
table(ss$relapse_at_stage_1)
# 0   1
# 339   62

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
# FUNCTION: PRS-stratified survival plot
# ========================
# data: any subset of ss (must have time_to_event, relapse_at_stage_1, PRS_Effect, PRS_Z)
# prs_type: "effect" or "z"
# coding: NA (default) -> tertiles (Low/Mid/High)
#         "binary_high" -> High vs Low/Mid
#         "binary_low"  -> Low vs Mid/High
# cohort_label: optional string folded into the output filename, e.g. "seminoma_adjuvant"
plot_prs_survival <- function(data, prs_type = c("effect", "z"), coding = NA, cohort_label = NULL,
                              outdir = here("tgct", "figures")) {
    prs_type <- match.arg(prs_type)
    prs_var <- if (prs_type == "effect") "PRS_Effect" else "PRS_Z"

    data <- data |> filter(!is.na(.data[[prs_var]]), !is.na(time_to_event), !is.na(relapse_at_stage_1))

    tertile_breaks <- quantile(data[[prs_var]], probs = c(1 / 3, 2 / 3), na.rm = TRUE)
    data$prs_group <- cut(
        data[[prs_var]],
        breaks = c(-Inf, tertile_breaks, Inf),
        labels = c("Low", "Mid", "High")
    )

    if (is.na(coding)) {
        legend_labs <- c("Low", "Mid", "High")
        legend_title <- "PRS tertile"
        file_suffix <- "tertile"
    } else if (coding == "binary_high") {
        data$prs_group <- if_else(data$prs_group == "High", "High", "Low/Mid")
        legend_labs <- c("High", "Low/Mid")
        legend_title <- "PRS group"
        file_suffix <- "2level_high"
    } else if (coding == "binary_low") {
        data$prs_group <- if_else(data$prs_group == "Low", "Low", "Mid/High")
        legend_labs <- c("Low", "Mid/High")
        legend_title <- "PRS group"
        file_suffix <- "2level_low"
    } else {
        stop("coding must be NA, 'binary_high', or 'binary_low'")
    }

    print(table(data$prs_group))

    # build Surv() directly in the formula (rather than a separate surv_obj
    # variable) so everything the formula needs lives inside `data` --
    # survminer's internals re-evaluate pieces of the formula outside this
    # function's scope, and a standalone local surv_obj isn't visible there
    fit <- survfit(Surv(time_to_event, as.numeric(as.character(relapse_at_stage_1))) ~ prs_group, data = data)

    p <- ggsurvplot(
        fit,
        data = data,
        risk.table = TRUE,
        conf.int = TRUE,
        pval = TRUE,
        xlab = "Time",
        ylab = "Relapse-free survival probability",
        legend.title = legend_title,
        legend.labs = legend_labs
    )

    filename <- paste0(
        "tgct_relapse_survival_prs_", prs_type,
        if (!is.null(cohort_label)) paste0("_", cohort_label) else "",
        "_", file_suffix, ".png"
    )
    png(file.path(outdir, filename), width = 7, height = 6.5, units = "in", res = 300)
    print(p)
    dev.off()

    invisible(p)
}

# ========================
# STRATIFIED BY PRS
# ========================
plot_prs_survival(ss, prs_type = "effect", coding = NA)
plot_prs_survival(ss, prs_type = "effect", coding = "binary_high")
plot_prs_survival(ss, prs_type = "effect", coding = "binary_low")

plot_prs_survival(ss, prs_type = "z", coding = NA)
plot_prs_survival(ss, prs_type = "z", coding = "binary_high")
plot_prs_survival(ss, prs_type = "z", coding = "binary_low")

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
# Seminoma Nonseminoma
# 228         173

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

# ========================
# PRS DISTRIBUTION BY TREATMENT AND HISTOLOGY
# ========================
# data: any subset of ss (must have PRS_Effect/PRS_Z, histology_group, adjuvant_therapy_s1)
# prs_type: "effect" or "z"
plot_prs_distribution <- function(data, prs_type = c("effect", "z"), cohort_label = NULL,
                                  outdir = here("tgct", "figures")) {
    prs_type <- match.arg(prs_type)
    prs_var <- if (prs_type == "effect") "PRS_Effect" else "PRS_Z"

    data <- data |>
        filter(!is.na(.data[[prs_var]]), !is.na(histology_group), adjuvant_therapy_s1 %in% c("0", "1")) |>
        mutate(
            treatment = if_else(adjuvant_therapy_s1 == "1", "Adjuvant Therapy", "Surveillance"),
            histology = factor(
                if_else(histology_group == "Nonseminoma", "Non-Seminoma", "Seminoma"),
                levels = c("Non-Seminoma", "Seminoma")
            )
        )

    summary_table <- data |>
        group_by(histology, treatment) |>
        summarise(
            n = n(),
            mean = mean(.data[[prs_var]]),
            sd = sd(.data[[prs_var]]),
            median = median(.data[[prs_var]]),
            q1 = quantile(.data[[prs_var]], 0.25),
            q3 = quantile(.data[[prs_var]], 0.75),
            min = min(.data[[prs_var]]),
            max = max(.data[[prs_var]]),
            .groups = "drop"
        )
    print(summary_table)

    # Wilcoxon rank-sum test: Adjuvant Therapy vs Surveillance, within each histology facet
    p <- ggplot(data, aes(x = treatment, y = .data[[prs_var]], fill = treatment)) +
        geom_boxplot() +
        facet_wrap(~histology) +
        stat_compare_means(method = "wilcox.test", label = "p.format", label.y.npc = 0.97) +
        labs(
            title = "PRS Distribution: Adjuvant vs Surveillance by Histology",
            x = NULL,
            y = paste0("PRS (", prs_var, ")")
        ) +
        theme_bw() +
        theme(legend.position = "none")

    filename <- paste0(
        "tgct_prs_distribution_adjuvant_vs_surveillance",
        if (!is.null(cohort_label)) paste0("_", cohort_label) else "",
        "_", prs_type, ".png"
    )
    png(file.path(outdir, filename), width = 7, height = 6, units = "in", res = 300)
    print(p)
    dev.off()

    invisible(list(plot = p, summary = summary_table))
}

plot_prs_distribution(ss, prs_type = "effect")
plot_prs_distribution(ss, prs_type = "z")

# ========================
# RESTRICT DOWN TO SEMINOMA, ADJUVANT - PRS Stratified
# ========================
table(ss$adjuvant_therapy_s1)
# 0   1  99
# 289 111   1

sem_adj <- ss |> filter(histology_group == "Seminoma", adjuvant_therapy_s1 == "1")
dim(sem_adj)
# 64

plot_prs_survival(sem_adj, prs_type = "z", coding = NA, cohort_label = "seminoma_adjuvant")

# ========================
# RESTRICT DOWN TO SEMINOMA, SURVEILLANCE - PRS Stratified
# ========================
sem_surv <- ss |> filter(histology_group == "Seminoma", adjuvant_therapy_s1 == "0")
dim(sem_surv)
# 163

plot_prs_survival(sem_surv, prs_type = "z", coding = NA, cohort_label = "seminoma_surveillance")

# recreate the same PRS_Z Low vs Mid/High split plot_prs_survival() used
z_breaks <- quantile(ss$PRS_Z, probs = c(1/3, 2/3), na.rm = TRUE)
ss$prs_tertile_z <- cut(ss$PRS_Z, breaks = c(-Inf, z_breaks, Inf), labels = c("Low", "Mid", "High"))
ss$prs_group_low_z <- if_else(ss$prs_tertile_z == "Low", "Low", "Mid/High")

# events/censoring driving the tail of the "Low" curve (risk table drops 12 -> 4 between t=20 and t=30)
ss |>
    filter(prs_group_low_z == "Low", time_to_event > 15) |>
    select(ktid, time_to_event, relapse_at_stage_1, PRS_Z, PRS_Effect,
           date_of_diagnosis, date1_of_relapse, date_last_contact,
           testicular_tumor_histology, adjuvant_therapy_s1) |>
    arrange(time_to_event)


x <- ss |>
    filter(prs_group_low_z == "Low", relapse_at_stage_1 == "1") |>
    select(ktid, time_to_event, PRS_Z, date_of_diagnosis, date1_of_relapse,
           date_last_contact, testicular_tumor_histology, adjuvant_therapy_s1) |>
    arrange(time_to_event)



