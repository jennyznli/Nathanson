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
master <- read_xlsx(here("tgct", "ss", "20260820_tgct_pmbb_cases.xlsx"))

ss <- read_xlsx(here("tgct", "ss", "PRSTCStageIRelapse-AllPatientsAndData_DATA_2026-08-17_0843.xlsx")) |>
    merge_duplicates("ktid") |>
    mutate(
        empi = as.character(empi),
        hup_mrn = as.character(hup_mrn),
        date1_of_relapse = as.Date(date1_of_relapse, format = "%Y-%m-%d"),
        date2_of_relapse = as.Date(date2_of_relapse, format = "%Y-%m-%d"),
        date_of_diagnosis = as.Date(date_of_diagnosis, format = "%Y-%m-%d"),
        date_last_contact = as.Date(date_last_contact, format = "%Y-%m-%d"),
        relapse_at_stage_1 = as.factor(relapse_at_stage_1),
        stage_at_diagnosis = as.factor(stage_at_diagnosis),
        adjuvant_therapy_s1 = as.factor(adjuvant_therapy_s1),
        prs_2026 = as.numeric(prs_2026)
    )

# time-to-event: time to relapse for events, time to last contact for censored
# (previously used date1_of_relapse for everyone, which is NA for the ~393
# patients who never relapsed -- that silently dropped all censored patients
# from the survival fit instead of treating them as censored)
ss$time_to_event <- as.numeric(difftime(
    if_else(ss$relapse_at_stage_1 == "1", ss$date1_of_relapse, ss$date_last_contact),
    ss$date_of_diagnosis,
    units = "days"
)) / 365.25

# 1 = seminoma, 2/3 = nonseminoma
ss$histology_group <- factor(
    case_when(
        ss$testicular_tumor_histology == 1 ~ "Seminoma",
        ss$testicular_tumor_histology %in% c(2, 3) ~ "Nonseminoma",
        TRUE ~ NA_character_
    ),
    levels = c("Seminoma", "Nonseminoma")
)

dim(ss)
# 1184
dim(master)
# 827

# ========================
# PREP DATA
# ========================
# first pass: join new PRS by ktid
ss2 <- ss |> left_join(master, by = c("ktid" = "KTID")) |>
    filter(!is.na(PRS26_Effect)) |>
    mutate(
        PRS21_Z = as.numeric(PRS21_Z),
        PRS21_Effect = as.numeric(PRS21_Effect),
        PRS26_Z = as.numeric(PRS26_Z),
        PRS26_Effect = as.numeric(PRS26_Effect)
    )
dim(ss2)
# 1003

# ========================
# PRS TERTILES (computed once, reused everywhere below)
# ========================
effect26_breaks   <- quantile(ss2$PRS26_Effect, probs = c(1 / 3, 2 / 3), na.rm = TRUE)
z26_breaks        <- quantile(ss2$PRS26_Z,      probs = c(1 / 3, 2 / 3), na.rm = TRUE)
effect21_breaks <- quantile(ss2$PRS21_Effect, probs = c(1 / 3, 2 / 3), na.rm = TRUE)
z21_breaks      <- quantile(ss2$PRS21_Z,      probs = c(1 / 3, 2 / 3), na.rm = TRUE)

# label each sample's tertile per PRS type
ss2 <- ss2 |> mutate(
    PRS26_Effect_tertile = cut(PRS26_Effect, breaks = c(-Inf, effect26_breaks, Inf), labels = c("Low", "Mid", "High")),
    PRS26_Z_tertile       = cut(PRS26_Z,      breaks = c(-Inf, z26_breaks, Inf),      labels = c("Low", "Mid", "High")),
    PRS21_Effect_tertile = cut(PRS21_Effect, breaks = c(-Inf, effect21_breaks, Inf), labels = c("Low", "Mid", "High")),
    PRS21_Z_tertile       = cut(PRS21_Z,      breaks = c(-Inf, z21_breaks, Inf),      labels = c("Low", "Mid", "High"))
)

table(ss2$PRS26_Effect_tertile, useNA = "ifany")
table(ss2$PRS26_Z_tertile, useNA = "ifany")
table(ss2$PRS21_Effect_tertile, useNA = "ifany")
table(ss2$PRS21_Z_tertile, useNA = "ifany")
# Low  Mid High
# 335  334  334

# ========================
# TERTILE CONCORDANCE
# ========================
### WITHIN A YEAR, Z vs. EFFECT ###
table(PRS26_Effect = ss2$PRS26_Effect_tertile, PRS26_Z = ss2$PRS26_Z_tertile)
# PRS26_Z
# PRS26_Effect Low Mid High
# Low  324  11    0
# Mid   11 311   12
# High   0  12  322

table(PRS21_Effect = ss2$PRS21_Effect_tertile, PRS21_Z = ss2$PRS21_Z_tertile)
# PRS21_Z
# PRS21_Effect Low Mid High
# Low  323  12    0
# Mid   12 309   13
# High   0  13  321

### WITHIN SAME METHOD, 2021 vs. 2026 ###
table(PRS21_Z = ss2$PRS21_Z_tertile, PRS26_Z = ss2$PRS26_Z_tertile)
# PRS26_Z
# PRS21_Z Low Mid High
# Low  203 108   24
# Mid  104 118  112
# High  28 108  198

table(PRS21_Effect = ss2$PRS21_Effect_tertile, PRS26_Effect = ss2$PRS26_Effect_tertile)
# PRS26_Effect
# PRS21_Effect Low Mid High
# Low  206 105   24
# Mid  104 126  104
# High  25 103  206

# % exact agreement
mean(ss2$PRS26_Effect_tertile == ss2$PRS26_Z_tertile, na.rm = TRUE)
# 0.95
mean(ss2$PRS21_Effect_tertile == ss2$PRS21_Z_tertile, na.rm = TRUE)
# 0.9501496

mean(ss2$PRS26_Z_tertile == ss2$PRS21_Z_tertile, na.rm = TRUE)
# 0.5174477
mean(ss2$PRS26_Effect_tertile == ss2$PRS21_Effect_tertile, na.rm = TRUE)
# 0.5363908

# ========================
# INITIAL TIME DISTRIBUTION
# ========================
# distribution of time to event, overlaid by censored vs relapse
survival <- ss2 |>
    mutate(status = if_else(as.character(relapse_at_stage_1) == "1", "Relapse", "Censored"))

table(survival$status)
# Censored  Relapse
# 444      100

p_dist <- ggplot(survival, aes(x = time_to_event, fill = status)) +
    geom_histogram(position = "identity", alpha = 0.5, binwidth = 1) +
    labs(
        title = "Distribution of Time to Event",
        x = "Time to event (years)",
        y = "Count",
        fill = NULL
    ) +
    theme_bw()
png(here("tgct", "figures", "tgct_time_to_event_distribution.png"), width = 5, height = 4, units = "in", res = 300)
print(p_dist)
dev.off()

# ========================
# PRS DISTRIBUTION BY TREATMENT AND HISTOLOGY
# ========================
# prs_type: "effect" or "z"
# just for 2026
plot_prs_distribution <- function(data, prs_type = c("effect", "z"), cohort_label = NULL,
                                  outdir = here("tgct", "figures")) {
    prs_type <- match.arg(prs_type)
    prs_var <- if (prs_type == "effect") "PRS26_Effect" else "PRS26_Z"

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
    png(file.path(outdir, filename), width = 5, height = 4, units = "in", res = 300)
    print(p)
    dev.off()

    invisible(list(plot = p, summary = summary_table))
}

plot_prs_distribution(survival, prs_type = "effect")
plot_prs_distribution(survival, prs_type = "z")

# ========================
# SURVIVAL - overall test
# ========================
surv_obj <- Surv(
    time = survival$time_to_event,
    event = as.numeric(as.character(survival$relapse_at_stage_1))
)

# Overall relapse-free survival
fit <- survfit(surv_obj ~ 1, data = survival)

p <- ggsurvplot(
    fit,
    data = survival,
    risk.table = TRUE,
    conf.int = TRUE,
    xlab = "Time",
    ylab = "Relapse-free survival probability"
)
png(here("tgct", "figures", "tgct_relapse_survival_overall.png"), width = 7, height = 6.5, units = "in", res = 300)
print(p)
dev.off()

# ========================
# PRS-stratified survival plot
# ========================
# data: any subset of survival (must have time_to_event, relapse_at_stage_1, and the chosen prs_var)
# breaks: fixed 2-value quantile vector, e.g. z26_breaks
#         from the "PRS TERTILES" section above -- always supplied by the
#         caller so tertile cutpoints stay fixed across every plot/subgroup
#         instead of being redefined per call
# prs_type: "effect26" (PRS26_Effect), "z26" (PRS26_Z),
#           "effect21" (PRS21_Effect), "z21" (PRS21_Z)
# coding: NA (default) -> tertiles (Low/Mid/High)
#         "binary_high" -> High vs Low/Mid
#         "binary_low"  -> Low vs Mid/High
# cohort_label: optional string folded into the output filename, e.g. "seminoma_adjuvant"
plot_prs_survival <- function(data, breaks, prs_type = c("effect26", "z26", "effect21", "z21"), coding = NA,
                               cohort_label = NULL, outdir = here("tgct", "figures")) {
    prs_type <- match.arg(prs_type)
    prs_var <- switch(prs_type,
        effect26   = "PRS26_Effect",
        z26        = "PRS26_Z",
        effect21 = "PRS21_Effect",
        z21      = "PRS21_Z"
    )

    data[[prs_var]] <- as.numeric(data[[prs_var]])
    data <- data |> filter(!is.na(.data[[prs_var]]), !is.na(time_to_event), !is.na(relapse_at_stage_1))

    data$prs_group <- cut(
        data[[prs_var]],
        breaks = c(-Inf, breaks, Inf),
        labels = c("Low", "Mid", "High")
    )

    if (is.na(coding)) {
        legend_labs <- c("Low", "Mid", "High")
        legend_title <- "PRS tertile"
        file_suffix <- "tertile"
    } else if (coding == "binary_high") {
        legend_labs <- c("High", "Low/Mid")
        data$prs_group <- factor(
            if_else(data$prs_group == "High", "High", "Low/Mid"),
            levels = legend_labs
        )
        legend_title <- "PRS group"
        file_suffix <- "2level_high"
    } else if (coding == "binary_low") {
        legend_labs <- c("Low", "Mid/High")
        data$prs_group <- factor(
            if_else(data$prs_group == "Low", "Low", "Mid/High"),
            levels = legend_labs
        )
        legend_title <- "PRS group"
        file_suffix <- "2level_low"
    } else {
        stop("coding must be NA, 'binary_high', or 'binary_low'")
    }

    print(table(data$prs_group))

    fit <- survfit(Surv(time_to_event, as.numeric(as.character(relapse_at_stage_1))) ~ prs_group, data = data)

    p <- ggsurvplot(
        fit,
        data = data,
        risk.table = TRUE,
        ylim = c(0.5, 1.0),
        xlim = c(0, 20),
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
# 2026 effect
plot_prs_survival(survival, effect26_breaks, prs_type = "effect26", coding = NA)
plot_prs_survival(survival, effect26_breaks, prs_type = "effect26", coding = "binary_high")
plot_prs_survival(survival, effect26_breaks, prs_type = "effect26", coding = "binary_low")

# 2026 z
plot_prs_survival(survival, z26_breaks, prs_type = "z26", coding = NA)
plot_prs_survival(survival, z26_breaks, prs_type = "z26", coding = "binary_high")
plot_prs_survival(survival, z26_breaks, prs_type = "z26", coding = "binary_low")

# 2021 effect
plot_prs_survival(survival, effect21_breaks, prs_type = "effect21", coding = NA)
plot_prs_survival(survival, effect21_breaks, prs_type = "effect21", coding = "binary_high")
plot_prs_survival(survival, effect21_breaks, prs_type = "effect21", coding = "binary_low")

# 2021 z
plot_prs_survival(survival, z21_breaks, prs_type = "z21", coding = NA)
plot_prs_survival(survival, z21_breaks, prs_type = "z21", coding = "binary_high")
plot_prs_survival(survival, z21_breaks, prs_type = "z21", coding = "binary_low")

# ========================
# SUBANALYSES - SEMINOMA
# ========================
fit_histology <- survfit(surv_obj ~ histology_group, data = survival)

p_histology <- ggsurvplot(
    fit_histology,
    data = survival,
    risk.table = TRUE,
    conf.int = FALSE,
    pval = TRUE,
    xlab = "Time",
    ylab = "Relapse-free survival probability",
    legend.title = "Histology",
    legend.labs = c("Seminoma", "Nonseminoma")
)

png(here("tgct", "figures", "tgct_relapse_survival_seminoma.png"), width = 7, height = 6.5, units = "in", res = 300)
print(p_histology)
dev.off()


table(survival$adjuvant_therapy_s1)
# 0   1  99
# 316 137   6

# ========================
# RESTRICT DOWN TO SEMINOMA, ADJUVANT - PRS Stratified
# ========================
sem_adj <- survival |> filter(histology_group == "Seminoma", adjuvant_therapy_s1 == "1")
dim(sem_adj)
# 75

plot_prs_survival(sem_adj, z26_breaks, prs_type = "z26", coding = NA, cohort_label = "seminoma_adjuvant")

# ========================
# RESTRICT DOWN TO SEMINOMA, SURVEILLANCE - PRS Stratified
# ========================
sem_surv <- survival |> filter(histology_group == "Seminoma", adjuvant_therapy_s1 == "0")
dim(sem_surv)
# 166

plot_prs_survival(sem_surv, z26_breaks, prs_type = "z26", coding = NA, cohort_label = "seminoma_surveillance")

# ========================
# RESTRICT DOWN TO NONSEMINOMA, ADJUVANT - PRS Stratified
# ========================
nonsem_adj <- survival |> filter(histology_group == "Nonseminoma", adjuvant_therapy_s1 == "1")
dim(nonsem_adj)
# 62

plot_prs_survival(nonsem_adj, z26_breaks, prs_type = "z26", coding = NA, cohort_label = "nonseminoma_adjuvant")

# ========================
# RESTRICT DOWN TO SEMINOMA, SURVEILLANCE - PRS Stratified
# ========================
nonsem_surv <- survival |> filter(histology_group == "Nonseminoma", adjuvant_therapy_s1 == "0")
dim(nonsem_surv)
# 149

plot_prs_survival(nonsem_surv, z26_breaks, prs_type = "z26", coding = NA, cohort_label = "nonseminoma_surveillance")


