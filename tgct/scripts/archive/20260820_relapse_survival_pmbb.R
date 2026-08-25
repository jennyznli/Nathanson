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
    # select(ktid, empi, dob, stage_at_diagnosis,
    #        relapse_at_stage_1, time_to_event,
    #        date_of_diagnosis, date_last_contact,
    #        date1_of_relapse, date2_of_relapse,
    #        adjuvant_therapy_s1, histology_group,
    #        GenoSource3,
    #        PRS26_Z, PRS26_Effect, prs_2026,
    #        PRS21_Z, PRS21_Effect,
    #        PRS26_) |>
    mutate(
        PRS21_Z = as.numeric(PRS21_Z),
        PRS21_Effect = as.numeric(PRS21_Effect),
        PRS26_Z = as.numeric(PRS26_Z),
        PRS26_Effect = as.numeric(PRS26_Effect)
    )
dim(ss2)
# 827

# ========================
# TERTILE CONCORDANCE
# ========================
### WITHIN A YEAR, Z vs. EFFECT ###
table(Effect = ss2$PRS26_Effect_tertile, Z = ss2$PRS26_Z_tertile)
# Z
# Effect Low Mid High
# Low  267  10    0
# Mid   10 255   11
# High   0  11  265

table(Effect = ss2$PRS21_Effect_tertile, Z = ss2$PRS21_Z_tertile)
# Z
# Effect Low Mid High
# Low  266  11    0
# Mid   11 254   11
# High   0  11  265


### 2021 vs 2026 - WITHIN SAME METHOD ###
table(Effect = ss2$PRS21_Z_tertile, PRS2026 = ss2$PRS26_Z_tertile)
# PRS2026
# Effect Low Mid High
# Low  174  85   18
# Mid   83 105   88
# High  20  86  170

table(Effect = ss2$PRS21_Effect_tertile, PRS2026 = ss2$PRS26_Effect_tertile)
# PRS2026
# Effect Low Mid High
# Low  168  89   20
# Mid   86  99   91
# High  23  88  165

# these are slightly in agreement ...?

# % exact agreement (landed in the same Low/Mid/High bucket)
mean(ss2$PRS26_Effect_tertile == ss2$PRS26_Z_tertile, na.rm = TRUE)
# 0.95
mean(ss2$PRS21_Effect_tertile == ss2$PRS21_Z_tertile, na.rm = TRUE)
# 0.946924

mean(ss2$PRS26_Z_tertile == ss2$PRS21_Z_tertile, na.rm = TRUE)
# 0.67
mean(ss2$PRS26_Effect_tertile == ss2$PRS21_Effect_tertile, na.rm = TRUE)
# 0.6868198

# ========================
# DISTRIBUTION PLOT
# ========================
dist <- ss2 %>%
    filter(stage_at_diagnosis %in% c(1, 2, 3, 4, 5)) %>%
    mutate(stage_label = case_when(
        stage_at_diagnosis == 1 ~ "Stage I",
        stage_at_diagnosis == 5 ~ "Stage IS",
        stage_at_diagnosis == 2 ~ "Stage II",
        stage_at_diagnosis == 3 ~ "Stage III",
        stage_at_diagnosis == 4 ~ "Stage II/III or Metastatic"
    ),
    stage_label = factor(stage_label,
                         levels = c("Stage I", "Stage IS", "Stage II",
                                    "Stage III", "Stage II/III or Metastatic")))
table(dist$stage_at_diagnosis)
# 1   2   3   4   5   6  99
# 368 106 103  66   8   0   0

# ── Distribution plot (density) ───────────────────────────────────────────────
p <- ggplot(dist, aes(x = PRS26_Effect, color = stage_label)) +
    geom_density(linewidth = 1) +
    labs(
        x     = "PRS (PRS26_Effect)",
        y     = "Density",
        color = "Stage",
        title = "PRS Distribution by Stage at Diagnosis"
    ) +
    theme_minimal(base_size = 14)
png(here("tgct", "figures", "tgct_pmbb_distribution.png"), width = 10, height = 6, units = "in", res = 300)
print(p)
dev.off()

# ── Box plot for visual comparison ────────────────────────────────────────────
p <- ggplot(dist, aes(x = stage_label, y = PRS26_Effect, fill = stage_label)) +
    geom_boxplot(alpha = 0.7) +
    labs(
        x     = "Stage at Diagnosis",
        y     = "PRS",
        title = "PRS by Stage at Diagnosis"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 20, hjust = 1))

png(here("tgct", "figures", "tgct_pmbb_stage_prs.png"), width = 5, height = 6, units = "in", res = 300)
print(p)
dev.off()

# ── Summary table ─────────────────────────────────────────────────────────────
dist %>%
    group_by(stage_label) %>%
    summarise(
        n        = n(),
        mean_prs = round(mean(PRS26_Effect), 4),
        sd_prs   = round(sd(PRS26_Effect), 4),
        median   = round(median(PRS26_Effect), 4),
        min      = round(min(PRS26_Effect), 4),
        max      = round(max(PRS26_Effect), 4)
    )
# stage_label                    n mean_prs sd_prs median   min   max
# <fct>                      <int>    <dbl>  <dbl>  <dbl> <dbl> <dbl>
#     1 Stage I                      368     4.03   1.16   4.01 0.612  7.24
# 2 Stage IS                       8     4.33   1.19   4.79 2.40   5.85
# 3 Stage II                     106     4.03   1.28   3.87 1.12   7.88
# 4 Stage III                    103     4.07   1.20   4.00 1.04   6.75
# 5 Stage II/III or Metastatic    66     3.72   1.09   3.72 0.819  5.77

# ========================
# OVERALL SUMMARY
# ========================
summary(ss2$PRS26_Effect)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.2653  3.1055  3.8810  3.8823  4.6822  7.9593

summary(ss2$PRS26_Z)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 9.688 128.881 165.520 165.109 202.602 361.455

table(ss2$adjuvant_therapy_s1, useNA = "ifany")
# 0    1   99 <NA>
# 315  137    6  543

table(ss2$histology_group, useNA = "ifany")
# Seminoma Nonseminoma        <NA>
#     360         484         157

table(ss2$relapse_at_stage_1, useNA = "ifany")
# 0    1 <NA>
# 443  100  458

# filter out by all except for dates of relapse (not required)
ss3 <- ss2 |> filter(!is.na(histology_group), !is.na(relapse_at_stage_1),
                     !is.na(relapse_at_stage_1), !is.na(time_to_event),
                     !is.na(adjuvant_therapy_s1), !is.na(date_of_diagnosis),
                     !is.na(date_last_contact))
dim(ss3)
# 377

table(ss3$relapse_at_stage_1)
# 0   1
# 313  64

# ========================
# PRS DISTRIBUTION BY TREATMENT AND HISTOLOGY
# ========================
# data: must have PRS26_Effect/PRS26_Z, histology_group, adjuvant_therapy_s1
# prs_type: "effect" or "z"
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

# TODO: `therapy` is never defined anywhere in this script -- these calls
# will error as-is. Left disabled until it's clear which subset/object this
# was meant to reference (ss3? survival? something filtered further?).
# plot_prs_distribution(therapy, prs_type = "effect")
# plot_prs_distribution(therapy, prs_type = "z")

# ========================
# INITIAL TIME DISTRIBUTION
# ========================
# distribution of time to event, overlaid by censored vs relapse
survival <- ss3 |>
    mutate(status = if_else(as.character(relapse_at_stage_1) == "1", "Relapse", "Censored"))

p_dist <- ggplot(survival, aes(x = time_to_event, fill = status)) +
    geom_histogram(position = "identity", alpha = 0.5, binwidth = 1) +
    labs(
        title = "Distribution of Time to Event",
        x = "Time to event (years)",
        y = "Count",
        fill = NULL
    ) +
    theme_bw()
png(here("tgct", "figures", "tgct_pmbb_time_to_event_distribution.png"), width = 5, height = 4, units = "in", res = 300)
print(p_dist)
dev.off()

# ========================
# SURVIVAL - overall test
# ========================
# turn event into 1 - no relapse, 2 - relapse
# 1 is censored, 2 is event
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
png(here("tgct", "figures", "tgct_pmbb_relapse_survival_overall.png"), width = 7, height = 6.5, units = "in", res = 300)
print(p)
dev.off()

# ========================
# PRS-stratified survival plot
# ========================
# data: any subset of survival (must have time_to_event, relapse_at_stage_1, and the chosen prs_var)
# breaks: fixed 2-value quantile vector, e.g. effect_breaks/z_breaks/breaks_2026
#         from the "PRS TERTILES" section above -- always supplied by the
#         caller so tertile cutpoints stay fixed across every plot/subgroup
#         instead of being redefined per call
# prs_type: "effect" (PRS26_Effect), "z" (PRS26_Z), "2026" (prs_2026),
#           "effect21" (PRS21_Effect), or "z21" (PRS21_Z)
# coding: NA (default) -> tertiles (Low/Mid/High)
#         "binary_high" -> High vs Low/Mid
#         "binary_low"  -> Low vs Mid/High
# cohort_label: optional string folded into the output filename, e.g. "seminoma_adjuvant"
plot_prs_survival <- function(data, breaks, prs_type = c("effect", "z", "2026", "effect21", "z21"), coding = NA,
                               cohort_label = NULL, outdir = here("tgct", "figures")) {
    prs_type <- match.arg(prs_type)
    prs_var <- switch(prs_type,
        effect   = "PRS26_Effect",
        z        = "PRS26_Z",
        "2026"   = "prs_2026",
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

    # build Surv() directly in the formula (rather than a separate surv_obj
    # variable) so everything the formula needs lives inside `data` --
    # survminer's internals re-evaluate pieces of the formula outside this
    # function's scope, and a standalone local surv_obj isn't visible there
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
plot_prs_survival(survival, effect_breaks, prs_type = "effect", coding = NA)
plot_prs_survival(survival, effect_breaks, prs_type = "effect", coding = "binary_high")
plot_prs_survival(survival, effect_breaks, prs_type = "effect", coding = "binary_low")

plot_prs_survival(survival, z_breaks, prs_type = "z", coding = NA)
plot_prs_survival(survival, z_breaks, prs_type = "z", coding = "binary_high")
plot_prs_survival(survival, z_breaks, prs_type = "z", coding = "binary_low")

plot_prs_survival(survival, breaks_2026, prs_type = "2026", coding = NA)
plot_prs_survival(survival, breaks_2026, prs_type = "2026", coding = "binary_high")
plot_prs_survival(survival, breaks_2026, prs_type = "2026", coding = "binary_low")

plot_prs_survival(survival, effect21_breaks, prs_type = "effect21", coding = NA)
plot_prs_survival(survival, effect21_breaks, prs_type = "effect21", coding = "binary_high")
plot_prs_survival(survival, effect21_breaks, prs_type = "effect21", coding = "binary_low")

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
# RESTRICT DOWN TO SEMINOMA, ADJUVANT - PRS Stratified
# ========================
table(survival$adjuvant_therapy_s1)
# 0   1  99
# 289 111   1

sem_adj <- survival |> filter(histology_group == "Seminoma", adjuvant_therapy_s1 == "1")
dim(sem_adj)
# 64

plot_prs_survival(sem_adj, z_breaks, prs_type = "z", coding = NA, cohort_label = "seminoma_adjuvant")

# ========================
# RESTRICT DOWN TO SEMINOMA, SURVEILLANCE - PRS Stratified
# ========================
sem_surv <- survival |> filter(histology_group == "Seminoma", adjuvant_therapy_s1 == "0")
dim(sem_surv)
# 163

plot_prs_survival(sem_surv, z_breaks, prs_type = "z", coding = NA, cohort_label = "seminoma_surveillance")

# same PRS26_Z Low vs Mid/High split, using the tertiles fixed above
survival$prs_tertile_z <- cut(survival$PRS26_Z, breaks = c(-Inf, z_breaks, Inf), labels = c("Low", "Mid", "High"))
survival$prs_group_low_z <- if_else(survival$prs_tertile_z == "Low", "Low", "Mid/High")

# events/censoring driving the tail of the "Low" curve (risk table drops 12 -> 4 between t=20 and t=30)
survival |>
    filter(prs_group_low_z == "Low", time_to_event > 15) |>
    select(ktid, time_to_event, relapse_at_stage_1, PRS26_Z, PRS26_Effect,
           date_of_diagnosis, date1_of_relapse, date_last_contact,
           testicular_tumor_histology, adjuvant_therapy_s1) |>
    arrange(time_to_event)

x <- survival |>
    filter(prs_group_low_z == "Low", relapse_at_stage_1 == "1") |>
    select(ktid, time_to_event, PRS26_Z, date_of_diagnosis, date1_of_relapse,
           date_last_contact, testicular_tumor_histology, adjuvant_therapy_s1) |>
    arrange(time_to_event)

