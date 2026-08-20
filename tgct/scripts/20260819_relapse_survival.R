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
master <- read_xlsx(here("tgct", "ss", "20260818_tgct_master_cleaned.xlsx"))
# master <- read_excel(here("tgct", "ss", "20260730_tgct_master_w_import_info_jenny.xlsx"), sheet = "data", na = "NA")

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
# 1000

# ========================
# PREP DATA
# ========================
# first pass: join new PRS by ktid
ss2 <- ss |> left_join(master, by = c("ktid" = "KTID"))
table(ss2$GenoSource3)
# what if i just filtered down to PMBB ...

matched_ktid <- ss2 |> filter(!is.na(PRS_Effect))
missing_ktid <- ss2 |> filter(is.na(PRS_Effect)) |> select(all_of(names(ss)))
dim(missing_ktid)
# 202

# sanity check -- duplications
missing_ktid |> count(hup_mrn) |> filter(n > 1)
master |> count(MRN) |> filter(n > 1)
# all just NA values

# second pass: only the missing rows, joined by MRN instead
matched_mrn <- missing_ktid |> left_join(master, by = c("hup_mrn" = "MRN"), na_matches = "never")

ss3 <- bind_rows(matched_ktid, matched_mrn)
dim(ss3)
# 1184

# who's still missing after both passes?
still_missing <- ss3 |> filter(is.na(PRS_Effect))
dim(still_missing)
# 184
# View(still_missing)
# most of these don't have PRS from before, so good

# make numeric and filter out those w/o PRS
ss4 <- ss3 |> mutate(
        PRS_Z = as.numeric(PRS_Z),
        PRS_Effect = as.numeric(PRS_Effect)
    ) |>
    filter(GenoSource3 == "PMBB") |>
    filter(!is.na(PRS_Effect), !is.na(PRS_Z)) |>
    select(histology_group, relapse_at_stage_1,
           time_to_event, adjuvant_therapy_s1,
           date_of_diagnosis, date_last_contact,
           date1_of_relapse, date2_of_relapse,
           PRS_Z, PRS_Effect, prs_2026)
dim(ss4)
# 827

# ========================
# PRS TERTILES (computed once, reused everywhere below)
# ========================
# computed on ss4 (everyone with a valid PRS score) rather than on the later
# complete-case survival cohort -- PRS is fixed regardless of who has full
# clinical follow-up, so tertile cutpoints shouldn't depend on which subset
# happens to have complete histology/relapse/date fields. fixed cutpoints
# also mean "Low/Mid/High" means the same PRS value in every plot and every
# subgroup below, instead of being recomputed (and redefined) per call
effect_breaks <- quantile(ss4$PRS_Effect, probs = c(1 / 3, 2 / 3), na.rm = TRUE)
z_breaks      <- quantile(ss4$PRS_Z,      probs = c(1 / 3, 2 / 3), na.rm = TRUE)
breaks_2026   <- quantile(ss4$prs_2026,   probs = c(1 / 3, 2 / 3), na.rm = TRUE)

effect_breaks
z_breaks
breaks_2026

# label each sample's tertile per PRS type -- carries through ss5/survival/
# sem_adj/sem_surv automatically since it's just a column on ss4
ss4 <- ss4 |> mutate(
    PRS_Effect_tertile = cut(PRS_Effect, breaks = c(-Inf, effect_breaks, Inf), labels = c("Low", "Mid", "High")),
    PRS_Z_tertile       = cut(PRS_Z,      breaks = c(-Inf, z_breaks, Inf),      labels = c("Low", "Mid", "High")),
    prs_2026_tertile    = cut(prs_2026,   breaks = c(-Inf, breaks_2026, Inf),   labels = c("Low", "Mid", "High"))
)

table(ss4$PRS_Effect_tertile, useNA = "ifany")
table(ss4$PRS_Z_tertile, useNA = "ifany")
table(ss4$prs_2026_tertile, useNA = "ifany")

# ========================
# TERTILE CONCORDANCE
# ========================
# do the three PRS methods agree on which risk group (Low/Mid/High) a sample
# falls in? PRS_Effect and PRS_Z come from the same genotypes just weighted
# differently, so those two should track closely; prs_2026 is a separate
# calculation, so it's the one worth actually scrutinizing here
table(Effect = ss4$PRS_Effect_tertile, Z       = ss4$PRS_Z_tertile)
# effect mostly agrees with Z
table(Effect = ss4$PRS_Effect_tertile, PRS2026 = ss4$prs_2026_tertile)
table(Z      = ss4$PRS_Z_tertile,      PRS2026 = ss4$prs_2026_tertile)
# these really don't agree with each other ...

# % exact agreement (landed in the same Low/Mid/High bucket)
mean(ss4$PRS_Effect_tertile == ss4$PRS_Z_tertile, na.rm = TRUE)
# 0.95
mean(ss4$PRS_Effect_tertile == ss4$prs_2026_tertile, na.rm = TRUE)
# 0.40
mean(ss4$PRS_Z_tertile == ss4$prs_2026_tertile, na.rm = TRUE)
# 0.40...

# weighted kappa: Low/Mid/High is ordinal, so a Low-vs-Mid mismatch should
# count as less disagreement than a Low-vs-High mismatch -- plain %
# agreement or an unweighted kappa treats both misses the same way
library(irr)
kappa2(ss4[, c("PRS_Effect_tertile", "PRS_Z_tertile")], weight = "squared")
kappa2(ss4[, c("PRS_Effect_tertile", "prs_2026_tertile")], weight = "squared")
kappa2(ss4[, c("PRS_Z_tertile", "prs_2026_tertile")], weight = "squared")

# how many samples land in the same tertile across all three at once?
ss4 <- ss4 |> mutate(
    tertile_agree_all3 = (PRS_Effect_tertile == PRS_Z_tertile) & (PRS_Z_tertile == prs_2026_tertile)
)
table(ss4$tertile_agree_all3, useNA = "ifany")
mean(ss4$tertile_agree_all3, na.rm = TRUE)

# ========================
# DISTRIBUTION PLOT
# ========================
ss4 <- ss4 %>%
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
table(ss4$stage_at_diagnosis)
# 1   2   3   4   5   6  99
# 428 128 115  74   8   0   0

# ── Distribution plot (density) ───────────────────────────────────────────────
p <- ggplot(ss3, aes(x = PRS_Effect, color = stage_label)) +
    geom_density(linewidth = 1) +
    labs(
        x     = "PRS (PRS_Effect)",
        y     = "Density",
        color = "Stage",
        title = "PRS Distribution by Stage at Diagnosis"
    ) +
    theme_minimal(base_size = 14)
png(here("tgct", "figures", "tgct_distribution.png"), width = 10, height = 6, units = "in", res = 300)
print(p)
dev.off()

# ── Box plot for visual comparison ────────────────────────────────────────────
p <- ggplot(ss3, aes(x = stage_label, y = PRS_Effect, fill = stage_label)) +
    geom_boxplot(alpha = 0.7) +
    labs(
        x     = "Stage at Diagnosis",
        y     = "PRS",
        title = "PRS by Stage at Diagnosis"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 20, hjust = 1))

png(here("tgct", "figures", "tgct_stage_prs.png"), width = 5, height = 6, units = "in", res = 300)
print(p)
dev.off()

# ── Summary table ─────────────────────────────────────────────────────────────
ss3 %>%
    group_by(stage_label) %>%
    summarise(
        n        = n(),
        mean_prs = round(mean(PRS_Effect), 4),
        sd_prs   = round(sd(PRS_Effect), 4),
        median   = round(median(PRS_Effect), 4),
        min      = round(min(PRS_Effect), 4),
        max      = round(max(PRS_Effect), 4)
    )
# stage_label                    n mean_prs sd_prs median   min   max
# <fct>                      <int>    <dbl>  <dbl>  <dbl> <dbl> <dbl>
#     1 Stage I                      429     3.85   1.14   3.84 0.409  6.97
# 2 Stage IS                       8     4.23   1.21   4.70 2.13   5.34
# 3 Stage II                     128     3.89   1.31   3.82 0.879  7.96
# 4 Stage III                    115     3.88   1.21   3.84 1.05   6.81
# 5 Stage II/III or Metastatic    74     3.56   1.10   3.52 0.713  5.5

# ========================
# OVERALL SUMMARY
# ========================
summary(ss4$PRS_Effect)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.2653  3.1055  3.8810  3.8823  4.6822  7.9593

summary(ss4$PRS_Z)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 9.688 128.881 165.520 165.109 202.602 361.455

table(ss4$adjuvant_therapy_s1, useNA = "ifany")
# 0    1   99 <NA>
# 315  137    6  543

table(ss4$histology_group, useNA = "ifany")
# Seminoma Nonseminoma        <NA>
#     360         484         157

table(ss4$relapse_at_stage_1, useNA = "ifany")
# 0    1 <NA>
# 443  100  458


colSums(is.na(examine))
# histology_group  relapse_at_stage_1       time_to_event adjuvant_therapy_s1
# 157                 458                 468                 543
# date_of_diagnosis   date_last_contact    date1_of_relapse    date2_of_relapse
# 132                 135                 905                1000

# filter out by all except for dates of relapse (not required)
ss5 <- ss4 |> filter(!is.na(histology_group), !is.na(relapse_at_stage_1),
                     !is.na(relapse_at_stage_1), !is.na(time_to_event),
                     !is.na(adjuvant_therapy_s1), !is.na(date_of_diagnosis),
                     !is.na(date_last_contact))
dim(ss5)
# 438

table(ss5$relapse_at_stage_1)
# 0   1
# 364  74

# ========================
# PRS DISTRIBUTION BY TREATMENT AND HISTOLOGY
# ========================
# data: must have PRS_Effect/PRS_Z, histology_group, adjuvant_therapy_s1
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
    png(file.path(outdir, filename), width = 5, height = 4, units = "in", res = 300)
    print(p)
    dev.off()

    invisible(list(plot = p, summary = summary_table))
}

# TODO: `therapy` is never defined anywhere in this script -- these calls
# will error as-is. Left disabled until it's clear which subset/object this
# was meant to reference (ss5? survival? something filtered further?).
# plot_prs_distribution(therapy, prs_type = "effect")
# plot_prs_distribution(therapy, prs_type = "z")

# ========================
# INITIAL TIME DISTRIBUTION
# ========================
# distribution of time to event, overlaid by censored vs relapse
survival <- ss5 |>
    mutate(status = if_else(as.character(relapse_at_stage_1) == "1", "Relapse", "Censored"))
dim(survival)
# 438

table(survival$status)
# Censored  Relapse
# 364       74

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
png(here("tgct", "figures", "tgct_relapse_survival_overall.png"), width = 7, height = 6.5, units = "in", res = 300)
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
# prs_type: "effect" (PRS_Effect), "z" (PRS_Z), or "2026" (prs_2026)
# coding: NA (default) -> tertiles (Low/Mid/High)
#         "binary_high" -> High vs Low/Mid
#         "binary_low"  -> Low vs Mid/High
# cohort_label: optional string folded into the output filename, e.g. "seminoma_adjuvant"
plot_prs_survival <- function(data, breaks, prs_type = c("effect", "z", "2026"), coding = NA,
                               cohort_label = NULL, outdir = here("tgct", "figures")) {
    prs_type <- match.arg(prs_type)
    prs_var <- switch(prs_type,
        effect = "PRS_Effect",
        z      = "PRS_Z",
        "2026" = "prs_2026"
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

# same PRS_Z Low vs Mid/High split, using the tertiles fixed above
survival$prs_tertile_z <- cut(survival$PRS_Z, breaks = c(-Inf, z_breaks, Inf), labels = c("Low", "Mid", "High"))
survival$prs_group_low_z <- if_else(survival$prs_tertile_z == "Low", "Low", "Mid/High")

# events/censoring driving the tail of the "Low" curve (risk table drops 12 -> 4 between t=20 and t=30)
survival |>
    filter(prs_group_low_z == "Low", time_to_event > 15) |>
    select(ktid, time_to_event, relapse_at_stage_1, PRS_Z, PRS_Effect,
           date_of_diagnosis, date1_of_relapse, date_last_contact,
           testicular_tumor_histology, adjuvant_therapy_s1) |>
    arrange(time_to_event)

x <- survival |>
    filter(prs_group_low_z == "Low", relapse_at_stage_1 == "1") |>
    select(ktid, time_to_event, PRS_Z, date_of_diagnosis, date1_of_relapse,
           date_last_contact, testicular_tumor_histology, adjuvant_therapy_s1) |>
    arrange(time_to_event)

