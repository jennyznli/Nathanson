library(readr)
library(readxl)
library(dplyr)
library(fuzzyjoin)
library(stringr)
library(survival)
library(survminer)
library(dplyr)
library(lubridate)

Stage_I_relapse <- read_csv("Downloads/PRSTCStageIRelapse-StageIPatients_DATA_2026-08-10_1212.csv")

#No score_sums
Stage_I_relapse <- Stage_I_relapse %>%
  filter(!is.na(score_sum))

Seminoma_relapse <- Stage_I_relapse %>%
  filter(testicular_tumor_histology == 1)

Nonseminoma_relapse <- Stage_I_relapse %>%
  filter(testicular_tumor_histology == 2 | testicular_tumor_histology == 3)

# ── Calculate tertile cutpoints from the full dataset ─────────────────────────
cutpoints <- quantile(Stage_I_relapse$score_sum, probs = c(1/3, 2/3), na.rm = TRUE)

# ── Helper function using fixed cutpoints ─────────────────────────────────────
run_km <- function(data, title_text, cuts) {
  data <- data %>%
    mutate(
      event = ifelse(!is.na(date1_of_relapse), 1, 0),
      end_date = case_when(
        !is.na(date1_of_relapse) ~ date1_of_relapse,
        TRUE                     ~ date_last_contact
      ),
      time_years = as.numeric(difftime(end_date, date_of_diagnosis, units = "days")) / 365.25,
      score_tertile = case_when(
        score_sum <= cuts[1] ~ "T1 (Low)",
        score_sum <= cuts[2] ~ "T2 (Mid)",
        TRUE                 ~ "T3 (High)"
      ),
      score_tertile = factor(score_tertile,
                             levels = c("T1 (Low)", "T2 (Mid)", "T3 (High)"))
    ) %>%
    filter(!is.na(time_years) & time_years >= 0)
  
  km_fit <- survfit(Surv(time_years, event) ~ score_tertile, data = data)
  tertile_counts <- data %>% count(score_tertile)
  tertile_labels <- paste0(tertile_counts$score_tertile, " (n=", tertile_counts$n, ")")
  
  p <- ggsurvplot(
    km_fit,
    data        = data,
    pval        = TRUE,
    pval.method = TRUE,
    risk.table  = FALSE,
    conf.int    = FALSE,
    ylim        = c(0.5, 1.0),
    xlim        = c(0, 20),
    break.time.by = 5,
    palette     = c("#2E86AB", "#A23B72", "#F18F01"),
    xlab        = "Time (years)",
    ylab        = "Relapse-Free Survival Probability",
    title       = title_text,
    legend.title = "Score Tertile",
    legend.labs  = tertile_labels,
    ggtheme      = theme_minimal(base_size = 14)
  )
  print(p)
  
  cat("\n--- Pairwise comparisons:", title_text, "---\n")
  pw <- pairwise_survdiff(Surv(time_years, event) ~ score_tertile,
                          data = data, p.adjust.method = "BH")
  print(pw)
  
  return(pw)
}

# ── Now all subgroups use the same cutpoints ──────────────────────────────────
run_km(Stage_I_relapse, "RFS by Score Sum Tertile (Stage I All)", cutpoints)
run_km(Seminoma_relapse, "RFS by Score Sum Tertile (Stage I Seminoma)", cutpoints)
run_km(Nonseminoma_relapse, "RFS by Score Sum Tertile (Stage I Non-Seminoma)", cutpoints)

Seminoma_surv <- Seminoma_relapse %>% filter(adjuvant_therapy_s1 == 0)
Nonseminoma_surv <- Nonseminoma_relapse %>% filter(adjuvant_therapy_s1 == 0)

run_km(Seminoma_surv, "RFS by Score Sum Tertile (Stage I Seminoma, Surveillance Only)", cutpoints)
run_km(Nonseminoma_surv, "RFS by Score Sum Tertile (Stage I Non-Seminoma, Surveillance Only)", cutpoints)

Seminoma_adj <- Seminoma_relapse %>% filter(adjuvant_therapy_s1 == 1)
Nonseminoma_adj <- Nonseminoma_relapse %>% filter(adjuvant_therapy_s1 == 1)

run_km(Seminoma_adj, "RFS by Score Sum Tertile (Stage I Seminoma, Adjuvant Therapy)", cutpoints)
run_km(Nonseminoma_adj, "RFS by Score Sum Tertile (Stage I Non-Seminoma, Adjuvant Therapy)", cutpoints)

All_stages <- read_csv("Downloads/PRSTCStageIRelapse-AllPatientsAndData_DATA_2026-08-10_1250.csv")

#No score_sums
All_stages <- All_stages %>%
  filter(!is.na(score_sum))

library(ggplot2)
library(dplyr)

# ── Filter and label stages ───────────────────────────────────────────────────
All_stages_filtered <- All_stages %>%
  filter(stage_at_diagnosis %in% c(1, 2, 3, 4, 5)) %>%
  filter(!is.na(score_sum)) %>%
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

# ── Distribution plot (density) ───────────────────────────────────────────────
ggplot(All_stages_filtered, aes(x = score_sum, color = stage_label)) +
  geom_density(linewidth = 1) +
  labs(
    x     = "PRS (score_sum)",
    y     = "Density",
    color = "Stage",
    title = "PRS Distribution by Stage at Diagnosis"
  ) +
  theme_minimal(base_size = 14)

# ── Box plot for visual comparison ────────────────────────────────────────────
ggplot(All_stages_filtered, aes(x = stage_label, y = score_sum, fill = stage_label)) +
  geom_boxplot(alpha = 0.7) +
  labs(
    x     = "Stage at Diagnosis",
    y     = "PRS (score_sum)",
    title = "PRS by Stage at Diagnosis"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))

# ── Summary table ─────────────────────────────────────────────────────────────
All_stages_filtered %>%
  group_by(stage_label) %>%
  summarise(
    n        = n(),
    mean_prs = round(mean(score_sum), 4),
    sd_prs   = round(sd(score_sum), 4),
    median   = round(median(score_sum), 4),
    min      = round(min(score_sum), 4),
    max      = round(max(score_sum), 4)
  )

# ── Statistical tests ────────────────────────────────────────────────────────

# Kruskal-Wallis (non-parametric, overall test across all groups)
kruskal.test(score_sum ~ stage_label, data = All_stages_filtered)

# Pairwise Wilcoxon tests with BH correction
pairwise.wilcox.test(All_stages_filtered$score_sum,
                     All_stages_filtered$stage_label,
                     p.adjust.method = "BH")

# If you prefer parametric (ANOVA + Tukey):
aov_result <- aov(score_sum ~ stage_label, data = All_stages_filtered)
summary(aov_result)
TukeyHSD(aov_result)

Stage_I_relapse %>%
  mutate(
    treatment = ifelse(adjuvant_therapy_s1 == 0,
                       "Surveillance", "Adjuvant Therapy"),
    histology = case_when(
      testicular_tumor_histology == 1 ~ "Seminoma",
      testicular_tumor_histology %in% c(2, 3) ~ "Non-Seminoma"
    )
  ) %>%
  filter(!is.na(score_sum), !is.na(histology)) %>%
  ggplot(aes(x = treatment, y = score_sum, fill = treatment)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~ histology) +
  labs(
    x     = "",
    y     = "PRS (score_sum)",
    title = "PRS Distribution: Adjuvant vs Surveillance by Histology"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

Stage_I_relapse %>%
  mutate(treatment = ifelse(adjuvant_therapy_s1 == 0,
                            "Surveillance", "Adjuvant Therapy")) %>%
  filter(!is.na(score_sum)) %>%
  group_by(treatment) %>%
  summarise(
    n      = n(),
    mean   = round(mean(score_sum), 4),
    sd     = round(sd(score_sum), 4),
    median = round(median(score_sum), 4),
    q1     = round(quantile(score_sum, 0.25), 4),
    q3     = round(quantile(score_sum, 0.75), 4),
    min    = round(min(score_sum), 4),
    max    = round(max(score_sum), 4)
  )

Stage_I_relapse %>%
  mutate(
    treatment = ifelse(adjuvant_therapy_s1 == 0,
                       "Surveillance", "Adjuvant Therapy"),
    histology = case_when(
      testicular_tumor_histology == 1 ~ "Seminoma",
      testicular_tumor_histology %in% c(2, 3) ~ "Non-Seminoma"
    )
  ) %>%
  filter(!is.na(score_sum), !is.na(histology)) %>%
  group_by(histology, treatment) %>%
  summarise(
    n      = n(),
    mean   = round(mean(score_sum), 4),
    sd     = round(sd(score_sum), 4),
    median = round(median(score_sum), 4),
    q1     = round(quantile(score_sum, 0.25), 4),
    q3     = round(quantile(score_sum, 0.75), 4),
    min    = round(min(score_sum), 4),
    max    = round(max(score_sum), 4),
    .groups = "drop"
  )
