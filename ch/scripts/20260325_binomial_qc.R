# ========================
# Binomial for germline QC
# ========================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
library(data.table, quietly = T)

# ========================
# READ IN DATA
# ========================
all_ch <- read_excel(file.path("ch", "data", "ch_wl_art.xlsx"))
dim(all_ch)
# 230

cov <- read.csv(file.path("ch", "data", "pmbb_brca12_cov_df.csv"))
all_ids <- read.csv(file.path("ch", "data", "ch_psm_matched4_case_control_ids.csv"))$x
cov <- cov %>% filter(person_id %in% all_ids)
dim(cov)
# 3004

# ========================
# BINOMIAL GERMLINE FILTER
# ========================
# run binomial test on every variant
# H0: VAF = 0.5 (germline heterozygous)
# p < 0.01 → consistent with germline → flag
all_ch_clean <- all_ch_clean %>%
    mutate(
        binom_p = mapply(function(alt, depth) {
            if (is.na(alt) || is.na(depth) || depth == 0) return(NA_real_)
            binom.test(alt, depth, p = 0.5, alternative = "two.sided")$p.value
        }, Sample.AltDepth, Sample.Depth),
        binom_fail = !is.na(binom_p) & binom_p < 0.01
    )

cat("Variants failing binomial test:", sum(all_ch_clean$binom_fail, na.rm = TRUE), "\n")
cat("Variants passing:", sum(!all_ch_clean$binom_fail, na.rm = TRUE), "\n")

# ========================
# TET2 + CBL: age association rescue
# for variants failing binomial test, check if including them
# improves or weakens age association — if improves, exempt them
# ========================
tet2_cbl_mis <- all_ch_clean %>%
    filter(Gene %in% c("TET2", "CBL"), wl.exception, binom_fail)

# for each failing variant site, test age association with vs without
rescue_test <- function(vkey) {
    # with the failing variant included
    with_var <- all_ch_clean %>%
        filter(variant_id == vkey) %>%
        distinct(Sample.ID) %>%
        mutate(carrier = 1)

    df_with <- all_samples_cov %>%
        left_join(with_var, by = "Sample.ID") %>%
        mutate(carrier = ifelse(is.na(carrier), 0, carrier))

    # without (passing only)
    without_var <- all_ch_clean %>%
        filter(variant_id == vkey, !binom_fail) %>%
        distinct(Sample.ID) %>%
        mutate(carrier = 1)

    df_without <- all_samples_cov %>%
        left_join(without_var, by = "Sample.ID") %>%
        mutate(carrier = ifelse(is.na(carrier), 0, carrier))

    tryCatch({
        m_with    <- glm(carrier ~ Age, data = df_with,    family = binomial())
        m_without <- glm(carrier ~ Age, data = df_without, family = binomial())

        p_with    <- coef(summary(m_with))["Age",    "Pr(>|z|)"]
        p_without <- coef(summary(m_without))["Age", "Pr(>|z|)"]

        data.frame(
            p_with    = p_with,
            p_without = p_without,
            # improved = adding failing variants strengthens age association
            rescue    = p_with < p_without
        )
    }, error = function(e) {
        data.frame(p_with = NA_real_, p_without = NA_real_, rescue = FALSE)
    })
}

tet2_cbl_rescue <- tet2_cbl_mis %>%
    distinct(variant_id, Gene) %>%
    rowwise() %>%
    mutate(res = list(rescue_test(variant_id))) %>%
    unnest(res) %>%
    ungroup()

print(tet2_cbl_rescue %>% select(variant_id, Gene, p_with, p_without, rescue))

# variants to exempt from binomial removal
rescue_ids <- tet2_cbl_rescue %>% filter(rescue) %>% pull(variant_id)
cat("TET2/CBL variants rescued:", length(rescue_ids), "\n")

# ========================
# ALL OTHER GENES: remove sites where ALL variants fail binomial
# ========================
site_summary <- all_ch_clean %>%
    filter(!is.na(variant_id)) %>%
    group_by(variant_id, Gene) %>%
    summarise(
        n_total    = n(),
        n_fail     = sum(binom_fail, na.rm = TRUE),
        all_fail   = all(binom_fail, na.rm = TRUE),
        .groups    = "drop"
    ) %>%
    filter(n_total >= FREQ_THRESHOLD)  # only check recurrent variants

cat("\nRecurrent variant sites where ALL fail binomial test:\n")
site_summary %>% filter(all_fail) %>% arrange(desc(n_total)) %>% print()

germline_ids <- site_summary %>%
    filter(all_fail) %>%
    pull(variant_id)

cat("Variant sites flagged as recurrent germline:", length(germline_ids), "\n")

# ========================
# APPLY FILTERS
# ========================
all_ch_clean <- all_ch_clean %>%
    mutate(
        binom_exempt = variant_id %in% rescue_ids,
        germline_flag = case_when(
            variant_id %in% germline_ids                    ~ "germline_recurrent",
            binom_fail & !binom_exempt                      ~ "germline_possible",
            TRUE                                            ~ "pass"
        )
    )

# remove: recurrent germline sites + individual failing variants (except rescued)
all_ch_final <- all_ch_clean %>%
    filter(germline_flag == "pass" | binom_exempt)

cat("\nVariants before germline filter:", nrow(all_ch_clean), "\n")
cat("Variants after germline filter: ", nrow(all_ch_final), "\n")
cat("Removed:", nrow(all_ch_clean) - nrow(all_ch_final), "\n")

# write out
write.csv(all_ch_final, file.path("ch", "data", "ch_wl_final_clean.csv"), row.names = FALSE)
write.csv(
    all_ch_clean %>% filter(germline_flag != "pass" & !binom_exempt),
    file.path("ch", "data", "ch_germline_flagged.csv"),
    row.names = FALSE
)

