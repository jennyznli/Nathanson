# ========================
# Age association for artifacts
# ========================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
library(data.table, quietly = T)
library(ggVennDiagram)
library(UpSetR)

# ========================
# SET GLOBALS
# ========================
GNOMAD_THRESHOLD  <- 0.005
FREQ_THRESHOLD    <- 4
CLUSTER_THRESHOLD <- 50
MIN_AD_THRESHOLDS <- 3:8

AGE_SIG      <- 0.10
FREEZE_SIG   <- 0.05
GERMLINE_SIG <- 0.05

n_f2 <- 760
n_f3 <- 2244

extra_covs <- "Smoke_History + Sequenced_gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6"

# ========================
# READ IN DATA
# ========================
all_ch <- read_excel(file.path("ch", "data", "ch_seq_wl_vars.xlsx"))
all_ch$gnomAD.MAX_AF <- as.numeric(all_ch$gnomAD.MAX_AF)
all_ch <- all_ch %>% mutate(
        VAF_Strata = case_when(
            Sample.AltFrac >= 0.10              ~ "10+",
            Sample.AltFrac >= 0.02              ~ "2-10",
            TRUE                                ~ NA_character_
        )
    )

cov <- read.csv(file.path("ch", "data", "pmbb_brca12_cov_df.csv"), row.names = 1)
all_ids <- read.csv(file.path("ch", "data", "ch_psm_matched4_case_control_ids.csv"))$x
cov <- cov %>% filter(person_id %in% all_ids)
table(cov$Batch)
# 1    2
# 760 2244
dim(cov)
# 3004

all_ch <- all_ch %>% left_join(cov %>% select(person_id, Sample_age), by = c("Sample.ID" = "person_id"))
dim(all_ch)
# 774 84

# ========================
# 1. FLAG - HIGH FREQ VARIANTS
# ========================
variant_counts <- all_ch %>%
    group_by(variant_id) %>%
    summarise(
        n_carriers          = n_distinct(Sample.ID),
        batch_1             = n_distinct(Sample.ID[Batch == 1]),
        batch_2             = n_distinct(Sample.ID[Batch == 2]),
        Gene                = first(Gene),
        Chr                 = first(Chr),
        Start               = first(Start),
        variant_category    = first(variant_category),
        Variant.Consequence = first(Variant.Consequence),
        .groups             = "drop"
    ) %>%
    arrange(desc(n_carriers)) %>%
    mutate(flag_high_freq = as.integer(n_carriers >= FREQ_THRESHOLD))

dim(variant_counts)
# 481

high_freq <- variant_counts %>% filter(flag_high_freq == 1)
high_freq_vars <- unique(high_freq$variant_id)
high_freq_var <- all_ch %>% filter(variant_id %in% high_freq_vars)

write_xlsx(
    list(summary = variant_counts, removed = high_freq_var),
    file.path("ch", "data", "ch_high_freq_removed.xlsx")
)

# ========================
# 2. FLAG - FREEZE ENRICHMENT
# ========================
high_freq2 <- high_freq %>%
    mutate(
        p_freeze = mapply(function(b1, b2) {
            fisher.test(matrix(
                c(b1, n_f2 - b1,
                  b2, n_f3 - b2),
                nrow = 2
            ))$p.value
        }, batch_1, batch_2),
        p_freeze_adj = p.adjust(p_freeze, method = "BH"),
        flag_freeze_enriched = p_freeze_adj < FREEZE_SIG
    )

freeze_enriched <- high_freq2 %>% filter(flag_freeze_enriched)
freeze_enriched_var <- freeze_enriched %>% pull(variant_id)
freeze_removed_var <- all_ch %>% filter(variant_id %in% freeze_enriched_var)

write_xlsx(
    list(summary = freeze_enriched, removed = freeze_removed_var),
    file.path("ch", "data", "ch_freeze_enriched_removed.xlsx")
)

# ========================
# 3. FLAG - AGE ASSOCIATION
# ========================
run_age_tests <- function(data, variants) {
    variants %>%
        rowwise() %>%
        mutate(res = list(artifact_test_age(variant_id, data, extra_covs = extra_covs))) %>%
        unnest(res) %>%
        ungroup() %>%
        mutate(
            p_age_adj = p.adjust(p_age, method = "BH"),
            flag_age_associated = !is.na(p_age) & !separation & p_age < AGE_SIG
        )
}

# test age association per variant category
artifact_test_age <- function(vkey, vars, extra_covs = NULL) {
    carriers <- vars %>%
        filter(variant_id == vkey) %>%
        distinct(Sample.ID) %>%
        mutate(carrier = 1)

    df <- cov %>%
        left_join(carriers, by = c("person_id" = "Sample.ID")) %>%
        mutate(carrier = ifelse(is.na(carrier), 0, carrier))

    base    <- "carrier ~ Sample_age + as.factor(Batch)"
    formula <- if (!is.null(extra_covs)) {
        as.formula(paste(base, "+", extra_covs))
    } else {
        as.formula(base)
    }

    tryCatch({
        m <- withCallingHandlers(
            glm(formula, data = df, family = binomial()),
            warning = function(w) {
                if (grepl("fitted probabilities numerically 0 or 1", conditionMessage(w)))
                    invokeRestart("muffleWarning")
            }
        )
        separation <- any(fitted(m) %in% c(0, 1)) | !m$converged
        s  <- coef(summary(m))
        ci <- confint.default(m)
        data.frame(
            OR         = exp(coef(m)["Sample_age"]),
            CI_lo      = exp(ci["Sample_age", 1]),
            CI_hi      = exp(ci["Sample_age", 2]),
            p_age      = s["Sample_age", "Pr(>|z|)"],
            separation = separation
        )
    }, error = function(e) {
        warning(sprintf("Model failed for %s: %s", vkey, conditionMessage(e)))
        data.frame(OR = NA_real_, CI_lo = NA_real_, CI_hi = NA_real_,
                   p_age = NA_real_, separation = NA)
    })
}

age_association <- high_freq2 %>%
    rowwise() %>%
    mutate(res = list(artifact_test_age(variant_id, all_ch, extra_covs = extra_covs))) %>%
    unnest(res) %>% ungroup() %>%
    mutate(
        p_age_adj = p.adjust(p_age, method = "BH"),
        flag_age_associated = !is.na(p_age) & !separation & p_age < AGE_SIG
    )

not_age_associated <- age_association %>% filter(!flag_age_associated)
not_age_associated_var <- not_age_associated %>% pull(variant_id)
age_removed_var <- all_ch %>% filter(variant_id %in% not_age_associated_var)

write_xlsx(
    list(summary = not_age_associated, removed = age_removed_var),
    file.path("ch", "data", "ch_age_association_removed.xlsx")
)

### TEST IF OVERALL ASSOCIATION IMPROVES ###
all_ch2 <- all_ch %>%
    left_join(
        age_removed_var %>% select(Sample.ID, variant_id) %>% mutate(flag_not_age_associated_var = TRUE),
        by = c("Sample.ID", "variant_id")
    ) %>%
    mutate(flag_not_age_associated_var = ifelse(is.na(flag_not_age_associated_var), FALSE, flag_not_age_associated_var))

age_config <- list(
    list(name = "before_age_association", flags = c()),
    list(name = "after_age_association",  flags = c("flag_not_age_associated_var"))
)

age_res <- test_age_association(age_config, all_ch2, cov)
age_res
# name stratum                       flags n_carriers n_total       OR    CI_lo    CI_hi        p_age separation
# Sample_age...1 before_age_association     all                                    556    3004 1.013259 1.006386 1.020179 1.487781e-04      FALSE
# Sample_age...2  after_age_association     all flag_not_age_associated_var        413    3004 1.018569 1.010947 1.026250 1.582579e-06      FALSE

# ========================
# 4. FLAG - GERMLINE TESTING
# ========================
binomial_germline_test <- function(alt_count, total_depth, p_germline = 0.5) {
    if (is.na(alt_count) || is.na(total_depth) || total_depth == 0) {
        return(list(p_binom = NA_real_, flag_germline_ind = NA))
    }
    p <- binom.test(x = alt_count, n = total_depth,
                    p = p_germline, alternative = "two.sided")$p.value
    list(p_binom = p, flag_germline_ind = p >= GERMLINE_SIG)
}

# flag variants to test
germline_test_vars <- all_ch %>%
    filter(
        Gene %in% c("TET2", "CBL") |
            variant_id %in% high_freq2$variant_id |
            Sample.AltFrac > 0.30
    ) %>%
    distinct(variant_id, Sample.ID, Sample_age, Sample.AltDepth, Sample.Depth, Sample.AltFrac,
             gnomAD.MAX_AF, ProteinChange_1L, variant_category, Gene) %>%
    mutate(
        pmap_dfr(
            list(Sample.AltDepth, Sample.Depth),
            ~ binomial_germline_test(..1, ..2)
        )
    )

# variant category proportions failing binomial test
germline_summary <- germline_test_vars %>%
    group_by(variant_id) %>%
    summarise(
        n_tested          = sum(!is.na(p_binom)),
        n_fail_binom      = sum(flag_germline_ind, na.rm = TRUE),
        prop_fail_binom   = n_fail_binom / n_tested,
        flag_germline_var = prop_fail_binom > 0.50,
        .groups           = "drop"
    )
fail_majority_germline <- germline_summary %>% filter(flag_germline_var)

# from all vars tested, remove those individuals' variants failing binomial
# as well as entire variant categories whose majority failed
germline_fail_vars <- germline_test_vars %>%
    filter(variant_id %in% fail_majority_germline$variant_id) %>%
    filter(flag_germline_ind)

germline_removed_vars <- all_ch %>% semi_join(germline_fail_vars, by = c("Sample.ID", "variant_id"))

cat(sprintf("Variants tested for germline: %d\n", nrow(germline_summary)))   # 93
cat(sprintf("Individuals flagged germline: %d\n", length(unique(germline_fail_vars$Sample.ID)))) # 25
cat(sprintf("Variants flagged germline: %d\n", nrow(germline_removed_vars))) # 25

write_xlsx(
    list(summary = fail_majority_germline, removed = germline_removed_vars),
    file.path("ch", "data", "ch_germline_removed.xlsx")
)

### PLOT ###
# colored by TET2 missense
tet2_vars <- germline_test_vars %>% filter(Gene == "TET2")
p <- ggplot(tet2_vars,
                 aes(x = Sample_age, y = Sample.AltFrac,
                     color = ProteinChange_1L,
                     size  = !flag_germline_ind)) +
    geom_point(alpha = 0.85) +
    scale_color_manual(
        values = setNames(
            colorRampPalette(c("#4B0082", "#C0392B", "#E67E22", "#F1C40F"))(n_distinct(tet2_vars$ProteinChange_1L)),
            sort(unique(tet2_vars$ProteinChange_1L))
        ),
        name = "Protein change"
    ) +
    scale_size_manual(
        values = c("TRUE" = 2, "FALSE" = 4),
        labels = c("TRUE" = "p \u2265 0.05", "FALSE" = "p < 0.05"),
        name   = "Binomial test"
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Baseline age", y = "Variant allele fraction (VAF %)") +
    theme_minimal(base_size = 11) +
    theme(
        panel.grid.minor = element_blank(),
        legend.position  = "right",
        legend.title     = element_text(face = "bold", size = 9),
        legend.text      = element_text(size = 9)
    ) +
    guides(
        color = guide_legend(order = 2, override.aes = list(size = 3)),
        size  = guide_legend(order = 1)
    )

ggsave(file.path("ch", "figures", "germline_vaf_age_tet2_missense.pdf"),
       p, width = 9, height = 6)

# colored by type
p <- ggplot(germline_test_vars,
                 aes(x = Sample_age, y = Sample.AltFrac,
                     color = variant_category,
                     size  = !flag_germline_ind)) +
    geom_point(alpha = 0.85) +
    scale_size_manual(
        values = c("TRUE" = 2, "FALSE" = 4),
        labels = c("TRUE" = "p \u2265 0.05", "FALSE" = "p < 0.05"),
        name   = "Binomial test"
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Baseline age", y = "Variant allele fraction (VAF %)") +
    theme_minimal(base_size = 11) +
    theme(
        panel.grid.minor = element_blank(),
        legend.position  = "right",
        legend.title     = element_text(face = "bold", size = 9),
        legend.text      = element_text(size = 9)
    ) +
    guides(
        color = guide_legend(order = 2, override.aes = list(size = 3)),
        size  = guide_legend(order = 1)
    )

ggsave(file.path("ch", "figures", "germline_vaf_age_category.pdf"), p, width = 9, height = 6)

# colored by gene
p <- ggplot(germline_test_vars,
            aes(x = Sample_age, y = Sample.AltFrac,
                color = Gene,
                size  = !flag_germline_ind)) +
    geom_point(alpha = 0.85) +
    scale_size_manual(
        values = c("TRUE" = 2, "FALSE" = 4),
        labels = c("TRUE" = "p \u2265 0.05", "FALSE" = "p < 0.05"),
        name   = "Binomial test"
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Baseline age", y = "Variant allele fraction (VAF %)") +
    theme_minimal(base_size = 11) +
    theme(
        panel.grid.minor = element_blank(),
        legend.position  = "right",
        legend.title     = element_text(face = "bold", size = 9),
        legend.text      = element_text(size = 9)
    ) +
    guides(
        color = guide_legend(order = 2, override.aes = list(size = 3)),
        size  = guide_legend(order = 1)
    )

ggsave(file.path("ch", "figures", "germline_vaf_age_gene.pdf"), p, width = 9, height = 6)

# colored by MAX AF
germline_test_vars$gnomAD.MAX_AF[germline_test_vars$gnomAD.MAX_AF == 0] <- NA
p <- ggplot(germline_test_vars,
            aes(x = Sample_age, y = Sample.AltFrac,
                color = gnomAD.MAX_AF,
                size  = !flag_germline_ind)) +
    geom_point(alpha = 0.85) +
    scale_color_gradient(low = "blue", high = "red", trans = "log10") +
    scale_size_manual(
        values = c("TRUE" = 2, "FALSE" = 4),
        labels = c("TRUE" = "p \u2265 0.05", "FALSE" = "p < 0.05"),
        name   = "Binomial test"
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Baseline age", y = "Variant allele fraction (VAF %)") +
    theme_minimal(base_size = 11) +
    theme(
        panel.grid.minor = element_blank(),
        legend.position  = "right",
        legend.title     = element_text(face = "bold", size = 9),
        legend.text      = element_text(size = 9)
    ) +
    guides(
        color = guide_colorbar(order = 2),
        size  = guide_legend(order = 1)
    )

ggsave(file.path("ch", "figures", "germline_vaf_age_gnomad.pdf"), p, width = 9, height = 6)

### TEST OVERALL AGE ASSOCIATION BEFORE/AFTER REMOVING
all_ch3 <- all_ch2 %>%
    left_join(
        germline_fail_vars %>% select(Sample.ID, variant_id) %>% mutate(flag_germline_var = TRUE),
        by = c("Sample.ID", "variant_id")
    ) %>%
    mutate(flag_germline_var = ifelse(is.na(flag_germline_var), FALSE, flag_germline_var))

germline_config <- list(
    list(name = "before_germline", flags = c()),
    list(name = "after_germline",  flags = c("flag_germline_var", "flag_not_age_associated_var"))
)

germline_res <- test_age_association(germline_config, all_ch3, cov)
germline_res
# name stratum                                          flags n_carriers n_total       OR    CI_lo    CI_hi        p_age separation
# Sample_age...1 before_germline     all                                                       556    3004 1.013259 1.006386 1.020179 1.487781e-04      FALSE
# Sample_age...2  after_germline     all flag_germline_var, flag_not_age_associated_var        393    3004 1.020680 1.012879 1.028540 1.700100e-07      FALSE

# ### TEST INDIVIDUAL VARIANTS - AGE ASSOCIATION BEFORE/AFTER REMOVING
# # variants to test = those that went through germline testing
# germline_test_variant_ids <- germline_test_vars %>% distinct(variant_id)
#
# # before
# germline_before_age <- germline_test_variant_ids %>%
#     rowwise() %>%
#     mutate(res = list(artifact_test_age(variant_id, all_ch2, extra_covs = extra_covs))) %>%
#     unnest(res) %>% ungroup() %>%
#     mutate(
#         p_age_adj           = p.adjust(p_age, method = "BH"),
#         flag_age_associated = !is.na(p_age) & !separation & p_age < AGE_SIG,
#         stage               = "before"
#     )
#
# # after
# germline_after_age <- germline_test_variant_ids %>%
#     rowwise() %>%
#     mutate(res = list(artifact_test_age(variant_id, all_ch2 %>% filter(!flag_germline_var), extra_covs = extra_covs))) %>%
#     unnest(res) %>% ungroup() %>%
#     mutate(
#         p_age_adj           = p.adjust(p_age, method = "BH"),
#         flag_age_associated = !is.na(p_age) & !separation & p_age < AGE_SIG,
#         stage               = "after"
#     )
#
# # combine and compare
# germline_age_comparison <- bind_rows(germline_before_age, germline_after_age) %>%
#     left_join(germline_summary %>% select(variant_id, flag_germline_var, prop_fail_binom),
#               by = "variant_id")
#
# # summary
# germline_age_comparison %>%
#     group_by(stage) %>%
#     summarise(
#         n_variants       = n(),
#         n_age_associated = sum(flag_age_associated, na.rm = TRUE),
#         prop_associated  = n_age_associated / n_variants,
#         .groups          = "drop"
#     )
#
# # specifically: which variants changed status before -> after?
# germline_age_comparison %>%
#     select(variant_id, stage, flag_age_associated, p_age, p_age_adj) %>%
#     pivot_wider(names_from = stage,
#                 values_from = c(flag_age_associated, p_age, p_age_adj)) %>%
#     filter(flag_age_associated_before != flag_age_associated_after)
#
# # how many failed to converge?
# germline_before_age %>%
#     filter(separation) %>%
#     nrow()
#
# germline_after_age %>%
#     filter(separation) %>%
#     nrow()
#
# # what's the carrier count for converging vs non-converging?
# germline_before_age %>%
#     left_join(variant_counts %>% select(variant_id, n_carriers), by = "variant_id") %>%
#     group_by(separation) %>%
#     summarise(
#         n          = n(),
#         median_carriers = median(n_carriers),
#         max_carriers    = max(n_carriers),
#         .groups = "drop"
#     )
#
# # separation     n median_carriers max_carriers
# # <lgl>      <int>           <dbl>        <int>
# #     1 FALSE         46               2          102
# # 2 TRUE          47               1            1

# ========================
# 5. FLAG - CLUSTERS
# ========================
all_ch_clustered <- all_ch3 %>%
    mutate(.row_idx = row_number()) %>%
    group_by(Sample.ID, Gene, Chr) %>%
    arrange(Sample.ID, Gene, Chr, Start) %>%
    mutate(
        near_prev  = c(FALSE, diff(Start) <= CLUSTER_THRESHOLD),
        in_cluster = near_prev | lead(near_prev, default = FALSE)
    ) %>%
    ungroup()

cluster_all <- all_ch_clustered %>%
    filter(in_cluster) %>%
    group_by(Sample.ID, Gene, Chr) %>%
    arrange(Start) %>%
    mutate(
        new_cluster = !near_prev,
        cluster_id  = paste(Sample.ID, Gene, Chr, cumsum(new_cluster), sep = "_")
    ) %>%
    ungroup() %>%
    select(-new_cluster, -.row_idx)

cluster_rep <- cluster_all %>%
    group_by(cluster_id) %>%
    arrange(!grepl("\\bC[0-9]", Existing.variation), Start) %>%
    slice(1) %>%
    ungroup()

cluster_not_rep_var <- cluster_all %>%
    filter(!(variant_id %in% cluster_rep$variant_id)) %>%
    pull(variant_id)

cluster_not_rep_var <- cluster_all %>%
    filter(!(variant_id %in% cluster_rep$variant_id)) %>%
    pull(variant_id)

# ========================
# JOIN BACK TO ALL_CH & VAF STRATA
# ========================
all_ch2 <- all_ch %>%
    mutate(
        flag_freeze_enriched    = variant_id %in% freeze_enriched_vars,
        flag_not_age_associated = variant_id %in% not_age_associated_var,
        flag_germline_ind       = paste(Sample.ID, variant_id) %in% paste(germline_ind_pairs$Sample.ID, germline_ind_pairs$variant_id),
        flag_cluster            = variant_id %in% cluster_not_rep_var
    )

for (thresh in MIN_AD_THRESHOLDS) {
    all_ch2[[paste0("minad", thresh)]] <- all_ch2$Sample.AltDepth < thresh
}

write_xlsx(all_ch2, file.path("ch", "data", "ch_seq_wl_art_flags_vars.xlsx"))

# ========================
# FLAG SUMMARY TABLE
# ========================
flag_cols <- c(
    "flag_freeze_enriched", "flag_not_age_associated",
    "flag_germline_ind", "flag_germline_var",
    "flag_cluster", "flag_gnomad",
    paste0("minad", MIN_AD_THRESHOLDS)
)

flag_summary <- do.call(rbind, lapply(flag_cols, function(f) {
    col      <- as.logical(all_ch2[[f]])
    batch_1  <- all_ch2$Batch == 1
    batch_2  <- all_ch2$Batch == 2
    n_total  <- nrow(all_ch2)
    n_f2     <- sum(batch_1)
    n_f3     <- sum(batch_2)

    data.frame(
        Flag          = f,
        N_Flagged     = sum(col, na.rm = TRUE),
        Pct_Overall   = round(100 * mean(col, na.rm = TRUE), 2),
        N_Freeze2     = sum(col & batch_1, na.rm = TRUE),
        Pct_Freeze2   = round(100 * sum(col & batch_1, na.rm = TRUE) / n_f2, 2),
        N_Freeze3     = sum(col & batch_2, na.rm = TRUE),
        Pct_Freeze3   = round(100 * sum(col & batch_2, na.rm = TRUE) / n_f3, 2),
        stringsAsFactors = FALSE
    )
}))

# Add total row at bottom for context
totals <- data.frame(
    Flag        = "TOTAL (rows)",
    N_Flagged   = nrow(all_ch2),
    Pct_Overall = 100.00,
    N_Freeze2   = sum(all_ch2$Batch == 1),
    Pct_Freeze2 = 100.00,
    N_Freeze3   = sum(all_ch2$Batch == 2),
    Pct_Freeze3 = 100.00,
    stringsAsFactors = FALSE
)
flag_summary <- rbind(flag_summary, totals)

print(flag_summary, row.names = FALSE)
# Flag N_Flagged Pct_Overall N_Freeze2 Pct_Freeze2 N_Freeze3 Pct_Freeze3
# flag_freeze_enriched       149       19.25         0        0.00       149       22.01
# flag_not_age_associated       188       24.29        10       10.31       178       26.29
# flag_germline_ind        18        2.33         3        3.09        15        2.22
# flag_germline_var         5        0.65         3        3.09         2        0.30
# flag_cluster       131       16.93         9        9.28       122       18.02
# flag_gnomad        81       10.47        28       28.87        53        7.83
# minad3         0        0.00         0        0.00         0        0.00
# minad4       321       41.47        38       39.18       283       41.80
# minad5       460       59.43        58       59.79       402       59.38
# minad6       521       67.31        68       70.10       453       66.91
# minad7       572       73.90        72       74.23       500       73.86
# minad8       617       79.72        73       75.26       544       80.35
# TOTAL (rows)       774      100.00        97      100.00       677      100.00
# TOTAL (rows)       774      100.00        97      100.00       677      100.00
#

write_xlsx(flag_summary, file.path("ch", "data", "ch_flag_summary.xlsx"))

# ========================
# UPSET
# ========================
all_ch <- read_excel(file.path("ch", "data", "ch_wl_art_flag_vars.xlsx"))

overlap_flags <- c(
    "flag_freeze_enriched", "flag_not_age_associated",
    "flag_germline_ind", "flag_germline_var",
    "flag_cluster", "flag_gnomad"
)

flag_mat <- all_ch2 %>%
    select(all_of(overlap_flags)) %>%
    mutate(across(everything(), as.integer))

### UPSET PLOT
pdf(file.path("ch", "figures", "qc_flag_overlap_upset.pdf"), width = 10, height = 5)
upset(
    as.data.frame(flag_mat),
    sets           = overlap_flags,
    order.by       = "freq",
    sets.bar.color = "#378ADD",
    main.bar.color = "#1D9E75",
    text.scale     = 1.2,
    mb.ratio       = c(0.6, 0.4)
)
dev.off()

all_ch2 <- all_ch2 %>%
    mutate(n_flags = rowSums(across(all_of(overlap_flags), as.integer), na.rm = TRUE))

cat("\nRows by number of flags triggered:\n")
print(table(all_ch2$n_flags))
# 0   1   2   3
# 359 266 141   8

write_xlsx(all_ch2 %>% filter(n_flags == 1), file.path("ch", "data", "ch_flag1_examine.xlsx"))
write_xlsx(all_ch2 %>% filter(n_flags >= 2), file.path("ch", "data", "ch_flag2_examine.xlsx"))

write_xlsx(all_ch3, file.path("ch", "data", "ch_wl_art_flag_vars.xlsx"))




# ========================
# OVERLAP VENN DIAGRAMS
# ========================
var_f2 <- sort(variant_counts %>% filter(batch_1 != 0) %>% pull(variant_id))
var_f3 <- sort(variant_counts %>% filter(batch_2 != 0) %>% pull(variant_id))

p_venn <- ggVennDiagram(
    list("Freeze 2" = var_f2, "Freeze 3" = var_f3),
    label = "count", label_alpha = 0
) +
    scale_fill_gradient(low = "white", high = "#1D9E75") +
    scale_color_manual(values = c("grey40", "grey40")) +
    labs(title = "Variant overlap between freezes") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
          legend.position = "none")

ggsave(file.path("ch", "figures", "qc_variant_overlap_venn.pdf"), p_venn, width = 6, height = 4)

gene_f2 <- sort(unique(all_ch %>% filter(Batch == 1) %>% pull(Gene)))
gene_f3 <- sort(unique(all_ch %>% filter(Batch == 2) %>% pull(Gene)))

p_venn_gene <- ggVennDiagram(
    list("Freeze 2" = gene_f2, "Freeze 3" = gene_f3),
    label = "count", label_alpha = 0
) +
    scale_fill_gradient(low = "white", high = "#378ADD") +
    scale_color_manual(values = c("grey40", "grey40")) +
    labs(title = "Gene overlap between freezes") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
          legend.position = "none")

ggsave(file.path("ch", "figures", "qc_gene_overlap_venn.pdf"), p_venn_gene, width = 6, height = 4)


