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
# READ IN DATA
# ========================
all_ch <- read_excel(file.path("ch", "data", "ch_seq_wl_vars.xlsx"))
dim(all_ch)
# 774

cov <- read.csv(file.path("ch", "data", "pmbb_brca12_cov_df.csv"), row.names = 1)
all_ids <- read.csv(file.path("ch", "data", "ch_psm_matched4_case_control_ids.csv"))$x
cov <- cov %>% filter(person_id %in% all_ids)
table(cov$Batch)
dim(cov)
# 3004

# ========================
# SET GLOBALS
# ========================
GNOMAD_THRESHOLD  <- 0.001
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
write_xlsx(variant_counts, file.path("ch", "data", "ch_wl_var_counts.xlsx"))

high_freq      <- variant_counts %>% filter(flag_high_freq == 1)
high_freq_vars <- unique(high_freq$variant_id)
cat(sprintf("High frequency variants to test: %d\n", nrow(high_freq)))
# 18

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
        p_freeze_adj         = p.adjust(p_freeze, method = "BH"),
        flag_freeze_enriched = p_freeze_adj < FREEZE_SIG
    )

freeze_enriched_vars <- high_freq2 %>% filter(flag_freeze_enriched) %>% pull(variant_id)
cat(sprintf("Freeze-enriched variants: %d\n", length(freeze_enriched_vars)))
# 2

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

# ========================
# 3. FLAG - AGE ASSOCIATION
# ========================
artifact_test_age <- function(vkey, cov, extra_covs = NULL) {
    carriers <- all_ch %>%
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
        separation <- any(fitted(m) %in% c(0, 1))
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

high_freq_base <- high_freq2 %>%
    rowwise() %>%
    mutate(res = list(artifact_test_age(variant_id, cov))) %>%
    unnest(res) %>% ungroup() %>%
    mutate(
        p_age_adj           = p.adjust(p_age, method = "BH"),
        flag_age_associated = !is.na(p_age) & !separation & p_age < AGE_SIG
    )

high_freq_full <- high_freq2 %>%
    rowwise() %>%
    mutate(res = list(artifact_test_age(variant_id, cov, extra_covs = extra_covs))) %>%
    unnest(res) %>% ungroup() %>%
    mutate(
        p_age_adj           = p.adjust(p_age, method = "BH"),
        flag_age_associated = !is.na(p_age) & !separation & p_age < AGE_SIG
    )

not_age_associated_var <- high_freq_full %>% filter(!flag_age_associated) %>% pull(variant_id)
write_xlsx(high_freq_full, file.path("ch", "data", "ch_high_freq_counts.xlsx"))

# ========================
# 4. FLAG - GERMLINE TESTING
# ========================
binomial_germline_test <- function(alt_count, total_depth, p_germline = 0.5) {
    if (is.na(alt_count) || is.na(total_depth) || total_depth == 0) {
        return(data.frame(p_binom = NA_real_, flag_germline_ind = NA))
    }
    p <- binom.test(x = alt_count, n = total_depth,
                    p = p_germline, alternative = "two.sided")$p.value
    data.frame(
        p_binom           = p,
        flag_germline_ind = p >= GERMLINE_SIG
    )
}

germline_test_vars <- all_ch %>%
    filter(
        (Gene %in% c("TET2", "CBL") & variant_category == "missense") |
            variant_id %in% high_freq2$variant_id
    ) %>%
    distinct(variant_id, Sample.ID, Sample.AltDepth, Sample.Depth) %>%
    rowwise() %>%
    mutate(res = list(binomial_germline_test(Sample.AltDepth, Sample.Depth))) %>%
    unnest(res) %>%
    ungroup()

germline_summary <- germline_test_vars %>%
    group_by(variant_id) %>%
    summarise(
        n_tested        = sum(!is.na(p_binom)),
        n_fail_binom    = sum(flag_germline_ind, na.rm = TRUE),
        prop_fail_binom = n_fail_binom / n_tested,
        flag_germline_var = prop_fail_binom > 0.5,
        .groups         = "drop"
    )

write_xlsx(germline_test_vars, file.path("ch", "data", "ch_germline_ind_df.xlsx"))
write_xlsx(germline_summary, file.path("ch", "data", "ch_germline_vars_df.xlsx"))

germline_ind_pairs <- germline_test_vars %>%
    filter(flag_germline_ind) %>%
    select(Sample.ID, variant_id)
germline_var <- germline_summary %>% filter(flag_germline_var) %>% pull(variant_id)

cat(sprintf("Variants tested for germline: %d\n", nrow(germline_summary)))
# 41
cat(sprintf("Individuals flagged germline: %d\n",  nrow(germline_ind_pairs)))
# 18
cat(sprintf("Variants flagged germline:    %d\n",  sum(germline_summary$flag_germline_var, na.rm = TRUE)))
# 5

# ========================
# 5. FLAG - CLUSTERS
# ========================
all_ch_clustered <- all_ch %>%
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

# ========================
# JOIN BACK TO ALL_CH & VAF STRATA
# ========================
all_ch2 <- all_ch %>%
    mutate(
        flag_freeze_enriched    = variant_id %in% freeze_enriched_vars,
        flag_not_age_associated = variant_id %in% not_age_associated_var,
        flag_germline_ind       = paste(Sample.ID, variant_id) %in% paste(germline_ind_pairs$Sample.ID, germline_ind_pairs$variant_id),
        flag_germline_var       = variant_id %in% germline_var,
        flag_cluster            = variant_id %in% cluster_not_rep_var,
        flag_gnomad             = !is.na(gnomAD.MAX_AF) & gnomAD.MAX_AF > GNOMAD_THRESHOLD,
        VAF_Strata = case_when(
            Sample.AltFrac >= 0.02 & Sample.AltFrac < 0.10 ~ "2-10",
            Sample.AltFrac >= 0.10                          ~ "10+",
            TRUE                                            ~ NA_character_
        )
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

# ========================
# DEFINE FILTER CONFIGURATIONS
# ========================
full_qc <- c("flag_freeze_enriched", "flag_not_age_associated",
             "flag_germline_ind", "flag_germline_var", "flag_cluster", "flag_gnomad")

filter_configs <- list(
    list(name = "no_filter",               flags = character(0)),
    list(name = "covs",                    flags = character(0), formula_extra = extra_covs),
    list(name = "flag_freeze_enriched",    flags = c("flag_freeze_enriched")),
    list(name = "flag_not_age_associated", flags = c("flag_not_age_associated")),
    list(name = "flag_germline_ind",       flags = c("flag_germline_ind")),
    list(name = "flag_germline_var",       flags = c("flag_germline_var")),
    list(name = "flag_cluster",            flags = c("flag_cluster")),
    list(name = "flag_gnomad",             flags = c("flag_gnomad")),
    list(name = "full_qc_no_minad",        flags = full_qc)
)

minad_configs <- list(
    list(name = "minad3",      flags = c("minad3")),
    list(name = "minad4",      flags = c("minad4")),
    list(name = "minad5",      flags = c("minad5")),
    list(name = "minad6",      flags = c("minad6")),
    list(name = "minad7",      flags = c("minad7")),
    list(name = "minad8",      flags = c("minad8")),
    list(name = "minad3_full", flags = c("minad3", full_qc)),
    list(name = "minad4_full", flags = c("minad4", full_qc)),
    list(name = "minad5_full", flags = c("minad5", full_qc)),
    list(name = "minad6_full", flags = c("minad6", full_qc)),
    list(name = "minad7_full", flags = c("minad7", full_qc)),
    list(name = "minad8_full", flags = c("minad8", full_qc))
)

cov_configs <- list(
    list(name = "full_covs",       flags = full_qc,                formula_extra = extra_covs),
    list(name = "minad3_full_covs",  flags = c("minad3",  full_qc), formula_extra = extra_covs),
    list(name = "minad4_full_covs",  flags = c("minad4",  full_qc), formula_extra = extra_covs),
    list(name = "minad5_full_covs",  flags = c("minad5",  full_qc), formula_extra = extra_covs),
    list(name = "minad6_full_covs",  flags = c("minad6",  full_qc), formula_extra = extra_covs),
    list(name = "minad7_full_covs",  flags = c("minad7",  full_qc), formula_extra = extra_covs),
    list(name = "minad8_full_covs",  flags = c("minad8",  full_qc), formula_extra = extra_covs)
)

# ========================
# AGE ASSOCIATION TESTER
# ========================
test_age_association <- function(config, all_ch, cov,
                                 default_formula_extra = NULL,
                                 stratify_by = NULL) {
    if (!is.null(stratify_by)) {
        strata_levels <- sort(unique(na.omit(all_ch[[stratify_by]])))
        strata_list   <- setNames(as.list(strata_levels), strata_levels)
    } else {
        strata_list <- list("all" = NULL)
    }

    results <- lapply(config, function(cfg) {
        nm        <- cfg$name
        flags     <- cfg$flags
        fml_extra <- if (!is.null(cfg$formula_extra)) cfg$formula_extra else default_formula_extra

        lapply(names(strata_list), function(stratum) {
            stratum_val <- strata_list[[stratum]]

            ch_filtered <- all_ch
            if (!is.null(stratum_val)) {
                ch_filtered <- ch_filtered %>% filter(.data[[stratify_by]] == stratum_val)
            }
            if (length(flags) > 0) {
                remove_mask <- Reduce(`|`, lapply(flags, function(f) {
                    col <- ch_filtered[[f]]
                    if (is.logical(col)) col else as.logical(col)
                }))
                remove_mask[is.na(remove_mask)] <- FALSE
                ch_filtered <- ch_filtered[!remove_mask, ]
            }

            carriers <- ch_filtered %>%
                distinct(Sample.ID) %>%
                mutate(carrier = 1L)

            df <- cov %>%
                left_join(carriers, by = c("person_id" = "Sample.ID")) %>%
                mutate(carrier = as.integer(ifelse(is.na(carrier), 0L, carrier)))

            n_carriers <- sum(df$carrier)
            n_total    <- nrow(df)

            cat(sprintf("[%s | stratum: %s] %d carriers / %d total\n",
                        nm, stratum, n_carriers, n_total))

            base_rhs <- if (!is.null(stratify_by) && stratify_by == "Batch") {
                "Sample_age"
            } else {
                "Sample_age + as.factor(Batch)"
            }
            rhs     <- if (!is.null(fml_extra)) paste(base_rhs, "+", fml_extra) else base_rhs
            formula <- as.formula(paste("carrier ~", rhs))

            tryCatch({
                m <- withCallingHandlers(
                    glm(formula, data = df, family = binomial()),
                    warning = function(w) {
                        if (grepl("fitted probabilities numerically 0 or 1", conditionMessage(w)))
                            invokeRestart("muffleWarning")
                    }
                )
                separation <- any(fitted(m) %in% c(0, 1))
                s  <- coef(summary(m))
                ci <- confint.default(m)
                data.frame(
                    name       = nm, stratum = stratum,
                    flags      = paste(flags, collapse = ", "),
                    n_carriers = n_carriers, n_total = n_total,
                    OR         = exp(coef(m)["Sample_age"]),
                    CI_lo      = exp(ci["Sample_age", 1]),
                    CI_hi      = exp(ci["Sample_age", 2]),
                    p_age      = s["Sample_age", "Pr(>|z|)"],
                    separation = separation,
                    stringsAsFactors = FALSE
                )
            }, error = function(e) {
                warning(sprintf("Model failed for [%s | %s]: %s", nm, stratum, conditionMessage(e)))
                data.frame(name = nm, stratum = stratum,
                           flags = paste(flags, collapse = ", "),
                           n_carriers = n_carriers, n_total = n_total,
                           OR = NA_real_, CI_lo = NA_real_, CI_hi = NA_real_,
                           p_age = NA_real_, separation = NA,
                           stringsAsFactors = FALSE)
            })
        }) %>% bind_rows()
    }) %>% bind_rows()

    results
}

# ========================
# PLOT FUNCTION
# ========================
plot_age_association <- function(results, title = "Age association under different QC filters",
                                 n_total = NULL) {
    stratified <- length(unique(results$stratum)) > 1

    plot_df <- results %>%
        filter(!is.na(OR)) %>%
        mutate(
            name    = factor(name, levels = rev(unique(name))),
            sig     = !is.na(p_age) & p_age < 0.05,
            stratum = factor(stratum)
        )

    subtitle <- if (!is.null(n_total)) sprintf("n = %d total individuals", n_total) else NULL

    p <- ggplot(plot_df, aes(x = OR, y = name, color = if (stratified) stratum else sig)) +
        geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
        geom_errorbar(aes(xmin = CI_lo, xmax = CI_hi, group = if (stratified) stratum else NULL),
                      width = 0.25, linewidth = 0.6, orientation = "y",
                      position = if (stratified) position_dodge(width = 0.6) else "identity") +
        geom_point(size = 3,
                   position = if (stratified) position_dodge(width = 0.6) else "identity") +
        scale_x_continuous(name = "Odds ratio per year of age (95% CI)",
                           breaks = scales::pretty_breaks(n = 6)) +
        labs(title = title, subtitle = subtitle, y = NULL) +
        theme_minimal(base_size = 11) +
        theme(
            plot.title         = element_text(face = "bold", hjust = 0.5),
            plot.subtitle      = element_text(hjust = 0.5, size = 9, color = "grey40"),
            panel.grid.minor   = element_blank(),
            panel.grid.major.y = element_blank(),
            legend.position    = "bottom"
        )

    if (stratified) {
        p <- p +
            scale_color_brewer(palette = "Set2") +
            guides(color = guide_legend(title = "stratum"))
    } else {
        p <- p +
            scale_color_manual(
                values = c("FALSE" = "grey50", "TRUE" = "#D85A30"),
                labels = c("FALSE" = "p \u2265 0.05", "TRUE" = "p < 0.05"),
                name   = NULL
            )
    }
    p
}

# ========================
# RUN
# ========================
filter_res       <- test_age_association(filter_configs, all_ch2, cov)
minad_res        <- test_age_association(minad_configs,  all_ch2, cov)
cov_res          <- test_age_association(cov_configs,    all_ch2, cov)

filter_vaf_res   <- test_age_association(filter_configs, all_ch2, cov, stratify_by = "VAF_Strata")
minad_vaf_res    <- test_age_association(minad_configs,  all_ch2, cov, stratify_by = "VAF_Strata")
cov_vaf_res      <- test_age_association(cov_configs,    all_ch2, cov, stratify_by = "VAF_Strata")

filter_freeze_res <- test_age_association(filter_configs, all_ch2, cov, stratify_by = "Batch")
minad_freeze_res  <- test_age_association(minad_configs,  all_ch2, cov, stratify_by = "Batch")
cov_freeze_res    <- test_age_association(cov_configs,    all_ch2, cov, stratify_by = "Batch")

# ========================
# SAVE RESULTS
# ========================
write_xlsx(filter_res,        file.path("ch", "data", "ch_filter_age_sensitivity.xlsx"))
write_xlsx(minad_res,         file.path("ch", "data", "ch_minad_age_sensitivity.xlsx"))
write_xlsx(cov_res,           file.path("ch", "data", "ch_cov_age_sensitivity.xlsx"))
write_xlsx(filter_vaf_res,    file.path("ch", "data", "ch_filter_age_sensitivity_vaf.xlsx"))
write_xlsx(minad_vaf_res,     file.path("ch", "data", "ch_minad_age_sensitivity_vaf.xlsx"))
write_xlsx(cov_vaf_res,       file.path("ch", "data", "ch_cov_age_sensitivity_vaf.xlsx"))
write_xlsx(filter_freeze_res, file.path("ch", "data", "ch_filter_age_sensitivity_batch.xlsx"))
write_xlsx(minad_freeze_res,  file.path("ch", "data", "ch_minad_age_sensitivity_batch.xlsx"))
write_xlsx(cov_freeze_res,    file.path("ch", "data", "ch_cov_age_sensitivity_batch.xlsx"))

# ========================
# PLOTS
# ========================
save_assoc_plot <- function(results, filename, title, n_total) {
    p <- plot_age_association(results, title = title, n_total = n_total)
    n_names   <- n_distinct(results$name)
    n_strata  <- n_distinct(results$stratum)
    width  <- if (n_strata > 1) 5 + 1.5 * n_strata else 7
    height <- 0.4 * n_names + 2
    ggsave(filename, p, width = width, height = height, limitsize = FALSE)
}

save_assoc_plot(filter_res,        file.path("ch", "figures", "qc_filter_age_sensitivity.pdf"),
                "Age association — QC filters", nrow(cov))
save_assoc_plot(minad_res,         file.path("ch", "figures", "qc_minad_age_sensitivity.pdf"),
                "Age association — minAD thresholds", nrow(cov))
save_assoc_plot(cov_res,           file.path("ch", "figures", "qc_cov_age_sensitivity.pdf"),
                "Age association — covariate models", nrow(cov))

save_assoc_plot(filter_vaf_res,    file.path("ch", "figures", "qc_filter_age_sensitivity_vaf.pdf"),
                "Age association — QC filters by VAF stratum", nrow(cov))
save_assoc_plot(minad_vaf_res,     file.path("ch", "figures", "qc_minad_age_sensitivity_vaf.pdf"),
                "Age association — minAD by VAF stratum", nrow(cov))
save_assoc_plot(cov_vaf_res,       file.path("ch", "figures", "qc_cov_age_sensitivity_vaf.pdf"),
                "Age association — covariate models by VAF stratum", nrow(cov))

save_assoc_plot(filter_freeze_res, file.path("ch", "figures", "qc_filter_age_sensitivity_batch.pdf"),
                "Age association — QC filters by batch", nrow(cov))
save_assoc_plot(minad_freeze_res,  file.path("ch", "figures", "qc_minad_age_sensitivity_batch.pdf"),
                "Age association — minAD by batch", nrow(cov))
save_assoc_plot(cov_freeze_res,    file.path("ch", "figures", "qc_cov_age_sensitivity_batch.pdf"),
                "Age association — covariate models by batch", nrow(cov))

# ========================
# COMPARE VAF 10< vs. ALL VAF
# ========================
vaf_compare_res <- bind_rows(
    test_age_association(cov_configs, all_ch2, cov) %>% mutate(vaf_group = "All VAF"),
    test_age_association(cov_configs, all_ch2 %>% filter(Sample.AltFrac >= 0.10), cov) %>%
        mutate(vaf_group = "VAF ≥ 10%")
)

plot_df <- vaf_compare_res %>%
    filter(!is.na(OR)) %>%
    mutate(
        name      = factor(name, levels = rev(unique(name))),
        sig       = !is.na(p_age) & p_age < 0.05,
        vaf_group = factor(vaf_group, levels = c("All VAF", "VAF ≥ 10%"))
    )

p_vaf_compare <- ggplot(plot_df, aes(x = OR, y = name, color = vaf_group)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    geom_errorbar(aes(xmin = CI_lo, xmax = CI_hi, group = vaf_group),
                  width = 0.25, linewidth = 0.6, orientation = "y",
                  position = position_dodge(width = 0.6)) +
    geom_point(aes(shape = sig), size = 3, position = position_dodge(width = 0.6)) +
    scale_color_manual(
        values = c("All VAF" = "#378ADD", "VAF ≥ 10%" = "#D85A30"),
        name   = NULL
    ) +
    scale_shape_manual(
        values = c("FALSE" = 16, "TRUE" = 8),
        labels = c("FALSE" = "p \u2265 0.05", "TRUE" = "p < 0.05"),
        name   = NULL
    ) +
    scale_x_continuous(
        name   = "Odds ratio per year of age (95% CI)",
        breaks = scales::pretty_breaks(n = 6)
    ) +
    labs(
        title    = "Age association by VAF threshold",
        subtitle = sprintf("n = %d total individuals", nrow(cov)),
        y        = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
        plot.title         = element_text(face = "bold", hjust = 0.5),
        plot.subtitle      = element_text(hjust = 0.5, size = 9, color = "grey40"),
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position    = "bottom"
    )

ggsave(file.path("ch", "figures", "qc_filter_age_cov_compare.pdf"),
       p_vaf_compare,
       width  = 7,
       height = 0.4 * n_distinct(plot_df$name) + 2,
       limitsize = FALSE)


# ========================
# FINAL FILTER
# ========================
full_qc <- c("flag_freeze_enriched", "flag_not_age_associated",
             "flag_germline_ind", "flag_germline_var", "flag_cluster", "flag_gnomad")
# minad5

all_ch3 <- all_ch2 %>% filter(n_flags == 0, minad5)
dim(all_ch3)
# ========================
# TOP GENES BY FREEZE (POST-QC)
# ========================
top_genes <- all_ch3 %>%
    count(Gene) %>%
    slice_max(n, n = 20) %>%
    pull(Gene)

gene_counts <- all_ch3 %>%
    filter(Gene %in% top_genes) %>%
    mutate(
        Freeze = factor(Batch, levels = c(1, 2), labels = c("Freeze 2", "Freeze 3")),
        Gene   = factor(Gene, levels = top_genes)  # ordered by overall count
    ) %>%
    count(Gene, Freeze)

p_genes <- ggplot(gene_counts, aes(x = Gene, y = n, fill = Freeze)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_fill_manual(values = c("Freeze 2" = "#378ADD", "Freeze 3" = "#1D9E75")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
        title    = "Top genes after QC filtering",
        subtitle = sprintf("minAD ≥ 5 + full QC | n = %d variants", nrow(all_ch3)),
        x        = NULL, y = "Variant count", fill = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
        plot.title       = element_text(face = "bold", hjust = 0.5),
        plot.subtitle    = element_text(hjust = 0.5, size = 9, color = "grey40"),
        axis.text.x      = element_text(angle = 40, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        legend.position    = "top"
    )

ggsave(file.path("ch", "figures", "ch_top_genes_by_freeze.pdf"),
       p_genes, width = 9, height = 5)
write_xlsx(all_ch3, file.path("ch", "data", "ch3.xlsx"))

