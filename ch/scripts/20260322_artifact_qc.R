# ========================
# Age association for artifacts
# ========================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
library(data.table, quietly = T)
library(ggVennDiagram)

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
# 1    2
# 760 2244
dim(cov)
# 3004

# ========================
# SET GLOBALS
# ========================
GNOMAD_THRESHOLD <- 0.001
FREQ_THRESHOLD <- 4
CLUSTER_THRESHOLD <- 50
MIN_AD_THRESHOLDS <- 3:10

AGE_SIG <- 0.10
FREEZE_SIG <- 0.05
GERMLINE_SIG <- 0.05

n_f2 <- 760
n_f3 <- 2244

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
# 481 unique variants

write_xlsx(variant_counts, file.path("ch", "data", "ch_wl_var_counts.xlsx"))

high_freq <- variant_counts %>% filter(flag_high_freq == 1)
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
        p_freeze_adj = p.adjust(p_freeze, method = "BH"),
        flag_freeze_enriched = p_freeze_adj < FREEZE_SIG
    )

freeze_enriched_vars <- high_freq2 %>% filter(flag_freeze_enriched) %>% pull(variant_id)

# ========================
# OVERLAP VENN DIAGRAMS
# ========================
var_f2 <- sort(variant_counts %>% filter(batch_1 != 0) %>% pull(variant_id))
var_f3 <- sort(variant_counts %>% filter(batch_2 != 0) %>% pull(variant_id))

p_venn <- ggVennDiagram(
    list("Freeze 2" = var_f2, "Freeze 3" = var_f3),
    label       = "count",
    label_alpha = 0
) +
    scale_fill_gradient(low = "white", high = "#1D9E75") +
    scale_color_manual(values = c("grey40", "grey40")) +
    labs(title = "Variant overlap between freezes") +
    theme(
        plot.title  = element_text(face = "bold", hjust = 0.5, size = 13),
        legend.position = "none"
    )

ggsave(file.path("ch", "figures", "qc_variant_overlap_venn.pdf"),
       p_venn, width = 6, height = 4)

### OVERLAP BTWN FREEZES - GENES ###
gene_f2 <- sort(unique(all_ch %>% filter(Batch == 1) %>% pull(Gene)))
gene_f3 <- sort(unique(all_ch %>% filter(Batch == 2) %>% pull(Gene)))

p_venn_gene <- ggVennDiagram(
    list("Freeze 2" = gene_f2, "Freeze 3" = gene_f3),
    label       = "count",
    label_alpha = 0
) +
    scale_fill_gradient(low = "white", high = "#378ADD") +
    scale_color_manual(values = c("grey40", "grey40")) +
    labs(title = "Gene overlap between freezes") +
    theme(
        plot.title      = element_text(face = "bold", hjust = 0.5, size = 13),
        legend.position = "none"
    )

ggsave(file.path("ch", "figures", "qc_gene_overlap_venn.pdf"),
       p_venn_gene, width = 6, height = 4)

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

    # build formula dynamically
    base <- "carrier ~ Sample_age + as.factor(Batch)"
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

# base model
high_freq_base <- high_freq2 %>%
    rowwise() %>%
    mutate(res = list(artifact_test_age(variant_id, cov))) %>%
    unnest(res) %>% ungroup() %>%
    mutate(
        p_age_adj = p.adjust(p_age, method = "BH"),
        flag_age_associated = !is.na(p_age) & !separation & p_age < AGE_SIG
    )

extra_covs <- "Smoke_History + Sequenced_gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6"

high_freq_full <- high_freq2 %>%
    rowwise() %>%
    mutate(res = list(artifact_test_age(variant_id, cov, extra_covs = extra_covs))) %>%
    unnest(res) %>% ungroup() %>%
    mutate(
        p_age_adj = p.adjust(p_age, method = "BH"),
        flag_age_associated = !is.na(p_age) & !separation & p_age < AGE_SIG
    )

not_age_associated_var <- high_freq_full %>% filter(!flag_age_associated) %>% pull(variant_id)

write_xlsx(high_freq_full, file.path("ch", "data", "ch_high_freq_counts.xlsx"))

# ========================
# 4. FLAG - GERMLINE TESTING
# ========================
# Tests whether observed VAF is consistent with germline het (expected = 0.5)
# Flags variants where we cannot reject germline origin at P < 0.01
binomial_germline_test <- function(alt_count, total_depth, p_germline = 0.5) {
    if (is.na(alt_count) || is.na(total_depth) || total_depth == 0) {
        return(data.frame(p_binom = NA_real_, flag_germline_binom = NA))
    }
    # two-sided: is VAF significantly different from 0.5?
    p <- binom.test(x = alt_count, n = total_depth,
                    p = p_germline, alternative = "two.sided")$p.value
    data.frame(
        p_binom            = p,
        flag_germline_ind = p >= GERMLINE_SIG
    )
}

# apply to TET2/CBL missense + high_freq variants
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

# summarize per variant: flag if MAJORITY of carriers fail binomial test
germline_summary <- germline_test_vars %>%
    group_by(variant_id) %>%
    summarise(
        n_tested          = sum(!is.na(p_binom)),
        n_fail_binom      = sum(flag_germline_ind, na.rm = TRUE),
        prop_fail_binom   = n_fail_binom / n_tested,
        flag_germline_var     = prop_fail_binom > 0.5,   # majority of carriers look germline
        .groups           = "drop"
    )


germline_ind <- germline_test_vars %>% filter(flag_germline_ind) %>% pull(variant_id)
germline_var <- germline_summary %>% filter(flag_germline_var) %>% pull(variant_id)

cat(sprintf("Variants tested for germline: %d\n",    nrow(germline_summary)))
# 41
cat(sprintf("Individuals flagged germline: %d\n",    sum(germline_test_vars$flag_germline_ind, na.rm = TRUE)))
# 18
cat(sprintf("Variants flagged germline: %d\n",    sum(germline_summary$flag_germline_var, na.rm = TRUE)))
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
        # new cluster starts when there's a gap > threshold to the previous row
        new_cluster = !near_prev,
        cluster_id  = paste(Sample.ID, Gene, Chr,
                            cumsum(new_cluster), sep = "_")
    ) %>%
    ungroup() %>%
    select(-new_cluster, -.row_idx)

cluster_rep <- cluster_all %>%
    group_by(cluster_id) %>%
    arrange(
        !grepl("\\bC[0-9]", Existing.variation),  # COSMIC rows sort first (FALSE < TRUE)
        Start
    ) %>%
    slice(1) %>%
    ungroup()

cluster_not_rep_var <- cluster_all %>% filter(!(variant_id %in% cluster_rep$variant_id)) %>%
    pull(variant_id)

# ========================
# JOIN BACK TO ALL_CH & VAF STRATA
# ========================
all_ch2 <- all_ch %>%
    mutate(
        flag_not_high_freq = ifelse(variant_id %in% high_freq_vars, FALSE, TRUE),
        flag_not_freeze_enriched = ifelse(variant_id %in% freeze_enriched_vars, FALSE, TRUE),
        flag_not_age_associated = ifelse(variant_id %in% not_age_associated_var, FALSE, TRUE),
        flag_germline_ind = ifelse(variant_id %in% germline_ind, FALSE, TRUE),
        flag_germline_var = ifelse(variant_id %in% germline_var, FALSE, TRUE),
        flag_cluster = ifelse(variant_id %in% cluster_not_rep_var, FALSE, TRUE),
    ) %>% mutate(
        flag_gnomad = ifelse(gnomAD.MAX_AF > GNOMAD_THRESHOLD, TRUE, FALSE)
    ) %>% mutate(
        VAF_Strata = case_when(
            Sample.AltFrac >= 0.02 & Sample.AltFrac < 0.10 ~ "2-10",
            Sample.AltFrac >= 0.10 ~ "10+",
            TRUE ~ NA_character_
        )
    )

### DEFINE MINAD THRESHOLDS
for (thresh in MIN_AD_THRESHOLDS) {
    col_name <- paste0("minad", thresh)
    all_ch2[[col_name]] <- all_ch2$Sample.AltDepth < thresh
}

write_xlsx(all_ch2, file.path("ch", "data", "ch_seq_wl_art_flags_vars.xlsx"))

# ========================
# DEFINE FILTER CONFIGURATIONS
# ========================
full_qc <- c("flag_not_high_freq", "flag_not_freeze_enriched", "flag_not_age_associated",
             "flag_germline_ind", "flag_germline_var", "flag_cluster", "flag_gnomad")

filter_configs <- list(
    # --- baseline: no filters ---
    list(name = "no_filter", flags = character(0)),
    list(name = "covs", flags = character(0), formula_extra = extra_covs),

    # --- single flags ---
    list(name = "flag_not_high_freq", flags = c("flag_not_high_freq")),
    list(name = "flag_not_freeze_enriched", flags = c("flag_not_freeze_enriched")),
    list(name = "flag_not_age_associated", flags = c("flag_not_age_associated")),
    list(name = "flag_germline_ind", flags = c("flag_germline_ind")),
    list(name = "flag_germline_var", flags = c("flag_germline_var")),
    list(name = "flag_cluster", flags = c("flag_cluster")),
    list(name = "flag_gnomad", flags = c("flag_gnomad")),

    # --- combined QC flags (no minAD) ---
    list(name = "full_qc_no_minad", flags = full_qc)
)

minad_configs <- list(
    # --- minAD thresholds alone ---
    list(name = "minad3",               flags = c("minad3")),
    list(name = "minad4",               flags = c("minad4")),
    list(name = "minad5",               flags = c("minad5")),
    list(name = "minad6",               flags = c("minad6")),
    list(name = "minad7",               flags = c("minad7")),

    # --- minAD + full QC combinations ---
    list(name = "minad3_full", flags = c("minad3", full_qc)),
    list(name = "minad4_full", flags = c("minad4", full_qc)),
    list(name = "minad5_full", flags = c("minad5", full_qc)),
    list(name = "minad6_full", flags = c("minad6", full_qc)),
    list(name = "minad7_full", flags = c("minad7", full_qc)),
    list(name = "minad8_full", flags = c("minad8", full_qc)),
    list(name = "minad9_full", flags = c("minad9", full_qc)),
    list(name = "minad9_full", flags = c("minad10", full_qc))
)

cov_configs <- list(
    # --- full QC + extra covariates ---
    list(name = "full_covs",
         flags = c(full_qc),
         formula_extra = extra_covs),
    list(name = "minad3_full_covs",
         flags = c("minad3", full_qc),
         formula_extra = extra_covs),
    list(name = "minad4_full_covs",
         flags = c("minad4", full_qc),
         formula_extra = extra_covs),
    list(name = "minad5_full_covs",
         flags = c("minad5", full_qc),
         formula_extra = extra_covs),
    list(name = "minad6_full_covs",
         flags = c("minad6", full_qc),
         formula_extra = extra_covs),
    list(name = "minad7_full_covs",
         flags = c("minad7", full_qc),
         formula_extra = extra_covs),
    list(name = "minad8_full_covs",
         flags = c("minad8", full_qc),
         formula_extra = extra_covs),
    list(name = "minad9_full_covs",
         flags = c("minad7", full_qc),
         formula_extra = extra_covs),
    list(name = "minad10_full_covs",
         flags = c("minad10", full_qc),
         formula_extra = extra_covs)
)

# ========================
# AGE ASSOCIATION TESTER
# ========================
# config: named list where each element has:
#   $name          : label for the row in results / plot
#   $flags         : character vector of column names in all_ch5; rows where ANY flag == TRUE are excluded
#   $formula_extra : (optional) string of extra RHS covariates
#
# stratify_by: (optional) column name in all_ch5 to stratify on (e.g. "VAF_Strata", "Batch")
#              runs the model separately for each level; results get a `stratum` column

test_age_association <- function(config, all_ch, cov,
                                 default_formula_extra = NULL,
                                 stratify_by = NULL) {

    # determine strata levels up front
    if (!is.null(stratify_by)) {
        strata_levels <- sort(unique(na.omit(all_ch[[stratify_by]])))
        strata_list   <- setNames(as.list(strata_levels), strata_levels)
    } else {
        strata_list <- list("all" = NULL)   # single pass, no subsetting
    }

    results <- lapply(config, function(cfg) {
        nm        <- cfg$name
        flags     <- cfg$flags
        fml_extra <- if (!is.null(cfg$formula_extra)) cfg$formula_extra else default_formula_extra

        lapply(names(strata_list), function(stratum) {
            stratum_val <- strata_list[[stratum]]

            # --- filter variants by QC flags ---
            ch_filtered <- all_ch
            if (!is.null(stratum_val)) {
                ch_filtered <- ch_filtered %>%
                    filter(.data[[stratify_by]] == stratum_val)
            }
            if (length(flags) > 0) {
                remove_mask <- Reduce(`|`, lapply(flags, function(f) {
                    col <- ch_filtered[[f]]
                    if (is.logical(col)) col else as.logical(col)
                }))
                remove_mask[is.na(remove_mask)] <- FALSE
                ch_filtered <- ch_filtered[!remove_mask, ]
            }

            # --- carrier status ---
            carriers <- ch_filtered %>%
                distinct(Sample.ID) %>%
                mutate(carrier = 1L)

            df <- cov %>%
                left_join(carriers, by = c("person_id" = "Sample.ID")) %>%
                mutate(carrier = as.integer(ifelse(is.na(carrier), 0L, carrier)))

            n_carriers <- sum(df$carrier)
            n_total    <- nrow(df)

            cat(sprintf("[%s | stratum: %s] %d carriers / %d total | flags: %s\n",
                        nm, stratum, n_carriers, n_total,
                        if (length(flags) == 0) "none" else paste(flags, collapse = ", ")))

            # --- build formula ---
            # if stratifying by Batch, drop as.factor(Batch) from base to avoid collinearity
            base_rhs <- if (!is.null(stratify_by) && stratify_by == "Batch") {
                "Sample_age"
            } else {
                "Sample_age + as.factor(Batch)"
            }
            rhs     <- if (!is.null(fml_extra)) paste(base_rhs, "+", fml_extra) else base_rhs
            formula <- as.formula(paste("carrier ~", rhs))

            # --- fit model ---
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
                    name        = nm,
                    stratum     = stratum,
                    flags       = paste(flags, collapse = ", "),
                    n_carriers  = n_carriers,
                    n_total     = n_total,
                    OR          = exp(coef(m)["Sample_age"]),
                    CI_lo       = exp(ci["Sample_age", 1]),
                    CI_hi       = exp(ci["Sample_age", 2]),
                    p_age       = s["Sample_age", "Pr(>|z|)"],
                    separation  = separation,
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
# PLOT FUNCTION (works for stratified and unstratified)
# ========================
plot_age_association <- function(results, title = "Age association under different QC filters",
                                 n_total = NULL) {
    stratified <- length(unique(results$stratum)) > 1

    plot_df <- results %>%
        filter(!is.na(OR)) %>%
        mutate(
            name   = factor(name, levels = rev(unique(name))),
            sig    = !is.na(p_age) & p_age < 0.05,
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
        scale_x_continuous(
            name   = "Odds ratio per year of age (95% CI)",
            breaks = scales::pretty_breaks(n = 6)
        ) +
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
            scale_color_brewer(palette = "Set2", name = unique(results$stratum[1]) |>
                                   (\(.) if (is.null(attr(results, "stratify_by"))) "stratum"
                                    else attr(results, "stratify_by"))()) +
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
extra_covs <- "Smoke_History + Sequenced_gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6"

### UNSTRATIFIED
filter_res <- test_age_association(
    config     = filter_configs,
    all_ch     = all_ch2,
    cov        = cov
)

minad_res <- test_age_association(
    config     = minad_configs,
    all_ch     = all_ch2,
    cov        = cov
)

cov_res <- test_age_association(
    config     = cov_configs,
    all_ch     = all_ch2,
    cov        = cov
)

### STRATIFIED - VAF
filter_vaf_res <- test_age_association(
    config      = filter_configs,
    all_ch      = all_ch2,
    cov         = cov,
    stratify_by = "VAF_Strata"
)

minad_vaf_res <- test_age_association(
    config      = minad_configs,
    all_ch      = all_ch2,
    cov         = cov,
    stratify_by = "VAF_Strata"
)

cov_vaf_res <- test_age_association(
    config      = cov_configs,
    all_ch      = all_ch2,
    cov         = cov,
    stratify_by = "VAF_Strata"
)

### STRATIFIED - FREEZE
filter_freeze_res <- test_age_association(
    config      = filter_configs,
    all_ch      = all_ch2,
    cov         = cov,
    stratify_by = "Batch"
)

minad_freeze_res <- test_age_association(
    config      = minad_configs,
    all_ch      = all_ch2,
    cov         = cov,
    stratify_by = "Batch"
)

cov_freeze_res <- test_age_association(
    config      = cov_configs,
    all_ch      = all_ch2,
    cov         = cov,
    stratify_by = "Batch"
)

# write_xlsx(age_assoc_results, file.path("ch", "data", "ch_filter_age_sensitivity.xlsx"))
# write_xlsx(age_assoc_vaf,     file.path("ch", "data", "ch_filter_age_sensitivity_vaf_strata.xlsx"))
# write_xlsx(age_assoc_batch,   file.path("ch", "data", "ch_filter_age_sensitivity_batch.xlsx"))

# ========================
# PLOTS
# ========================
p <- plot_age_association(filter_res, n_total = nrow(cov))
ggsave(file.path("ch", "figures", "qc_filter_age_sensitivity.pdf"),
       p, width = 7,
       height = 0.4 * n_distinct(age_assoc_results$name) + 2,
       limitsize = FALSE)

p <- plot_age_association(minad_res,  title   = "Age association by VAF stratum", n_total = nrow(cov))
ggsave(file.path("ch", "figures", "qc_filter_age_sensitivity_vaf_strata.pdf"),
       p, width = 5 * n_distinct(age_assoc_vaf$stratum),
       height = 0.4 * n_distinct(age_assoc_vaf$name) + 2,
       limitsize = FALSE)

p <- plot_age_association(cov_res, title   = "Age association by batch", n_total = nrow(cov))
ggsave(file.path("ch", "figures", "qc_filter_age_sensitivity_batch.pdf"),
       p, width = 4 * n_distinct(age_assoc_batch$stratum),
       height = 0.4 * n_distinct(age_assoc_batch$name) + 2,
       limitsize = FALSE)

### STRATIFIED - VAF
p <- plot_age_association(filter_vaf_res, n_total = nrow(cov))
ggsave(file.path("ch", "figures", "qc_filter_age_sensitivity.pdf"),
       p, width = 7,
       height = 0.4 * n_distinct(age_assoc_results$name) + 2,
       limitsize = FALSE)

p <- plot_age_association(minad_vaf_res,  title   = "Age association by VAF stratum", n_total = nrow(cov))
ggsave(file.path("ch", "figures", "qc_filter_age_sensitivity_vaf_strata.pdf"),
       p, width = 5 * n_distinct(age_assoc_vaf$stratum),
       height = 0.4 * n_distinct(age_assoc_vaf$name) + 2,
       limitsize = FALSE)

p <- plot_age_association(cov_vaf_res, title   = "Age association by batch", n_total = nrow(cov))
ggsave(file.path("ch", "figures", "qc_filter_age_sensitivity_batch.pdf"),
       p, width = 4 * n_distinct(age_assoc_batch$stratum),
       height = 0.4 * n_distinct(age_assoc_batch$name) + 2,
       limitsize = FALSE)

### STRATIFIED - FREEZE
p <- plot_age_association(filter_freeze_res, n_total = nrow(cov))
ggsave(file.path("ch", "figures", "qc_filter_age_sensitivity.pdf"),
       p, width = 7,
       height = 0.4 * n_distinct(age_assoc_results$name) + 2,
       limitsize = FALSE)

p <- plot_age_association(minad_freeze_res,  title   = "Age association by VAF stratum", n_total = nrow(cov))
ggsave(file.path("ch", "figures", "qc_filter_age_sensitivity_vaf_strata.pdf"),
       p, width = 5 * n_distinct(age_assoc_vaf$stratum),
       height = 0.4 * n_distinct(age_assoc_vaf$name) + 2,
       limitsize = FALSE)

p <- plot_age_association(cov_freeze_res, title   = "Age association by batch", n_total = nrow(cov))
ggsave(file.path("ch", "figures", "qc_filter_age_sensitivity_batch.pdf"),
       p, width = 4 * n_distinct(age_assoc_batch$stratum),
       height = 0.4 * n_distinct(age_assoc_batch$name) + 2,
       limitsize = FALSE)


# # ========================
# # PLOTTING FUNCTIONS
# # ========================
# plot_category_breakdown <- function(df, title = "CHIP variant categories", filename = NULL, width = 6, height = 5) {
#     cat_summary <- df %>%
#         filter(!is.na(variant_category), variant_category != "non_whitelist") %>%
#         count(variant_category, name = "N") %>%
#         mutate(
#             variant_category = factor(variant_category,
#                                       levels = c("splice", "lof", "missense", "exception", "manual_review")),
#             fill_color = case_when(
#                 variant_category == "manual_review" ~ "#C0392B",
#                 variant_category == "exception"     ~ "#BA7517",
#                 TRUE                                ~ "steelblue"
#             )
#         )
#
#     p <- ggplot(cat_summary, aes(x = variant_category, y = N, fill = variant_category)) +
#         geom_col(width = 0.65) +
#         geom_text(aes(label = paste0(N, "\n(", round(100 * N / sum(N), 1), "%)")),
#                   vjust = -0.4, size = 3.2, color = "gray30") +
#         scale_fill_manual(values = setNames(cat_summary$fill_color, cat_summary$variant_category)) +
#         scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
#         scale_x_discrete(labels = c(
#             "splice"        = "Splice",
#             "lof"           = "LoF / Frameshift",
#             "missense"      = "Missense",
#             "exception"     = "Gene-specific\nexception",
#             "manual_review" = "Manual review"
#         )) +
#         labs(title = title, x = NULL, y = "Number of variants",
#              subtitle = paste0("Total: ", nrow(df))) +
#         theme_minimal(base_size = 12) +
#         theme(legend.position = "none", panel.grid.major.x = element_blank(),
#               panel.grid.minor = element_blank(),
#               plot.title    = element_text(face = "bold", size = 13),
#               plot.subtitle = element_text(size = 10, color = "gray50"),
#               axis.text.x   = element_text(size = 11))
#
#     if (!is.null(filename))
#         ggsave(file.path("ch", "figures", filename), p, width = width, height = height)
#     p
# }
#
# plot_genes <- function(df, title = "CHIP variants by gene", top_n = NULL,
#                        fill = "lightblue", filename = NULL, width = 10, height = 5) {
#     gene_counts <- df %>%
#         filter(!is.na(Gene)) %>%
#         count(Gene, name = "N") %>%
#         arrange(desc(N))
#
#     if (!is.null(top_n)) gene_counts <- slice_max(gene_counts, N, n = top_n)
#
#     gene_counts <- gene_counts %>%
#         mutate(Gene = factor(Gene, levels = (Gene)))  # descending left to right
#
#     p <- ggplot(gene_counts, aes(x = Gene, y = N)) +
#         geom_col(width = 0.65, fill = fill) +
#         geom_text(aes(label = N), vjust = -0.4, size = 3.2, color = "gray30") +
#         scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
#         labs(title = title, x = NULL, y = "Number of variants",
#              subtitle = paste0("Total: ", sum(gene_counts$N),
#                                " | Genes: ", nrow(gene_counts))) +
#         theme_minimal(base_size = 12) +
#         theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
#               plot.title    = element_text(face = "bold", size = 13),
#               plot.subtitle = element_text(size = 10, color = "gray50"),
#               axis.text.x   = element_text(size = 9, angle = 45, hjust = 1))
#
#     if (!is.null(filename))
#         ggsave(file.path("ch", "figures", filename), p, width = width, height = height)
#     p
# }
#
# plot_gene_exceptions <- function(df, title = "Gene-specific exceptions", filename = NULL, width = 4, height = 4) {
#     gene_exc <- df %>%
#         filter(wl.exception) %>%
#         count(Gene, name = "N") %>%
#         # arrange(desc(N)) %>%
#         mutate(Gene = factor(Gene, levels = Gene[order(-N)]))
#
#     p <- ggplot(gene_exc, aes(x = Gene, y = N)) +
#         geom_col(width = 0.6, fill = "#BA7517") +
#         geom_text(aes(label = N), vjust = -0.4, size = 3.2, color = "gray30") +
#         scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
#         labs(title = title, x = NULL, y = "Number of variants") +
#         theme_minimal(base_size = 12) +
#         theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
#               plot.title = element_text(face = "bold", size = 12),
#               axis.text.x = element_text(size = 11, angle = 45, hjust = 1))
#
#     if (!is.null(filename))
#         ggsave(file.path("ch", "figures", filename), p, width = width, height = height)
#     p
# }
#
# # ========================
# # SUMMARY FIGUREts
# # ========================
# # write.csv(results, file.path("ch", "data", "ch_artifact_test_results.csv"), row.names = FALSE)
#
# # cat_summary <- all_ch_clean %>%
# #     filter(variant_category != "non_whitelist") %>%
# #     count(variant_category, name = "N") %>%
# #     mutate(
# #         variant_category = factor(variant_category,
# #                                   levels = c("splice", "lof", "missense",
# #                                              "exception", "manual_review")),
# #         fill_color = case_when(
# #             variant_category == "manual_review" ~ "#C0392B",
# #             variant_category == "exception"     ~ "#BA7517",
# #             TRUE                                ~ "steelblue"
# #         )
# #     )
#
# # rs <- read_tsv(file.path("ch", "data", "rs7705526_output.raw"))
# # all_samples_cov <- cov %>%
# #     rename(Sample.ID = person_id) %>%
# #     select(Sample.ID, Age = all_of("Sample_age"))
#
#
#
# # # ========================
# # # FOREST PLOT
# # # ========================
# # results_plot <- results %>%
# #     left_join(
# #         all_ch %>% select(variant_id, HGVSp) %>% distinct(),
# #         by = "variant_id"
# #     ) %>%
# #     mutate(
# #         label = ifelse(!is.na(HGVSp) & HGVSp != "" & HGVSp != ".",
# #                        paste0(Gene, ":", HGVSp),
# #                        variant_id),
# #         label = factor(label, levels = label[order(OR)])
# #     )
# #
# # p_artifact <- ggplot(results_plot, aes(x = OR, y = label, color = status)) +
# #     geom_vline(xintercept = 1, linetype = "dotted", color = "gray40") +
# #     geom_errorbar(aes(xmin = CI_lo, xmax = CI_hi), height = 0.3, linewidth = 0.6, orientation = "y") +
# #     geom_point(size = 2.5, shape = 15) +
# #     scale_color_manual(
# #         values = c("retained" = "steelblue", "artifact" = "firebrick"),
# #         labels = c("retained" = "Associated with age (retained)",
# #                    "artifact" = "Not associated (flagged)")
# #     ) +
# #     labs(
# #         title    = "Association of common variants with age",
# #         subtitle = paste0("Variants in ≥", FREQ_THRESHOLD, " individuals | P < 0.10 threshold"),
# #         x = "Odds ratio for association with age",
# #         y = NULL, color = NULL
# #     ) +
# #     theme_minimal() +
# #     theme(
# #         legend.position    = "bottom",
# #         axis.text.y        = element_text(size = 7),
# #         panel.grid.major.y = element_blank()
# #     )
# #
# # ggsave(file.path("ch", "figures", "artifact_age_association.pdf"),
# #        plot = p_artifact, width = 8,
# #        height = max(5, nrow(results_plot) * 0.3), dpi = 300)
