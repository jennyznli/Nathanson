# ========================
# Age association for artifacts
# binomial testing for germline variants
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
dim(cov)
# 3004

# ========================
# IDENTIFY HIGH FREQ VARIANTS
# ========================
FREQ_THRESHOLD <- 3

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
    mutate(potential_artifact = as.integer(n_carriers >= FREQ_THRESHOLD))
write_xlsx(variant_counts, file.path("ch", "data", "ch_var_counts.xlsx"))

high_freq <- variant_counts %>% filter(potential_artifact == 1)

cat(sprintf("High frequency variants to test: %d\n", nrow(high_freq)))
# 25

# ========================
# FREEZE ENRICHMENT FLAG (Fisher's exact test)
# ========================
batch_totals <- all_ch %>%
    group_by(Batch) %>%
    summarise(n_total = n_distinct(Sample.ID), .groups = "drop")

n_batch_1 <- batch_totals$n_total[batch_totals$Batch == 1]
n_batch_2 <- batch_totals$n_total[batch_totals$Batch == 2]

high_freq2 <- high_freq %>%
    mutate(
        p_freeze         = mapply(function(b1, b2) {
            fisher.test(matrix(
                c(b1, n_batch_1 - b1,
                  b2, n_batch_2 - b2),
                nrow = 2
            ))$p.value
        }, batch_1, batch_2),
        p_freeze_adj         = p.adjust(p_freeze, method = "BH"),
        flag_freeze_enriched = p_freeze_adj < 0.05
    )

cat(sprintf("High frequency variants to test: %d\n", nrow(high_freq2)))
cat(sprintf("  of which freeze-enriched:      %d\n", sum(high_freq2$flag_freeze_enriched, na.rm = TRUE)))
# 3

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
# AGE ASSOCIATION FLAG
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

# base model
high_freq_base <- high_freq2 %>%
    rowwise() %>%
    mutate(res = list(artifact_test_age(variant_id, cov))) %>%
    unnest(res) %>% ungroup() %>%
    mutate(
        p_age_adj           = p.adjust(p_age, method = "BH"),
        flag_age_associated = !is.na(p_age) & !separation & p_age < 0.10
    )

# with extra covariates
extra_covs <- "Smoke_History + Sequenced_gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6"

high_freq_full <- high_freq2 %>%
    rowwise() %>%
    mutate(res = list(artifact_test_age(variant_id, cov, extra_covs = extra_covs))) %>%
    unnest(res) %>% ungroup() %>%
    mutate(
        p_age_adj           = p.adjust(p_age, method = "BH"),
        flag_age_associated = !is.na(p_age) & !separation & p_age < 0.10
    )

write_xlsx(high_freq_full, file.path("ch", "data", "ch_high_freq_counts.xlsx"))

# ========================
# JOIN FLAGS BACK ONTO all_ch
# ========================
all_ch2 <- all_ch %>%
    left_join(
        high_freq_full %>% select(variant_id, potential_artifact, p_freeze, p_freeze_adj, flag_freeze_enriched,
                             p_age, p_age_adj, OR, CI_lo, CI_hi,
                             flag_age_associated, separation),
        by = "variant_id"
    ) %>%
    mutate(
        potential_artifact   = replace_na(potential_artifact,   0L),
        flag_age_associated  = replace_na(flag_age_associated,  FALSE),
        flag_freeze_enriched = replace_na(flag_freeze_enriched, FALSE),

        flag_artifact = case_when(
            potential_artifact & flag_freeze_enriched          ~ TRUE,
            potential_artifact & !is.na(p_age) & !flag_age_associated ~ TRUE,
            TRUE                          ~ FALSE
        )
    )

cat(sprintf("Variants flagged as freeze-enriched: %d\n", sum(all_ch2$flag_freeze_enriched)))
cat(sprintf("Variants flagged no age association: %d\n",
            sum(all_ch2$potential_artifact == 1 & !all_ch2$flag_age_associated, na.rm = TRUE)))
cat(sprintf("Variants flagged artifact (either):  %d\n", sum(all_ch2$flag_artifact, na.rm = TRUE)))
cat(sprintf("Unique individuals affected:         %d\n", n_distinct(all_ch2$Sample.ID[all_ch2$flag_artifact])))
# > cat(sprintf("Variants flagged as freeze-enriched: %d\n", sum(all_ch2$flag_freeze_enriched)))
# Variants flagged as freeze-enriched: 167
# > cat(sprintf("Variants flagged no age association: %d\n", sum(!all_ch2$flag_age_associated, na.rm = TRUE)))
# Variants flagged no age association: 206
# > cat(sprintf("Variants flagged artifact (either):  %d\n", sum(all_ch2$flag_artifact, na.rm = TRUE)))
# Variants flagged artifact (either):  253
# > cat(sprintf("Unique individuals affected:         %d\n", n_distinct(all_ch2$Sample.ID[all_ch2$flag_artifact])))
# Unique individuals affected:         230

write_xlsx(all_ch2, file.path("ch", "data", "ch_seq_wl_art_vars.xlsx"))

# ========================
# GERMLINE TESTING
# ========================
all_ch3 <- all_ch2 %>% mutate(
    high_gnomad = ifelse(gnomAD.MAX_AF > 0.001, TRUE, FALSE)
)
sum(all_ch3$high_gnomad)
# 81

# ========================
# GERMLINE TESTING (binomial test)
# ========================
all_ch3 <- all_ch3 %>% mutate(
    VAF_Strata = case_when(
        Sample.AltFrac >= 0.02 & Sample.AltFrac < 0.05 ~ "2-5",
        Sample.AltFrac >= 0.05 & Sample.AltFrac < 0.10 ~ "5-10",
        Sample.AltFrac >= 0.10 ~ "10+",
        TRUE ~ NA_character_
    )
)

# Tests whether observed VAF is consistent with germline het (expected = 0.5)
# Flags variants where we cannot reject germline origin at P < 0.01
# tet <- variant_counts %>% filter(Gene %in% c("TET2", "CBL"))
# examine_germ <- rbind(tet, high_freq)

binomial_germline_test <- function(alt_count, total_depth, p_germline = 0.5) {
    if (is.na(alt_count) || is.na(total_depth) || total_depth == 0) {
        return(data.frame(p_binom = NA_real_, flag_germline_binom = NA))
    }
    # two-sided: is VAF significantly different from 0.5?
    p <- binom.test(x = alt_count, n = total_depth,
                    p = p_germline, alternative = "two.sided")$p.value
    data.frame(
        p_binom            = p,
        flag_germline_binom = p >= 0.01   # fails to reject germline = suspicious
    )
}

# apply to TET2/CBL missense + high_freq variants
germline_test_vars <- all_ch3 %>%
    filter(
        (Gene %in% c("TET2", "CBL") & variant_category == "missense") |
            potential_artifact == 1
    ) %>%
    distinct(variant_id, Sample.ID, Sample.AltDepth, Sample.Depth) %>%   # adjust col names to match yours
    rowwise() %>%
    mutate(res = list(binomial_germline_test(Sample.AltDepth, Sample.Depth))) %>%
    unnest(res) %>%
    ungroup()

# summarize per variant: flag if MAJORITY of carriers fail binomial test
germline_summary <- germline_test_vars %>%
    group_by(variant_id) %>%
    summarise(
        n_tested          = sum(!is.na(p_binom)),
        n_fail_binom      = sum(flag_germline_binom, na.rm = TRUE),
        prop_fail_binom   = n_fail_binom / n_tested,
        flag_germline     = prop_fail_binom > 0.5,   # majority of carriers look germline
        .groups           = "drop"
    )

cat(sprintf("Variants tested for germline: %d\n",    nrow(germline_summary)))
# 48
cat(sprintf("Variants flagged germline:    %d\n",    sum(germline_summary$flag_germline, na.rm = TRUE)))
# 6

# join back onto all_ch3
all_ch4 <- all_ch3 %>%
    left_join(germline_summary %>% select(variant_id, n_fail_binom, prop_fail_binom, flag_germline),
              by = "variant_id") %>%
    mutate(
        flag_germline = replace_na(flag_germline, FALSE),
        # combined final artifact/germline flag
        flag_remove = flag_artifact | flag_germline | high_gnomad
    )

cat(sprintf("Variants flagged artifact:          %d\n", sum(all_ch4$flag_artifact,  na.rm = TRUE)))
cat(sprintf("Variants flagged germline (binom):  %d\n", sum(all_ch4$flag_germline,  na.rm = TRUE)))
cat(sprintf("Variants flagged high gnomAD AF:    %d\n", sum(all_ch4$high_gnomad,    na.rm = TRUE)))
cat(sprintf("Variants flagged for removal (any): %d\n", sum(all_ch4$flag_remove,    na.rm = TRUE)))
cat(sprintf("Unique individuals affected:        %d\n", n_distinct(all_ch4$Sample.ID[all_ch4$flag_remove])))

write_xlsx(all_ch4, file.path("ch", "data", "ch_seq_wl_art_germ_vars.xlsx"))

# ========================
# COMPLEX INDELS
# ========================
all_ch4 <- all_ch4 %>% mutate(.row_idx = row_number())

cluster_idx <- all_ch4 %>%
    filter(whitelist) %>%
    group_by(Sample.ID, Gene, Chr) %>%
    filter(n() > 1) %>%
    arrange(Sample.ID, Gene, Chr, Start) %>%
    mutate(near_prev = c(FALSE, diff(Start) <= 50)) %>%
    filter(near_prev | lead(near_prev, default = FALSE)) %>%
    ungroup() %>%
    pull(.row_idx)

all_ch4$cluster[all_ch4$.row_idx %in% cluster_idx] <- TRUE
# all_ch$manualreview[all_ch4$.row_idx %in% cluster_idx] <- TRUE
all_ch5 <- all_ch4 %>% select(-.row_idx)
write_xlsx(all_ch4, file.path("ch", "data", "ch_seq_wl_art_germ_cluster_vars.xlsx"))

# ========================
# MIN AD - VARIOUS THRESHOLDS
# ========================
all_ch_min4 <- all_ch5 %>% filter(Sample.AltDepth >= 4)
all_ch_min5 <- all_ch5 %>% filter(Sample.AltDepth >= 5)

# test age association
# overall
# vs diff strata

# ========================
# TRY DIFF FILTERS + EVALUATE
# ========================





# ========================
# REMOVE ARTIFACTS AND HIGH VAF - MANUAL
# ========================
cat("Variants before artifact removal:", nrow(all_ch), "\n")
cat("Variants after artifact removal: ", nrow(all_ch3), "\n")

write_xlsx(all_ch_clean, file.path("ch", "data", "ch_wl_art_germ.xlsx"))

plot_category_breakdown(all_ch_clean, title = "CHIP variant categories (artifact-filtered)",
                        filename = "ch_wl_art_categories.pdf")
plot_genes(all_ch_clean, title = "Variants per gene (artifact-filtered)",
           filename = "ch_wl_art_genes.pdf", width = 12)
plot_genes(all_ch_clean, top_n = 10, title = "Top 10 genes (artifact-filtered)",
           filename = "ch_wl_art_genes_top10.pdf", width = 7)

# ========================
# PLOTTING FUNCTIONS
# ========================
plot_category_breakdown <- function(df, title = "CHIP variant categories", filename = NULL, width = 6, height = 5) {
    cat_summary <- df %>%
        filter(!is.na(variant_category), variant_category != "non_whitelist") %>%
        count(variant_category, name = "N") %>%
        mutate(
            variant_category = factor(variant_category,
                                      levels = c("splice", "lof", "missense", "exception", "manual_review")),
            fill_color = case_when(
                variant_category == "manual_review" ~ "#C0392B",
                variant_category == "exception"     ~ "#BA7517",
                TRUE                                ~ "steelblue"
            )
        )

    p <- ggplot(cat_summary, aes(x = variant_category, y = N, fill = variant_category)) +
        geom_col(width = 0.65) +
        geom_text(aes(label = paste0(N, "\n(", round(100 * N / sum(N), 1), "%)")),
                  vjust = -0.4, size = 3.2, color = "gray30") +
        scale_fill_manual(values = setNames(cat_summary$fill_color, cat_summary$variant_category)) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
        scale_x_discrete(labels = c(
            "splice"        = "Splice",
            "lof"           = "LoF / Frameshift",
            "missense"      = "Missense",
            "exception"     = "Gene-specific\nexception",
            "manual_review" = "Manual review"
        )) +
        labs(title = title, x = NULL, y = "Number of variants",
             subtitle = paste0("Total: ", nrow(df))) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "none", panel.grid.major.x = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title    = element_text(face = "bold", size = 13),
              plot.subtitle = element_text(size = 10, color = "gray50"),
              axis.text.x   = element_text(size = 11))

    if (!is.null(filename))
        ggsave(file.path("ch", "figures", filename), p, width = width, height = height)
    p
}

plot_genes <- function(df, title = "CHIP variants by gene", top_n = NULL,
                       fill = "lightblue", filename = NULL, width = 10, height = 5) {
    gene_counts <- df %>%
        filter(!is.na(Gene)) %>%
        count(Gene, name = "N") %>%
        arrange(desc(N))

    if (!is.null(top_n)) gene_counts <- slice_max(gene_counts, N, n = top_n)

    gene_counts <- gene_counts %>%
        mutate(Gene = factor(Gene, levels = (Gene)))  # descending left to right

    p <- ggplot(gene_counts, aes(x = Gene, y = N)) +
        geom_col(width = 0.65, fill = fill) +
        geom_text(aes(label = N), vjust = -0.4, size = 3.2, color = "gray30") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
        labs(title = title, x = NULL, y = "Number of variants",
             subtitle = paste0("Total: ", sum(gene_counts$N),
                               " | Genes: ", nrow(gene_counts))) +
        theme_minimal(base_size = 12) +
        theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
              plot.title    = element_text(face = "bold", size = 13),
              plot.subtitle = element_text(size = 10, color = "gray50"),
              axis.text.x   = element_text(size = 9, angle = 45, hjust = 1))

    if (!is.null(filename))
        ggsave(file.path("ch", "figures", filename), p, width = width, height = height)
    p
}

plot_gene_exceptions <- function(df, title = "Gene-specific exceptions", filename = NULL, width = 4, height = 4) {
    gene_exc <- df %>%
        filter(wl.exception) %>%
        count(Gene, name = "N") %>%
        # arrange(desc(N)) %>%
        mutate(Gene = factor(Gene, levels = Gene[order(-N)]))

    p <- ggplot(gene_exc, aes(x = Gene, y = N)) +
        geom_col(width = 0.6, fill = "#BA7517") +
        geom_text(aes(label = N), vjust = -0.4, size = 3.2, color = "gray30") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
        labs(title = title, x = NULL, y = "Number of variants") +
        theme_minimal(base_size = 12) +
        theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
              plot.title = element_text(face = "bold", size = 12),
              axis.text.x = element_text(size = 11, angle = 45, hjust = 1))

    if (!is.null(filename))
        ggsave(file.path("ch", "figures", filename), p, width = width, height = height)
    p
}

# ========================
# SUMMARY FIGUREts
# ========================
# write.csv(results, file.path("ch", "data", "ch_artifact_test_results.csv"), row.names = FALSE)

# cat_summary <- all_ch_clean %>%
#     filter(variant_category != "non_whitelist") %>%
#     count(variant_category, name = "N") %>%
#     mutate(
#         variant_category = factor(variant_category,
#                                   levels = c("splice", "lof", "missense",
#                                              "exception", "manual_review")),
#         fill_color = case_when(
#             variant_category == "manual_review" ~ "#C0392B",
#             variant_category == "exception"     ~ "#BA7517",
#             TRUE                                ~ "steelblue"
#         )
#     )

# rs <- read_tsv(file.path("ch", "data", "rs7705526_output.raw"))
# all_samples_cov <- cov %>%
#     rename(Sample.ID = person_id) %>%
#     select(Sample.ID, Age = all_of("Sample_age"))



# # ========================
# # FOREST PLOT
# # ========================
# results_plot <- results %>%
#     left_join(
#         all_ch %>% select(variant_id, HGVSp) %>% distinct(),
#         by = "variant_id"
#     ) %>%
#     mutate(
#         label = ifelse(!is.na(HGVSp) & HGVSp != "" & HGVSp != ".",
#                        paste0(Gene, ":", HGVSp),
#                        variant_id),
#         label = factor(label, levels = label[order(OR)])
#     )
#
# p_artifact <- ggplot(results_plot, aes(x = OR, y = label, color = status)) +
#     geom_vline(xintercept = 1, linetype = "dotted", color = "gray40") +
#     geom_errorbar(aes(xmin = CI_lo, xmax = CI_hi), height = 0.3, linewidth = 0.6, orientation = "y") +
#     geom_point(size = 2.5, shape = 15) +
#     scale_color_manual(
#         values = c("retained" = "steelblue", "artifact" = "firebrick"),
#         labels = c("retained" = "Associated with age (retained)",
#                    "artifact" = "Not associated (flagged)")
#     ) +
#     labs(
#         title    = "Association of common variants with age",
#         subtitle = paste0("Variants in ≥", FREQ_THRESHOLD, " individuals | P < 0.10 threshold"),
#         x = "Odds ratio for association with age",
#         y = NULL, color = NULL
#     ) +
#     theme_minimal() +
#     theme(
#         legend.position    = "bottom",
#         axis.text.y        = element_text(size = 7),
#         panel.grid.major.y = element_blank()
#     )
#
# ggsave(file.path("ch", "figures", "artifact_age_association.pdf"),
#        plot = p_artifact, width = 8,
#        height = max(5, nrow(results_plot) * 0.3), dpi = 300)

# ========================
# EXAMINE ASXL1
# ========================
x <- read.csv(file.path("ch", "data", "ch_all_vars.csv"), row.names = 1)
y <- read_excel(file.path("ch", "data", "ch_wl_final.xlsx"))
z <- read_excel(file.path("ch", "data", "ch_wl_art.xlsx"))

View(x %>% filter(Gene == "ASXL1"))
View(y %>% filter(Gene == "ASXL1"))
View(z %>% filter(Gene == "ASXL1"))

