# ========================
# Age association for artifacts
# ========================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
library(data.table, quietly = T)
library(ggVennDiagram)
library(UpSetR)
library(ggrepel)

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
all_ch <- read_excel(file.path("ch", "data", "ch_seq_wl_vars.xlsx")) %>% filter(keep == 1)
dim(all_ch)
# 891

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

# ========================
# 1. HIGH FREQ VARIANTS
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
    mutate(flag_high_freq = (n_carriers >= FREQ_THRESHOLD))

dim(variant_counts)
# 607 unique variants

high_freq <- variant_counts %>% filter(flag_high_freq)
high_freq_var <- unique(high_freq$variant_id)

length(high_freq_var)
# 20

high_freq_vars <- all_ch %>% filter(variant_id %in% high_freq_var) %>% arrange(variant_id)
dim(high_freq_vars)
# 263

# ========================
# 2. FREEZE ENRICHMENT
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

length(freeze_enriched_var)
# 2

freeze_removed_var <- all_ch %>% filter(variant_id %in% freeze_enriched_var) %>% arrange(variant_id)
dim(freeze_removed_var)
# 149

# # ========================
# # 3. AGE ASSOCIATION
# # ========================
# run_age_tests <- function(data, variants) {
#     variants %>%
#         rowwise() %>%
#         mutate(res = list(artifact_test_age(variant_id, data, extra_covs = extra_covs))) %>%
#         unnest(res) %>%
#         ungroup() %>%
#         mutate(
#             p_age_adj = p.adjust(p_age, method = "BH"),
#             flag_age_associated = !is.na(p_age) & !separation & p_age < AGE_SIG
#         )
# }
#
# # test age association per variant category
# artifact_test_age <- function(vkey, vars, extra_covs = NULL) {
#     carriers <- vars %>%
#         filter(variant_id == vkey) %>%
#         distinct(Sample.ID) %>%
#         mutate(carrier = 1)
#
#     df <- cov %>%
#         left_join(carriers, by = c("person_id" = "Sample.ID")) %>%
#         mutate(carrier = ifelse(is.na(carrier), 0, carrier))
#
#     base    <- "carrier ~ Sample_age + as.factor(Batch)"
#     formula <- if (!is.null(extra_covs)) {
#         as.formula(paste(base, "+", extra_covs))
#     } else {
#         as.formula(base)
#     }
#
#     tryCatch({
#         m <- withCallingHandlers(
#             glm(formula, data = df, family = binomial()),
#             warning = function(w) {
#                 if (grepl("fitted probabilities numerically 0 or 1", conditionMessage(w)))
#                     invokeRestart("muffleWarning")
#             }
#         )
#         separation <- any(fitted(m) %in% c(0, 1)) | !m$converged
#         s  <- coef(summary(m))
#         ci <- confint.default(m)
#         data.frame(
#             OR         = exp(coef(m)["Sample_age"]),
#             CI_lo      = exp(ci["Sample_age", 1]),
#             CI_hi      = exp(ci["Sample_age", 2]),
#             p_age      = s["Sample_age", "Pr(>|z|)"],
#             separation = separation
#         )
#     }, error = function(e) {
#         warning(sprintf("Model failed for %s: %s", vkey, conditionMessage(e)))
#         data.frame(OR = NA_real_, CI_lo = NA_real_, CI_hi = NA_real_,
#                    p_age = NA_real_, separation = NA)
#     })
# }
#
# age_association <- high_freq2 %>%
#     rowwise() %>%
#     mutate(res = list(artifact_test_age(variant_id, all_ch, extra_covs = extra_covs))) %>%
#     unnest(res) %>% ungroup() %>%
#     mutate(
#         p_age_adj = p.adjust(p_age, method = "BH"),
#         flag_age_associated = !is.na(p_age) & !separation & p_age < AGE_SIG
#     )
#
# not_age_associated <- age_association %>% filter(!flag_age_associated)
# not_age_associated_var <- not_age_associated %>% pull(variant_id)
# length(not_age_associated_var)
# # 21
#
# age_removed_var <- all_ch %>% filter(variant_id %in% not_age_associated_var) %>% arrange(variant_id)
# dim(age_removed_var)
# # 265

# ========================
# 4. GERMLINE TESTING
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
dim(fail_majority_germline)
# 51 - most of these are singletons...

# from all vars tested, remove those individuals' variants failing binomial
# as well as entire variant categories whose majority failed
germline_fail_vars <- germline_test_vars %>%
    filter(variant_id %in% fail_majority_germline$variant_id) %>%
    filter(flag_germline_ind)

germline_removed_vars <- all_ch %>% inner_join(germline_fail_vars, by = c("Sample.ID", "variant_id"))
dim(germline_removed_vars)
# 84

### PLOT ###
# colored by TET2 missense
tet2_vars <- germline_test_vars %>% filter(Gene == "TET2")

p <- ggplot(tet2_vars,
            aes(x = Sample_age, y = Sample.AltFrac,
                color = ProteinChange_1L,
                size  = !flag_germline_ind)) +
    geom_point(alpha = 0.85) +
    geom_text_repel(
        data = tet2_vars %>% filter(flag_germline_ind),
        aes(label = ProteinChange_1L),
        size          = 3,
        fontface      = "bold",
        max.overlaps  = 20,
        box.padding   = 0.4,
        point.padding = 0.3,
        segment.color = "grey60",
        segment.size  = 0.3,
        show.legend   = FALSE
    ) +
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

# ========================
# 5. CLUSTERS
# ========================
# annotate with cluster cols
all_ch_clustered <- all_ch %>%
    mutate(.row_idx = row_number()) %>%
    group_by(Sample.ID, Gene, Chr) %>%
    arrange(Sample.ID, Gene, Chr, Start) %>%
    mutate(
        near_prev  = c(FALSE, diff(Start) <= CLUSTER_THRESHOLD),
        in_cluster = near_prev | lead(near_prev, default = FALSE)
    ) %>%
    ungroup()

# just rows that are in clusters
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

# summary by cluster
cluster_summary <- cluster_all %>%
    group_by(cluster_id) %>%
    summarise(
        Sample.ID        = first(Sample.ID),
        Gene             = first(Gene),
        Chr              = first(Chr),
        Start_min        = min(Start),
        Start_max        = max(Start),
        span_bp          = Start_max - Start_min,
        n_variants       = n(),
        variant_ids      = paste(variant_id, collapse = "; "),
        protein_changes  = paste(unique(ProteinChange_1L), collapse = "; "),
        # has_dbsnp        = any(grepl("\\brs[0-9]", Existing.variation)),
        # has_cosmic       = any(grepl("\\bCO[0-9]", Existing.variation)),
        mean_VAF         = round(mean(Sample.AltFrac, na.rm = TRUE), 4),
        .groups          = "drop"
    ) %>%
    arrange(Sample.ID, Chr, Start_min)
dim(cluster_summary)
# 99 unique clusters

# get one representative variant per cluster
cluster_rep <- cluster_all %>%
    group_by(cluster_id) %>%
    arrange(!grepl("\\bC[0-9]", Existing.variation), Start) %>%
    slice(1) %>%
    ungroup()

# filter out all other variants not representative
cluster_not_rep_var <- cluster_all %>%
    filter(!(variant_id %in% cluster_rep$variant_id)) %>%
    pull(variant_id)

# ========================
# JOIN BACK TO ALL_CH - ALL FLAGS
# ========================
all_ch2 <- all_ch %>%
    mutate(
        flag_high_freq          = variant_id %in% high_freq_var,
        flag_freeze_enriched    = variant_id %in% freeze_enriched_var,
        flag_not_age_associated = variant_id %in% not_age_associated_var,
        flag_germline_ind       = paste(Sample.ID, variant_id) %in%
            paste(germline_removed_vars$Sample.ID, germline_removed_vars$variant_id),
        flag_cluster            = variant_id %in% cluster_not_rep_var,
        flag_full               =  flag_freeze_enriched | flag_not_age_associated | flag_germline_ind | flag_cluster
    )
all_ch2 <- all_ch2 %>% mutate(n_flags = rowSums(across(c("flag_freeze_enriched",
                                      "flag_not_age_associated", "flag_germline_ind",
                                      "flag_cluster"), as.integer), na.rm = TRUE))

flags <- c("flag_high_freq", "flag_freeze_enriched", "flag_not_age_associated",
           "flag_germline_ind", "flag_cluster", "flag_full")

# summary of how many rows each flag removes
flag_summary <- tibble(
    flag                    = flags,
    n_rows_flagged          = c(
        sum(all_ch2$flag_high_freq),
        sum(all_ch2$flag_freeze_enriched),
        sum(all_ch2$flag_not_age_associated),
        sum(all_ch2$flag_germline_ind),
        sum(all_ch2$flag_cluster),
        sum(all_ch2$flag_full)
    ),
    n_individuals_flagged   = c(
        n_distinct(all_ch2$Sample.ID[all_ch2$flag_high_freq]),
        n_distinct(all_ch2$Sample.ID[all_ch2$flag_freeze_enriched]),
        n_distinct(all_ch2$Sample.ID[all_ch2$flag_not_age_associated]),
        n_distinct(all_ch2$Sample.ID[all_ch2$flag_germline_ind]),
        n_distinct(all_ch2$Sample.ID[all_ch2$flag_cluster]),
        n_distinct(all_ch2$Sample.ID[all_ch2$flag_full])
    ),
    n_variants_flagged      = c(
        n_distinct(all_ch2$variant_id[all_ch2$flag_high_freq]),
        n_distinct(all_ch2$variant_id[all_ch2$flag_freeze_enriched]),
        n_distinct(all_ch2$variant_id[all_ch2$flag_not_age_associated]),
        n_distinct(all_ch2$variant_id[all_ch2$flag_germline_ind]),
        n_distinct(all_ch2$variant_id[all_ch2$flag_cluster]),
        n_distinct(all_ch2$variant_id[all_ch2$flag_full])
    )
)
flag_summary
# flag                    n_rows_flagged n_individuals_flagged n_variants_flagged
# <chr>                            <int>                 <int>              <int>
#     1 flag_high_freq                     329                   287                 25
# 2 flag_freeze_enriched               149                   148                  2
# 3 flag_not_age_associated            265                   241                 21
# 4 flag_germline_ind                   84                    82                 51
# 5 flag_cluster                       147                    80                114
# 6 flag_full                          500                   387                180

cat("\nRows by number of flags triggered:\n")
print(table(all_ch2$n_flags))
# 0   1   2
# 499 355 145

# final clean dataset
all_ch_clean <- all_ch2 %>% filter(!flag_full)
cat(sprintf("Rows before filtering: %d\n", nrow(all_ch2)))
# 999
cat(sprintf("Rows after filtering:  %d\n", nrow(all_ch_clean)))
# 499
cat(sprintf("Rows removed:          %d\n", nrow(all_ch2) - nrow(all_ch_clean)))
# 500

write_xlsx(all_ch2, file.path("ch", "data", "ch_seq_wl_flags_vars.xlsx"))
write_xlsx(all_ch_clean, file.path("ch", "data", "ch_seq_wl_art_vars.xlsx"))
write_xlsx(flag_summary, file.path("ch", "data", "ch_flags_summary.xlsx"))

# ========================
# CHIP VARIANT COUNT PER INDIVIDUAL
# ========================
chip_per_ind <- all_ch_clean %>%
    group_by(Sample.ID) %>%
    summarise(n_chip = n_distinct(variant_id), .groups = "drop")

chip_counts <- chip_per_ind %>%
    count(n_chip) %>%
    mutate(n_chip = factor(n_chip))

p_chip <- ggplot(chip_counts, aes(x = n_chip, y = n)) +
    geom_col(fill = "#378ADD", width = 0.6) +
    geom_text(aes(label = n), vjust = -0.4, size = 3.5, fontface = "bold") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
        title    = "Number of CH variants per individual",
        subtitle = sprintf("n = %d individuals with CH", n_distinct(chip_per_ind$Sample.ID)),
        x        = "Number of CH variants",
        y        = "Number of individuals"
    ) +
    theme_minimal(base_size = 11) +
    theme(
        plot.title         = element_text(face = "bold", hjust = 0.5),
        plot.subtitle      = element_text(hjust = 0.5, size = 9, color = "grey40"),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank()
    )

ggsave(file.path("ch", "figures", "ch_variants_per_individual.pdf"),
       p_chip, width = 6, height = 5)

# ========================
# FLAGS UPSET PLOT
# =======================
upset_flags <- c("flag_freeze_enriched", "flag_not_age_associated",
           "flag_germline_ind", "flag_cluster")

flag_mat <- all_ch2 %>%
    select(all_of(upset_flags)) %>%
    mutate(across(everything(), as.integer))

pdf(file.path("ch", "figures", "qc_flag_overlap_upset.pdf"), width = 10, height = 5)
upset(
    as.data.frame(flag_mat),
    sets           = upset_flags,
    order.by       = "freq",
    sets.bar.color = "#378ADD",
    main.bar.color = "#1D9E75",
    text.scale     = 1.2,
    mb.ratio       = c(0.6, 0.4)
)
dev.off()

# ========================
# COMPILE RESULTS FOR MANUAL REVIEW
# ========================
flag_overview <- tibble(
    filter_step        = c("1. High frequency",
                           "2. Freeze enriched",
                           "3. Not age associated",
                           "4. Germline (majority fail)",
                           "5. Cluster (non-representative)"),
    n_variants_tested  = c(nrow(variant_counts),
                           nrow(high_freq2),
                           nrow(age_association),
                           nrow(germline_summary),
                           n_distinct(cluster_all$variant_id)),
    n_variants_flagged = c(sum(variant_counts$flag_high_freq),
                           sum(high_freq2$flag_freeze_enriched, na.rm = TRUE),
                           sum(!age_association$flag_age_associated &
                                   !age_association$separation, na.rm = TRUE),
                           sum(germline_summary$flag_germline_var, na.rm = TRUE),
                           length(cluster_not_rep_var)),
    n_rows_flagged     = c(sum(all_ch2$flag_high_freq),
                           sum(all_ch2$flag_freeze_enriched),
                           sum(all_ch2$flag_not_age_associated),
                           sum(all_ch2$flag_germline_ind),
                           sum(all_ch2$flag_cluster)),
    n_individuals      = c(n_distinct(all_ch2$Sample.ID[all_ch2$flag_high_freq]),
                           n_distinct(all_ch2$Sample.ID[all_ch2$flag_freeze_enriched]),
                           n_distinct(all_ch2$Sample.ID[all_ch2$flag_not_age_associated]),
                           n_distinct(all_ch2$Sample.ID[all_ch2$flag_germline_ind]),
                           n_distinct(all_ch2$Sample.ID[all_ch2$flag_cluster])),
    threshold          = c(sprintf("n_carriers >= %d", FREQ_THRESHOLD),
                           sprintf("p_freeze_adj < %.2f", FREEZE_SIG),
                           sprintf("p_age < %.2f (among high-freq variants)", AGE_SIG),
                           sprintf("prop_fail_binom > 0.50, p_binom >= %.2f", GERMLINE_SIG),
                           sprintf("within %d bp, same Sample/Gene/Chr", CLUSTER_THRESHOLD))
)

all_variants_master <- variant_counts %>%
    left_join(
        high_freq2 %>% select(variant_id, p_freeze, p_freeze_adj, flag_freeze_enriched),
        by = "variant_id"
    ) %>%
    left_join(
        age_association %>% select(variant_id, OR, CI_lo, CI_hi,
                                   p_age, p_age_adj, separation, flag_age_associated),
        by = "variant_id"
    ) %>%
    left_join(
        germline_summary %>% select(variant_id, n_tested, n_fail_binom,
                                    prop_fail_binom, flag_germline_var),
        by = "variant_id"
    ) %>%
    left_join(
        all_ch %>%
            group_by(variant_id) %>%
            summarise(
                across(c(whitelist, wl.mis, wl.lof, wl.splice, wl.exception), first),
                wl_manualreview = any(manualreview, na.rm = TRUE),
                .groups = "drop"
            ),
        by = "variant_id"
    ) %>%
    mutate(
        flag_cluster            = variant_id %in% cluster_not_rep_var,
        flag_high_freq          = flag_high_freq == 1,
        flag_freeze_enriched    = ifelse(is.na(flag_freeze_enriched), FALSE, flag_freeze_enriched),
        flag_not_age_associated = variant_id %in% not_age_associated_var,
        flag_germline_var       = ifelse(is.na(flag_germline_var), FALSE, flag_germline_var),
        n_flags                 = flag_high_freq + flag_freeze_enriched +
            flag_not_age_associated + flag_germline_var + flag_cluster,
        final_manual_review     = as.integer(n_flags > 0 | wl_manualreview)
    ) %>%
    arrange(desc(final_manual_review), desc(n_flags), desc(n_carriers))

write_xlsx(
    list(
        "flag_overview"           = flag_overview,
        "all_variants"            = all_variants_master,
        "whitelist_manual"        = all_ch %>% filter(manualreview),
        "high_freq_summary"       = age_association %>%
            arrange(desc(n_carriers)),
        "high_freq_vars"          = high_freq_vars,
        "freeze_enriched_vars"    = freeze_removed_var %>% arrange(variant_id),
        "not_age_association_vars" = age_removed_var %>% arrange(variant_id),
        "germline_summary" = germline_summary %>%
            arrange(desc(prop_fail_binom)),
        "germline_vars"          = germline_removed_vars %>% arrange(variant_id),
        "cluster_summary"         = cluster_summary,
        "cluster_rep_vars"        = cluster_rep,
        "cluster_not_rep_vars"    = cluster_all %>%
            filter(variant_id %in% cluster_not_rep_var)
    ),
    file.path("ch", "data", "ch_wl_artifact_review.xlsx")
)

cat("Master review sheet written.\n")
cat(sprintf("Total variants:     %d\n", nrow(all_variants_master)))
# 481
cat(sprintf("Flagged for review: %d\n", sum(all_variants_master$manual_review)))
# 182
cat(sprintf("Clean variants:     %d\n", sum(all_variants_master$manual_review == 0)))
# 425

# ========================
# OVERLAP VENN DIAGRAMS
# ========================
variant_counts2 <- all_ch_clean %>%
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
    mutate(flag_high_freq = (n_carriers >= FREQ_THRESHOLD))

var_f2 <- sort(variant_counts2 %>% filter(batch_1 != 0) %>% pull(variant_id))
var_f3 <- sort(variant_counts2 %>% filter(batch_2 != 0) %>% pull(variant_id))

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

gene_f2 <- sort(unique(all_ch_clean %>% filter(Batch == 1) %>% pull(Gene)))
gene_f3 <- sort(unique(all_ch_clean %>% filter(Batch == 2) %>% pull(Gene)))

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
# QC PLOTTING FUNCTION
# ========================
plot_batch_gene_freq <- function(df, label = "",
                                 sample_col = "Sample.ID",
                                 gene_col   = "Gene",
                                 batch_col  = "Batch",
                                 batch_labels = c("Freeze 2", "Freeze 3"),
                                 out_dir    = file.path("ch", "figures")) {

    # --- overlap summary ---
    f2_vars <- df %>% filter(.data[[batch_col]] == 1) %>% pull(variant_id) %>% unique()
    f3_vars <- df %>% filter(.data[[batch_col]] == 2) %>% pull(variant_id) %>% unique()
    cat(sprintf("[%s] Variants only in F2: %d\n",  label, sum(!f2_vars %in% f3_vars)))
    cat(sprintf("[%s] Variants only in F3: %d\n",  label, sum(!f3_vars %in% f2_vars)))
    cat(sprintf("[%s] Variants in both:   %d\n\n", label, sum(f2_vars %in% f3_vars)))

    # --- batch sizes ---
    batch_sizes <- df %>%
        group_by(.data[[batch_col]]) %>%
        summarise(n_total = n_distinct(.data[[sample_col]]), .groups = "drop")

    n_total_overall <- n_distinct(df[[sample_col]])
    subtitle_str <- sprintf("%s n = %d  |  %s n = %d  |  Overall n = %d",
                            batch_labels[1], batch_sizes$n_total[batch_sizes[[batch_col]] == 1],
                            batch_labels[2], batch_sizes$n_total[batch_sizes[[batch_col]] == 2],
                            n_total_overall)

    # --- batch carrier freq ---
    gene_freq <- df %>%
        group_by(.data[[gene_col]], .data[[batch_col]]) %>%
        summarise(n_carriers = n_distinct(.data[[sample_col]]), .groups = "drop") %>%
        left_join(batch_sizes, by = batch_col) %>%
        mutate(freq  = n_carriers / n_total,
               Batch = factor(.data[[batch_col]], labels = batch_labels))

    # --- overall carrier freq ---
    gene_freq_overall <- df %>%
        group_by(.data[[gene_col]]) %>%
        summarise(n_carriers = n_distinct(.data[[sample_col]]),
                  n_total    = n_total_overall,
                  freq       = n_carriers / n_total_overall,
                  Batch      = "Overall",
                  .groups    = "drop")

    # gene order by overall freq
    gene_order <- gene_freq_overall %>%
        arrange(desc(freq)) %>%
        pull(.data[[gene_col]])

    gene_freq <- bind_rows(gene_freq, gene_freq_overall) %>%
        mutate(Gene  = factor(.data[[gene_col]], levels = gene_order),
               Batch = factor(Batch, levels = c(batch_labels, "Overall")))

    batch_colors <- c(RColorBrewer::brewer.pal(3, "Set2")[1:2], "grey30") %>%
        setNames(c(batch_labels, "Overall"))

    # --- dot plot (horizontal) ---
    p_dot <- ggplot(gene_freq, aes(x = Gene, y = freq, color = Batch)) +
        geom_line(aes(group = Gene), color = "grey80", linewidth = 0.4) +
        geom_point(size = 2.5, alpha = 0.9) +
        scale_color_manual(values = batch_colors, name = NULL) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),
                           breaks = scales::pretty_breaks(n = 6),
                           expand = expansion(mult = c(0.02, 0.1))) +
        labs(title    = sprintf("Carrier frequency per gene by batch%s",
                                if (label != "") sprintf(" - %s", label) else ""),
             subtitle = subtitle_str,
             x = NULL, y = "Carrier frequency") +
        theme_minimal(base_size = 11) +
        theme(plot.title       = element_text(face = "bold", hjust = 0.5),
              plot.subtitle    = element_text(hjust = 0.5, size = 9, color = "grey40"),
              axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
              legend.position  = "bottom",
              panel.grid.minor = element_blank(),
              panel.grid.major.x = element_blank())

    slug <- if (label != "") paste0("_", gsub("\\s+", "_", tolower(label))) else ""
    ggsave(file.path(out_dir, sprintf("qc_carrier_freq_dot%s.pdf", slug)),
           p_dot, width = 14, height = 6)

    invisible(list(dot = p_dot, freq_table = gene_freq))
}

plot_batch_gene_freq(all_ch_clean, label = "post artifact")

controls <- cov %>% filter(BRCA12_Case == 0)
dim(controls)
# 2327

control_clean_vars <- all_ch_clean %>% filter(Sample.ID %in% controls$person_id)
plot_batch_gene_freq(control_clean_vars, label = "post artifact controls")




