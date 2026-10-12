# ============================================================
# PRS MERGING
# ============================================================
library(here)
source("R/config.R")

# ========================
# LOAD DATA & CONVERT
# ========================
sso <- read_tsv(here("tgct", "ss", "pennsamples_masterslist_062226.txt"))
# 1350

# flags
flags <- fread(here("PMBB", "4.0", "rgcname_pmbbid_metadata_flags_freeze4.csv"), header = TRUE)
imputed_ids <- read.csv(here("PMBB", "4.0", "PMBB-Release-2026-4.0_genetic_imputed.sample_list.txt"), header = FALSE)$V1

# 2026 version
prs_26e <- read.table(here("tgct", "data", "output", "gwas2026_effect_pos.txt"), header = TRUE)
prs_26z <- read.table(here("tgct", "data", "output", "gwas2026_z_pos.txt"), header = TRUE)

# 2021 version
prs_21e <- read.table(here("tgct", "data", "output", "gwas2021_effect.txt"), header = TRUE)
prs_21z <- read.table(here("tgct", "data", "output", "gwas2021_z.txt"), header = TRUE)

# ========================
# COMBINE ALL PRS INTO ONE PRS
# ========================
prs <- prs_26e |>
    select(-IID, -PHENO, -CNT, -CNT2) |>
    rename(PRS26_Effect = SCORESUM) |>
    left_join(prs_26z |> select(FID, SCORESUM) |> rename(PRS26_Z = SCORESUM), by = "FID") |>
    left_join(prs_21e |> select(FID, SCORESUM) |> rename(PRS21_Effect = SCORESUM), by = "FID") |>
    left_join(prs_21z |> select(FID, SCORESUM) |> rename(PRS21_Z = SCORESUM), by = "FID")

write_xlsx(prs, here("tgct", "ss", "20260820_tgct_prs_case_control.xlsx"))

# ========================
# CASE MASTERLIST FILTERING
# ========================
# CHECK IF THERE ARE MULTIPLE IDS WITHIN SAME ROW
# choose by priority and make a GenoSource2 new column
ss <- sso %>%
    left_join(flags, by = c("PMBBID2" = "PMBB_ID")) %>%
    mutate(
        # row-wise: does THIS row individually have multiple id columns filled in?
        GenoSource2 = case_when(
            !is.na(PMBBID2) ~ "PMBB", # checked first
            !is.na(GencoveID) ~ "Gencove", # checked next, so on if PMBB was false
            !is.na(ReplicationID) ~ "Replication",
            !is.na(DiscoveryID) ~ "Discovery"
        ),
        Final_ID = coalesce(PMBBID2, GencoveID, ReplicationID, DiscoveryID)
    )
dim(ss)
# 1350

sum(is.na(ss$Geno.Avail))
# 110
# only ones w/o Geno.Avail column are the TGCT Not KT ones..
# so in the original, there are so with TRUE, NA, and FALSE for Geno.Avail
# for NA - these are TGCT Not KT, but they're in PMBB
# so join them with PMBB IDs and flag them as available
not_kt <- ss |> filter(KTID == "TGCT_Not_KT")
sum(is.na(not_kt$PMBBID2))
# 0

dim(not_kt)
sum(not_kt$PMBBID2 %in% imputed_ids)
# 109 aren't in imputed data

# View(not_kt[!(not_kt$Final_ID %in% imputed_ids), ])
# this guy not in it - PMBB4565870091478
# just remove him and mark the rest as fine

ss <- ss |> filter(is.na(PMBBID2) | PMBBID2 != "PMBB4565870091478")
dim(ss)
# 1349

table(ss$Geno.Avail)
# FALSE  TRUE
# 216  1024

# so who are the 216 FALSE?
missing_geno <- ss |> filter(Geno.Avail == FALSE)

prs_ids <- prs_26e$IID
have_geno <- missing_geno |> filter(Final_ID %in% prs_ids)
dim(have_geno)
head(have_geno)
# 1

# put that one back
ss2 <- ss |> filter(Geno.Avail == TRUE) |>
    rbind(have_geno)
ss2$Geno.Avail <- TRUE
dim(ss2)
# 1025

# get the right IDs for writing ...
# for replication + discovery, fix their ID formats
ss2 <- ss2 %>%
  mutate(Final_ID2 = if_else(GenoSource2 %in% c("Replication"),
                             paste0(Final_ID, "_", Final_ID), Final_ID)) |>
    mutate(
        EMPI = as.character(EMPI),
        MRN = as.character(MRN),
        DOB = as.Date(DOB, format = "%m/%d/%y")
    )
dim(ss2)
# 1025
write_xlsx(ss2, here("tgct", "ss", "20260817_tgct_master.xlsx"))
# hand over to wenting

# View(ss2[is.na(ss2$Final_ID2),])
# PMBB_KT0039901_KT0039901 - manually add this flagged
# manually edit this one entry, read back in
manual <- read_excel(here("tgct", "ss", "20260817_tgct_master_edit.xlsx")) |>
    filter(KTID == "KT0039901")
ss2 <- ss2 |> filter(KTID != "KT0039901") |> rbind(manual)
dim(ss2)
# 1025

# ========================
# DEDUPLICATE FROM WENTING's SHEET
# ========================
# 20 of those are duplicates according to Wenting
# 54 don't have any match in redcap
# this is an old sheet I sent over, so it has more rows than ss2
ssw <- read_excel(here("tgct", "ss", "20260730_tgct_master_w_import_info_jenny.xlsx"), sheet = "data", na = "NA") |>
    select(-CNT, -IID, -PHENO, -CNT2, -SCORESUM) |>
    mutate(
        EMPI = as.character(EMPI),
        MRN = as.character(MRN),
        DOB = as.Date(DOB, format = "%m/%d/%y")
    )
dim(ssw)
# 1136
sum(is.na(ssw$DOB))
# 111
# View(ssw[is.na(ssw$DOB),])
# these are the TGCT_Not_KT people...and this guy KT1004401

ssw1 <- ssw |> filter(Final_ID2 %in% ss2$Final_ID2)
dim(ssw1)
# 1025

dup <- read_excel(here("tgct", "ss", "20260730_tgct_master_w_import_info_jenny.xlsx"), sheet = "duplicates", na = "NA") |>
    select(-CNT, -IID, -PHENO, -CNT2, -SCORESUM) |>
    mutate(
        EMPI = as.character(EMPI),
        MRN = as.character(MRN),
        DOB = as.Date(DOB, format = "%m/%d/%y")
    )
dim(dup)
# 40 entries - 20 people are duplicated

keep <- dup |>
    mutate(
        GenoSource3 = case_when(
            !is.na(PMBBID2)      ~ "PMBB",        # checked first
            !is.na(GencoveID)    ~ "Gencove",      # checked next, so on if PMBB was false
            !is.na(ReplicationID) ~ "Replication",
            !is.na(DiscoveryID)  ~ "Discovery",
            TRUE                 ~ NA_character_
        ),
        Final_ID = coalesce(PMBBID2, GencoveID, ReplicationID, DiscoveryID),
        priority_rank = match(GenoSource3, c("PMBB", "Gencove", "Replication", "Discovery"))
    ) |>
    group_by(REDCap_record_id) |> #  within-group priority comparison happens here
    mutate(
        Representative = !is.na(priority_rank) &
            rank(priority_rank, ties.method = "first", na.last = "keep") == 1
    ) |>
    mutate(
        Final_ID3 = Final_ID[Representative][1],
        GenoSource3 = GenoSource3[Representative][1]
    ) |>
    ungroup() |>
    select(-Representative, -priority_rank)

merged <- keep |>
    merge_duplicates("REDCap_record_id")
dim(merged)
# 20 merged

# remove the TGCT_Not_KT part of the second half the semicolon away from merged$KTID
merged <- merged |>
    mutate(
        KTID = sapply(strsplit(KTID, ";"), function(x) {
            x <- trimws(x)
            x <- x[x != "TGCT_Not_KT"]
            if (length(x) == 0) "TGCT_Not_KT" else paste(x, collapse = ";")
        })
    )

# these are all GenoAvailable now
merged$GenoAvailable <- TRUE

# remove the 40 duplicate IDs from the ss and add them back
dup_remove_ids <- keep$Final_ID2

# only 18 in the duplicate remove IDs ...
sum(ss2$Final_ID2 %in% dup_remove_ids)
# 14
dup_remove_ids[!(dup_remove_ids %in% ss2$Final_ID2)]
# [1] "KT69201_KT69201"                   "PMBB7849981120953"
# [3] "TSTC_KT31501_DB1_TSTC_KT31501_DB1" "PMBB5949308134471"
# [5] "PMBB4669070078607"                 "PMBB2755815558865"
# [7] "PMBB1233436025782"                 "PMBB9093337904007"
# [9] "PMBB1510264236523"                 "PMBB2756885606141"
# [11] "PMBB4112293723018"                 "PMBB7555736208392"
# [13] "PMBB9622336032181"                 "PMBB3640218039050"
# [15] "PMBB6103969212975"                 "KT62401_KT62401"
# [17] "PMBB5307202982640"                 "KT75301_KT75301"
# [19] "PMBB2247235279001"                 "PMBB5807319271199"
# [21] "PMBB6487687571593"                 "PMBB7157941302055"
# [23] "KT64601_KT64601"                   "PMBB2192862851834"
# [25] "KT25401_C06_KT25401_C06"           "PMBB9207827289760"

# I think there other 22 I already removed from before

dim(ss2)
# 1025
ss3 <- ss2 |> filter(!(Final_ID2 %in% dup_remove_ids))
dim(ss3)
# 1011

# now have to merge flags back in
merged2 <- merged |> left_join(flags |> select(-Freeze, -PMBB_consent), by = c("PMBBID2" = "PMBB_ID")) |>
    mutate(
        EMPI = as.character(EMPI),
        MRN = as.character(MRN),
        DOB = as.Date(DOB, format = "%m/%d/%y")
    )

# prepare for merging
setdiff(names(ss3), names(merged2))
setdiff(names(merged2), names(ss3))

# rename GenoAvailable to Geno.Avail in merged3 and remove RowNum column
merged3 <- merged2 |> rename(Geno.Avail = GenoAvailable) |> select(-RowNum, -REDCap_record_id)

# add a Final_ID3 col overall that doesn't have the duplicated sample IDs
ss3$Final_ID3 <- ss3$Final_ID2
# add a GenoSource3 col overall that doesn't have the duplicated sources
ss3$GenoSource3 <- ss3$GenoSource2

setdiff(names(ss3), names(merged3))
setdiff(names(merged3), names(ss3))
# we good

dim(ss3)
# 1011
dim(merged2)
# 20

ss4 <- bind_rows(ss3, merged3)
dim(ss4)
# 1031

# join PRS -- reuse the already-combined `prs` table (built above) instead of
# re-joining each of the 4 raw PRS files again here, so there's only one
# place that maps SCORESUM -> PRS26_Effect/PRS26_Z/PRS21_Effect/PRS21_Z
# discovery actually doesnt' ahve the duplicated IDs...
ss4 <- ss4 |> left_join(prs, by = c("Final_ID3" = "FID"))

sum(is.na(ss4$PRS26_Z))
# 24 without
# View(ss4[is.na(ss4$PRS26_Z),])
# these are mostly Gencove and discovery
# examine their files, but still can't find them ...

# filter out those w/o PRS
table(ss4$PMBB_consent)
# current        ended         null null;current
# 273            5          553            1

ss5 <- ss4 |> filter(!is.na(PRS26_Z)) |>
    filter(PMBB_consent != "ended" | is.na(PMBB_consent))
dim(ss5)
# 1003

# ========================
# FINAL STATS
# ========================
table(ss5$GenoSource3)
# Discovery     Gencove        PMBB Replication
# 64          50         827          62

table(ss5$Freeze)
# Freeze 1          Freeze 2          Freeze 3 Freeze 3;Freeze 2          Freeze 4
# 6                69               749                 1                 2

# clean up
ss6 <- ss5 |> select(-ReplicationID, -GencoveID, -RegeneronID, -PMBBID1, -PMBBID2, -Geno.Avail,
                     -GenoSource, -AncestryContinental, -rgc_sample_name, -ID1, -ID2, -ID3, -ID4, -ID5,
                     -LRC, -CREP, -Melanoma, -Testicular_CA, -Para_Pheo)
dim(ss6)
write_xlsx(ss6, here("tgct", "ss", "20260820_tgct_master_cleaned.xlsx"))
cases <- ss6

# ========================
# CHECK IF DUPLICATE SCORES ARE SAME
# ========================
dup <- dup |>  mutate(Final_ID2 = if_else(GenoSource2 %in% c("Replication"),
                                    paste0(Final_ID, "_", Final_ID), Final_ID))
check <- dup |> left_join(prs, by = c("Final_ID2" = "FID"))
write_xlsx(check, file.path("tgct", "data", "20260831_prs_dups.xlsx"))

# CHECK IF THERE ARE MULTIPLE IDS WITHIN SAME ROW
ss <- sso %>%
  left_join(flags, by = c("PMBBID2" = "PMBB_ID")) %>%
  mutate(
    # row-wise: does THIS row individually have multiple id columns filled in?
    GenoSource2 = case_when(
      !is.na(PMBBID2) ~ "PMBB", # checked first
      !is.na(GencoveID) ~ "Gencove", # checked next, so on if PMBB was false
      !is.na(ReplicationID) ~ "Replication",
      !is.na(DiscoveryID) ~ "Discovery"
    ),
    Final_ID = coalesce(PMBBID2, GencoveID, ReplicationID, DiscoveryID)
  )

# check if any row has more than one ID in the column....
ss <- sso %>%
  left_join(flags, by = c("PMBBID2" = "PMBB_ID")) %>%
  mutate(
    .row_id = row_number(),
    GenoSource2 = case_when(
      !is.na(PMBBID2) ~ "PMBB",
      !is.na(GencoveID) ~ "Gencove",
      !is.na(ReplicationID) ~ "Replication",
      !is.na(DiscoveryID) ~ "Discovery"
    ),
    Final_ID = coalesce(PMBBID2, GencoveID, ReplicationID, DiscoveryID),
    source_col = recode(GenoSource2,
                        PMBB = "PMBBID2", Gencove = "GencoveID",
                        Replication = "ReplicationID", Discovery = "DiscoveryID"
    )
  )

id_cols <- c("PMBBID2", "GencoveID", "ReplicationID", "DiscoveryID")

id_long <- ss %>%
  select(.row_id, all_of(id_cols)) %>%
  pivot_longer(all_of(id_cols), names_to = "id_col", values_to = "id_value") %>%
  filter(!is.na(id_value))

dupe_rows <- ss %>%
  select(.row_id, source_col, Final_ID) %>%
  inner_join(id_long, by = c("Final_ID" = "id_value"), suffix = c("", "_other")) %>%
  filter(.row_id != .row_id_other, source_col != id_col) %>%
  distinct(.row_id) %>%
  pull(.row_id)

ss <- ss %>%
  mutate(ID_in_other_row_diff_col = .row_id %in% dupe_rows)

# ========================
# PRS SUMMARY STATISTICS
# ========================
prs_cols <- c("PRS26_Effect", "PRS26_Z", "PRS21_Effect", "PRS21_Z")

prs_stats <- function(x) {
  tibble(
    n = sum(!is.na(x)),
    n_missing = sum(is.na(x)),
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    IQR = IQR(x, na.rm = TRUE),
    Q1 = quantile(x, 0.25, na.rm = TRUE),
    Q3 = quantile(x, 0.75, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
}

prs_summary <- map_dfr(prs_cols, function(col) {
  by_source <- cases |>
    group_by(GenoSource3) |>
    group_modify(~ prs_stats(.x[[col]])) |>
    rename(Source = GenoSource3)

  overall <- prs_stats(cases[[col]]) |> mutate(Source = "Overall", .before = 1)

  bind_rows(by_source, overall) |> mutate(PRS = col, .before = 1)
})

print(prs_summary)
write_xlsx(prs_summary, here("tgct", "data", "20260910_tgct_prs_summary_stats.xlsx"))

# ========================
# PRS DISTRIBUTIONS BY SOURCE
# ========================
table(cases$GenoSource3)
# Discovery     Gencove        PMBB Replication 
# 64          50         827          62 

dist_long <- cases |>
    select(GenoSource3, all_of(prs_cols)) |>
    pivot_longer(all_of(prs_cols), names_to = "PRS", values_to = "Score")

stat_test <- dist_long |>
  group_by(PRS) |>
  wilcox_test(Score ~ GenoSource3) |>          
  adjust_pvalue(method = "BH") |>               
  add_significance("p.adj") |>
  add_xy_position(x = "GenoSource3", fun = "max", step.increase = 0.12)

p_dist <- ggplot(dist_long, aes(x = GenoSource3, y = Score, fill = GenoSource3)) +
  geom_boxplot(outlier.shape = 21, outlier.alpha = 0.5) +
  facet_wrap(~PRS, scales = "free_y") +
  stat_pvalue_manual(stat_test, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE) +
  labs(x = "Genotyping Source", y = "PRS", fill = "Source") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_dist)
ggsave(here("tgct", "figures", "20260908_tgct_prs_distribution_by_source_sig.png"), plot = p_dist, width = 10, height = 8)

write_xlsx(stat_test, here("tgct", "data", "20260908_tgct_prs_pairwise_by_source.xlsx"))

# print out significant ones.. there are none
stat_test |>
  mutate(across(c(statistic, p, p.adj), \(x) round(x, 4))) |>
  select(PRS, group1, group2, n1, n2, p, p.adj, p.adj.signif) |>
  arrange(PRS, p.adj) |>
  print(n = Inf)

# ========================
# PAIRWISE SCATTERPLOTS
# ========================
library(patchwork)

prs_pairs <- combn(prs_cols, 2, simplify = FALSE)

for (pair in prs_pairs) {
    p <- ggplot(cases, aes(x = .data[[pair[1]]], y = .data[[pair[2]]], color = GenoSource3)) +
        geom_point(alpha = 0.6, size = 1.5) +
        stat_cor(aes(color = NULL), method = "pearson", label.sep = "\n") +
        labs(x = pair[1], y = pair[2], color = "Source") +
        theme_bw()
    print(p)
    fname <- paste0("20260901_tgct_prs_scatter_", pair[1], "_vs_", pair[2], ".png")
    ggsave(here("tgct", "figures", fname), plot = p, width = 6, height = 5)
}

prs_pairs <- list(
  c("PRS26_Effect", "PRS21_Effect"),
  c("PRS26_Effect", "PRS26_Z"),
  c("PRS21_Effect", "PRS21_Z"),
  c("PRS26_Z",      "PRS21_Z")
)

plots <- lapply(prs_pairs, function(pair) {
  ggplot(cases, aes(x = .data[[pair[1]]], y = .data[[pair[2]]], color = GenoSource3)) +
    geom_point(alpha = 0.6, size = 1.5) +
    stat_cor(aes(color = NULL), method = "pearson", label.sep = "\n") +
    labs(x = pair[1], y = pair[2], color = "Source") +
    theme_bw()
})

p_all <- wrap_plots(plots, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(p_all)

ggsave(
  here("tgct", "figures", "20260901_tgct_prs_scatter_all_pairs.png"),
  plot = p_all, width = 9, height = 9, dpi = 300
)

# ========================
# PRS TERTILES 
# ========================
effect26_breaks   <- quantile(cases$PRS26_Effect, probs = c(1 / 3, 2 / 3), na.rm = TRUE)
z26_breaks        <- quantile(cases$PRS26_Z,      probs = c(1 / 3, 2 / 3), na.rm = TRUE)
effect21_breaks <- quantile(cases$PRS21_Effect, probs = c(1 / 3, 2 / 3), na.rm = TRUE)
z21_breaks      <- quantile(cases$PRS21_Z,      probs = c(1 / 3, 2 / 3), na.rm = TRUE)

# label each sample's tertile per PRS type
cases <- cases |> mutate(
  PRS26_Effect_tertile = cut(PRS26_Effect, breaks = c(-Inf, effect26_breaks, Inf), labels = c("Low", "Mid", "High")),
  PRS26_Z_tertile       = cut(PRS26_Z,      breaks = c(-Inf, z26_breaks, Inf),      labels = c("Low", "Mid", "High")),
  PRS21_Effect_tertile = cut(PRS21_Effect, breaks = c(-Inf, effect21_breaks, Inf), labels = c("Low", "Mid", "High")),
  PRS21_Z_tertile       = cut(PRS21_Z,      breaks = c(-Inf, z21_breaks, Inf),      labels = c("Low", "Mid", "High"))
)

table(cases$PRS26_Effect_tertile, useNA = "ifany")
# table(cases$PRS26_Z_tertile, useNA = "ifany")
# table(cases$PRS21_Effect_tertile, useNA = "ifany")
# table(cases$PRS21_Z_tertile, useNA = "ifany")
# Low  Mid High
# 335  334  334

# ========================
# CHECK DUPLICATES
# ========================
check <- check |>
  filter(!is.na(PRS26_Effect)) |>
  group_by(REDCap_record_id) |>
  filter(n() >= 2) |>
  ungroup()
check |> group_by(REDCap_record_id) |> summarize(sources = paste(sort(GenoSource2), collapse = " vs ")) |> count(sources)

# # A tibble: 4 × 2
# sources                 n
# <chr>               <int>
# 1 Discovery vs PMBB       3
# 2 Gencove vs PMBB         7
# 3 PMBB vs PMBB            1
# 4 PMBB vs Replication     4

wide <- check |>
  distinct(REDCap_record_id, GenoSource2, .keep_all = TRUE) |>
  pivot_wider(id_cols = REDCap_record_id, names_from = GenoSource2, values_from = PRS26_Effect)

paired <- check |>
  group_by(REDCap_record_id) |>
  mutate(n_pmbb = sum(GenoSource2 == "PMBB")) |>
  filter(n_pmbb == 1) |>  
  summarize(
    pmbb_score   = PRS26_Effect[GenoSource2 == "PMBB"],
    other_score  = PRS26_Effect[GenoSource2 != "PMBB"],
    other_source = GenoSource2[GenoSource2 != "PMBB"],
    .groups = "drop"
  )

# paired_long <- paired |>
#   mutate(person = row_number()) |>
#   pivot_longer(c(pmbb_score, other_score), names_to = "source_type", values_to = "score") |>
#   mutate(source_type = if_else(source_type == "pmbb_score", "PMBB", other_source))
# 
# p_slope <- ggplot(paired_long, aes(x = source_type, y = score, group = person)) +
#   geom_line(alpha = 0.5) +
#   geom_point(aes(color = source_type), size = 2) +
#   facet_wrap(~other_source, scales = "free_x") +
#   labs(x = NULL, y = "PRS26_Effect") +
#   theme_bw()
# 
# ggsave(here("tgct", "figures", "20260908_tgct_prs_paired_slope_by_source_26effect.png"),
#        plot = p_slope, width = 8, height = 5)

# what if we graphed the tertile breaks as a line to see if any cross the line?
paired <- paired |>
  mutate(
    pmbb_tertile    = cut(pmbb_score,  breaks = c(-Inf, effect26_breaks, Inf), labels = c("Low", "Mid", "High")),
    other_tertile   = cut(other_score, breaks = c(-Inf, effect26_breaks, Inf), labels = c("Low", "Mid", "High")),
    crosses_tertile = pmbb_tertile != other_tertile,
    person          = row_number()
  )

cat(sum(paired$crosses_tertile), "of", nrow(paired), "duplicate pairs land in different tertiles\n")
paired |> select(REDCap_record_id, other_source, pmbb_score, pmbb_tertile, other_score, other_tertile, crosses_tertile) |>
  print(n = Inf)

paired_long <- paired |>
  pivot_longer(c(pmbb_score, other_score), names_to = "source_type", values_to = "score") |>
  mutate(source_type = if_else(source_type == "pmbb_score", "PMBB", other_source))

# p_slope_tertile <- ggplot(paired_long, aes(x = source_type, y = score, group = person)) +
#   geom_hline(yintercept = effect26_breaks, linetype = "dashed", color = "grey50") +
#   geom_line(aes(color = crosses_tertile), alpha = 0.7) +
#   geom_point(aes(color = crosses_tertile), size = 2) +
#   facet_wrap(~other_source, scales = "free_x") +
#   scale_color_manual(values = c(`FALSE` = "grey50", `TRUE` = "red"), name = "Crosses tertile?") +
#   labs(x = NULL, y = "PRS26_Effect") +
#   theme_bw()
# 
# ggsave(here("tgct", "figures", "20260908_tgct_prs_paired_slope_tertile_crossing_26effect.png"),
#        plot = p_slope_tertile, width = 8, height = 5)

# ========================
# PAIRWISE COMPARISON
# ========================
prs_cols <- c("PRS26_Effect", "PRS26_Z", "PRS21_Effect", "PRS21_Z")

# per-column tertile breaks computed earlier in the PRS TERTILES section
tertile_breaks <- list(
  PRS26_Effect = effect26_breaks,
  PRS26_Z      = z26_breaks,
  PRS21_Effect = effect21_breaks,
  PRS21_Z      = z21_breaks
)

for (col in prs_cols) {

  check_i <- dup |>
    left_join(prs, by = c("Final_ID2" = "FID")) |>
    filter(!is.na(.data[[col]])) |>
    group_by(REDCap_record_id) |>
    filter(n() >= 2) |>
    ungroup()

  paired_i <- check_i |>
    group_by(REDCap_record_id) |>
    mutate(n_pmbb = sum(GenoSource2 == "PMBB")) |>
    filter(n_pmbb == 1) |>
    summarize(
      pmbb_score   = .data[[col]][GenoSource2 == "PMBB"],
      pmbb_id      = Final_ID2[GenoSource2 == "PMBB"],
      other_score  = .data[[col]][GenoSource2 != "PMBB"],
      other_id     = Final_ID2[GenoSource2 != "PMBB"],
      other_source = GenoSource2[GenoSource2 != "PMBB"],
      .groups = "drop"
    )

  if (nrow(paired_i) == 0) next

  breaks_i <- tertile_breaks[[col]]

  paired_i <- paired_i |>
    mutate(
      pmbb_tertile    = cut(pmbb_score,  breaks = c(-Inf, breaks_i, Inf), labels = c("Low", "Mid", "High")),
      other_tertile   = cut(other_score, breaks = c(-Inf, breaks_i, Inf), labels = c("Low", "Mid", "High")),
      crosses_tertile = pmbb_tertile != other_tertile,
      person          = row_number()
    )

  cat(col, ":", sum(paired_i$crosses_tertile), "of", nrow(paired_i), "duplicate pairs cross tertiles\n")

  paired_long_i <- bind_rows(
    paired_i |> transmute(person, other_source, crosses_tertile, source_type = "PMBB",       score = pmbb_score,  id = pmbb_id),
    paired_i |> transmute(person, other_source, crosses_tertile, source_type = other_source, score = other_score, id = other_id)
  )

  p <- ggplot(paired_long_i, aes(x = source_type, y = score, group = person)) +
    geom_hline(yintercept = breaks_i, linetype = "dashed", color = "grey50") +
    geom_line(alpha = 0.5) +
    geom_point(aes(color = source_type, shape = crosses_tertile), size = 2) +
    geom_text_repel(aes(label = id), size = 2.3, show.legend = FALSE,
                    max.overlaps = Inf, segment.size = 0.2) +
    facet_wrap(~other_source, scales = "free_x") +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 8), name = "Crosses tertile?") +
    labs(x = NULL, y = col, title = paste("Paired PRS by GenoSource:", col)) +
    theme_bw()

  ggsave(here("tgct", "figures", paste0("20260908_tgct_prs_paired_slope_", col, ".png")),
         plot = p, width = 8, height = 6)
}


# ========================
# STANDARDIZE LIKE EMILIO
# ========================
plot_data <- cases |>
  mutate(PRS26_Z_std = as.numeric(scale(PRS26_Z))) |>   # mean 0, SD 1 across your own cohort
  filter(!is.na(PRS26_Z_tertile))

facet_means <- plot_data |>
  group_by(PRS26_Z_tertile, GenoSource3) |>
  summarize(mean_prs = mean(PRS26_Z_std), .groups = "drop")

p_dist_tecac <- ggplot(plot_data, aes(x = PRS26_Z_std)) +
  geom_histogram(aes(y = after_stat(density)), fill = "pink", color = "black", alpha = 0.6, bins = 30) +
  geom_density(color = "black", linewidth = 0.6) +
  geom_vline(data = facet_means, aes(xintercept = mean_prs), color = "blue", linetype = "dashed") +
  facet_grid(rows = vars(PRS26_Z_tertile, GenoSource3)) +   # tertile outer, source inner -- matches Emilio's layout
  labs(x = "PRS26_Z (standardized)", y = "density") +
  theme_bw()

ggsave(here("tgct", "figures", "20260908_tgct_prs_std_tertile_by_source.png"),
       plot = p_dist_tecac, width = 8, height = 12)

