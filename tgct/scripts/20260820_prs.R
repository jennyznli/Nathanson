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
prs_26e <- read.table(here("tgct", "data", "output", "gwas2026_effect.txt"), header = TRUE)
prs_26z <- read.table(here("tgct", "data", "output", "gwas2026_z.txt"), header = TRUE)

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
# but they all have KT IDS?

missing_geno_ids <- missing_geno$KTID
sum(prs_ids %in% missing_geno_ids)
# 0
# so actually probably missing ...

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
View(ss4[is.na(ss4$PRS26_Z),])
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
# COMPUTE TERTILES on CONTROLS ONLY
# ========================
# get all cases ever
# but what if some are duplicats...
# so literally take all IDs across PMBBID1, PMBBID2, GencoveID, ReplicationID, DiscoveryID and
# make it a unique vector
# but remember replication has to be like the existing one ID_ID with the underscore and doubled.
# pulled from sso (the untouched 1350-row masterlist) rather than ss6 -- ss6
# is the final deduped/filtered case cohort, so anyone dropped along the way
# (missing genotype, removed dup, etc.) would wrongly fall through to
# "control" if we sourced case IDs from ss6 instead
case_ids <- sso |>
    mutate(ReplicationID = if_else(!is.na(ReplicationID), paste0(ReplicationID, "_", ReplicationID), ReplicationID)) |>
    select(PMBBID1, PMBBID2, GencoveID, ReplicationID, DiscoveryID) |>
    unlist(use.names = FALSE) |>
    na.omit() |>
    unique()
length(case_ids)
# 2751

# mark case/controls
prs2 <- prs |> mutate(Case = if_else(FID %in% case_ids, 1, 0))
table(prs2$Case)
# 0     1
# 84282  1799

# computed on ss4 (everyone with a valid PRS score) rather than on the later
# complete-case survival cohort -- PRS is fixed regardless of who has full
# clinical follow-up, so tertile cutpoints shouldn't depend on which subset
# happens to have complete histology/relapse/date fields. fixed cutpoints
# also mean "Low/Mid/High" means the same PRS value in every plot and every
# subgroup below, instead of being recomputed (and redefined) per call
controls <- prs2 |> filter(Case == 0)

effect_breaks   <- quantile(controls$PRS26_Effect, probs = c(1 / 3, 2 / 3), na.rm = TRUE)
z_breaks        <- quantile(controls$PRS26_Z,      probs = c(1 / 3, 2 / 3), na.rm = TRUE)
effect21_breaks <- quantile(controls$PRS21_Effect, probs = c(1 / 3, 2 / 3), na.rm = TRUE)
z21_breaks      <- quantile(controls$PRS21_Z,      probs = c(1 / 3, 2 / 3), na.rm = TRUE)

effect_breaks
z_breaks
effect21_breaks
z21_breaks

# label each sample's tertile per PRS type -- carries through ss5/survival/
# sem_adj/sem_surv automatically since it's just a column on ss4
ss7 <- ss6 |> mutate(
    PRS26_Effect_tertile = cut(PRS26_Effect, breaks = c(-Inf, effect_breaks, Inf), labels = c("Low", "Mid", "High")),
    PRS26_Z_tertile       = cut(PRS26_Z,      breaks = c(-Inf, z_breaks, Inf),      labels = c("Low", "Mid", "High")),
    PRS21_Effect_tertile = cut(PRS21_Effect, breaks = c(-Inf, effect21_breaks, Inf), labels = c("Low", "Mid", "High")),
    PRS21_Z_tertile       = cut(PRS21_Z,      breaks = c(-Inf, z21_breaks, Inf),      labels = c("Low", "Mid", "High"))
)

table(ss7$PRS26_Effect_tertile, useNA = "ifany")
table(ss7$PRS26_Z_tertile, useNA = "ifany")
table(ss7$PRS21_Effect_tertile, useNA = "ifany")
table(ss7$PRS21_Z_tertile, useNA = "ifany")

# ========================
# PMBB RESTRICTED - COMPUTE TERTILES on CONTROLS ONLY
# ========================
pmbb_ids <- unique(flags$PMBB_ID)
pmbb_controls <- prs2 |> filter(Case == 0, FID %in% pmbb_ids)
dim(pmbb_controls)
# 69488

effect_breaks   <- quantile(pmbb_controls$PRS26_Effect, probs = c(1 / 3, 2 / 3), na.rm = TRUE)
z_breaks        <- quantile(pmbb_controls$PRS26_Z,      probs = c(1 / 3, 2 / 3), na.rm = TRUE)
effect21_breaks <- quantile(pmbb_controls$PRS21_Effect, probs = c(1 / 3, 2 / 3), na.rm = TRUE)
z21_breaks      <- quantile(pmbb_controls$PRS21_Z,      probs = c(1 / 3, 2 / 3), na.rm = TRUE)

effect_breaks
z_breaks
effect21_breaks
z21_breaks
# this is insanely high...not quite right??

pmbb_cases <- ss6 |> filter(Final_ID3 %in% pmbb_ids)
dim(pmbb_cases)
# 827

pmbb_cases <- pmbb_cases |> mutate(
        PRS26_Effect_tertile = cut(PRS26_Effect, breaks = c(-Inf, effect_breaks, Inf), labels = c("Low", "Mid", "High")),
        PRS26_Z_tertile       = cut(PRS26_Z,      breaks = c(-Inf, z_breaks, Inf),      labels = c("Low", "Mid", "High")),
        PRS21_Effect_tertile = cut(PRS21_Effect, breaks = c(-Inf, effect21_breaks, Inf), labels = c("Low", "Mid", "High")),
        PRS21_Z_tertile       = cut(PRS21_Z,      breaks = c(-Inf, z21_breaks, Inf),      labels = c("Low", "Mid", "High"))
    )

table(pmbb_cases$PRS26_Effect_tertile, useNA = "ifany")
# Low  Mid High
# 53  172  602

table(pmbb_cases$PRS26_Z_tertile, useNA = "ifany")
# Low  Mid High
# 52  177  598

table(pmbb_cases$PRS21_Effect_tertile, useNA = "ifany")
# Low  Mid High
# 68  228  531

table(pmbb_cases$PRS21_Z_tertile, useNA = "ifany")
# Low  Mid High
# 75  226  526

write_xlsx(pmbb_cases, here("tgct", "ss", "20260820_tgct_pmbb_cases.xlsx"))





