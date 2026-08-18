# ============================================================
# CNA ANALYSIS - BATCH PROCESSING
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

# prs scores
prsz <- read.table(here("tgct", "data", "tgct_prs_z.profile"), header = TRUE)
prse <- read.table(here("tgct", "data", "tgct_prs_effect.profile"), header = TRUE)

# ========================
# CHECK IF THERE ARE MULTIPLE IDS WITHIN SAME ROW
# ========================
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

prs_ids <- prsz$IID
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
  mutate(Final_ID2 = if_else(GenoSource2 %in% c("Replication", "Discovery"),
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

View(ss2[is.na(ss2$Final_ID2),])
# PMBB_KT0039901_KT0039901 - manually add this flagged
# manually edit this one entry, read back in
ss2 <- read_excel(here("tgct", "ss", "20260817_tgct_master_edit.xlsx"))
dim(ss2)
# 1025

# ========================
# DEDUPLICATE FROM WENTING's SHEET
# ========================
# 20 of those are duplicates according to Wenting
# 54 don't have any match in redcap
# this is an old sheet I sent over, so it has more rows than ss2
ssw <- read_excel(here("tgct", "ss", "20260730_tgct_master_w_import_info.xlsx"), sheet = "data", na = "NA") |>
    select(-CNT, -IID, -PHENO, -CNT2, -SCORESUM) |>
    mutate(
        EMPI = as.character(EMPI),
        MRN = as.character(MRN),
        DOB = as.Date(DOB, format = "%m/%d/%y")
    )
dim(ssw)
# 1136

ssw1 <- ssw |> filter(Final_ID2 %in% ss2$Final_ID2)
dim(ssw1)
# 1025

dup <- read_excel(here("tgct", "ss", "20260730_tgct_master_w_import_info.xlsx"), sheet = "duplicates", na = "NA") |>
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
dup_remove_ids[!(dup_remove_ids %in% ss2$Final_ID2)]
# [1] "PMBB7849981120953" "PMBB5949308134471" "PMBB4669070078607" "PMBB2755815558865"
# [5] "PMBB1233436025782" "PMBB9093337904007" "PMBB1510264236523" "PMBB2756885606141"
# [9] "PMBB4112293723018" "PMBB7555736208392" "PMBB9622336032181" "PMBB3640218039050"
# [13] "PMBB6103969212975" "KT62401_KT62401"   "PMBB5307202982640" "PMBB2247235279001"
# [17] "PMBB5807319271199" "PMBB6487687571593" "PMBB7157941302055" "KT64601_KT64601"
# [21] "PMBB2192862851834" "PMBB9207827289760"

# I think there other 22 I already removed from before

dim(ss2)
# 1025
ss3 <- ss2 |> filter(!(Final_ID2 %in% dup_remove_ids))
dim(ss3)
# 1007

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
# 1007
dim(merged2)
# 20

ss4 <- bind_rows(ss3, merged3)
dim(ss4)
# 1027

# join PRS
ss4 <- ss4 |> left_join(prsz |> select(FID, SCORESUM) |> rename(PRS_Z = SCORESUM), by = c("Final_ID3" = "FID")) |>
    left_join(prse |> select(FID, SCORESUM) |> rename(PRS_Effect = SCORESUM), by = c("Final_ID3" = "FID"))
sum(is.na(ss4$PRS_Effect))
# 23 without
View(ss4[is.na(ss4$PRS_Effect),])
# these are mostly Gencove
# examine their files, but still can't find them ...

# filter out those w/o PRS
table(ss4$PMBB_consent)
# current        ended         null null;current
# 273            5          553            1

ss5 <- ss4 |> filter(!is.na(PRS_Z)) |>
    filter(PMBB_consent != "ended" | is.na(PMBB_consent))
dim(ss5)
# 1000

# ========================
# FINAL STATS
# ========================
table(ss5$GenoSource3)
# Discovery     Gencove        PMBB Replication
# 61          72         827          62

table(ss5$Freeze)
# Freeze 1          Freeze 2          Freeze 3 Freeze 3;Freeze 2          Freeze 4
# 6                69               749                 1                 2

write_xlsx(ss5, here("tgct", "ss", "20260818_tgct_master_cleaned.xlsx"))



