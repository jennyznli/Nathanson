# ============================================================
# CNA ANALYSIS - BATCH PROCESSING
# ============================================================
source("config.R")

# ========================
# LOAD DATA & CONVERT
# ========================
meta  <- read_xlsx(here("tgct", "ss", "meta-hits-062526.xlsx"))
prs <- read_xlsx(here("tgct", "data", "tgct_pooled.xlsx"))
snps <- fread(here("tgct", "data", "emilio-scores.txt"))

# ssw <- read_excel(here("tgct", "ss", "20260730_tgct_master_w_import_info.xlsx"))
dup <- read_excel(here("tgct", "ss", "20260730_tgct_master_w_import_info.xlsx"), sheet = "duplicates", na = "NA") |>
    mutate(
        EMPI = as.character(EMPI),
        MRN = as.character(MRN),
        DOB = as.Date(DOB, format = "%m/%d/%y")
    )
# ssw <- read_excel(here("tgct", "ss", "20260730_tgct_master_w_import_info.xlsx"), sheet = "duplicates")

sso <- read_tsv(here("tgct", "ss", "pennsamples_masterslist_062226.txt"))
ss_old <- read_excel(here("tgct", "ss", "20260730_tgct_master.xlsx"))

imputed_ids <- read.csv(here("PMBB", "4.0", "PMBB-Release-2026-4.0_genetic_imputed.sample_list.txt"), header = FALSE)$V1

dim(sso)
# 1350
dim(dup)
# 40 - 20 people are duplicated

flags <- fread(here("PMBB", "4.0", "rgcname_pmbbid_metadata_flags_freeze4.csv"), header = TRUE)

# ========================
# DEDUPLICATE BY MATCHING NAMES
# ========================
# okay so for the old sheet, check if there are any with multiple IDs WITHIN the same row
# choose by priority and make a GenoSourceFinal new column
ss <- sso %>%
  left_join(flags %>% select(PMBB_ID, Freeze, PMBB_consent), by = c("PMBBID2" = "PMBB_ID")) %>%
  mutate(
    # row-wise: does THIS row individually have multiple id columns filled in?
    GenoSourceFinal = case_when(
      !is.na(PMBBID2) ~ "PMBB", # checked first
      !is.na(GencoveID) ~ "Gencove", # checked next, so on if PMBB was false
      !is.na(ReplicationID) ~ "Replication",
      !is.na(DiscoveryID) ~ "Discovery"
    ),
    Final_ID = coalesce(PMBBID2, GencoveID, ReplicationID, DiscoveryID)
  )
sum(is.na(ss$Geno.Avail))
# 110 TGCT Not Kts
# but 20 of those are duplicates according to Wenting
# 54 don't have any match in redcap

# only ones w/o Geno.Avail are the TGCT Not KT ones..
# chekc if they're in the PMBB and if so mark them as available
not_kt <- ss |> filter(KTID == "TGCT_Not_KT") |>
    inner_join(flags, by = c("Final_ID" = "PMBB_ID"))
sum(not_kt$Final_ID %in% imputed_ids)
# 109 - so one w/o genotyping data
# View(not_kt[!(not_kt$Final_ID %in% imputed_ids), ])
# this guy not in it? PMBB4565870091478
# just remove him and mark the rest as fine

ss2 <- ss |> filter(Final_ID != "PMBB4565870091478")
ss2$Geno.Avail <- TRUE
dim(ss2)
# 1135

# ========================
# DEDUPLICATE
# ========================
ss <- read_excel(here("tgct", "ss", "20260730_tgct_master_w_import_info.xlsx"), sheet = "data", na = "NA") |>
    mutate(
        EMPI = as.character(EMPI),
        MRN = as.character(MRN),
        DOB = as.Date(DOB, format = "%m/%d/%y")
    )

dim(ss)
# same

keep <- dup |>
    mutate(
        GenoSourceFinal = case_when(
            !is.na(PMBBID2)      ~ "PMBB",        # checked first
            !is.na(GencoveID)    ~ "Gencove",      # checked next, so on if PMBB was false
            !is.na(ReplicationID) ~ "Replication",
            !is.na(DiscoveryID)  ~ "Discovery",
            TRUE                 ~ NA_character_
        ),
        Final_ID = coalesce(PMBBID2, GencoveID, ReplicationID, DiscoveryID),
        priority_rank = match(GenoSourceFinal, c("PMBB", "Gencove", "Replication", "Discovery"))
    ) |>
    group_by(REDCap_record_id) |>                 # <- within-group priority comparison happens here
    mutate(
        Representative = !is.na(priority_rank) &
            rank(priority_rank, ties.method = "first", na.last = "keep") == 1
    ) |>
    ungroup()

# so basically i just wanna keep everything in the representative rows
# EXCEPT the columns GencoveID, ReplicationID, DiscoveryID, PMBBID1, PMBBID2 -- those should be merged across each group
merged_ids <- keep |>
    select(REDCap_record_id, GencoveID, ReplicationID, DiscoveryID, PMBBID1, PMBBID2) |>
    merge_duplicates("REDCap_record_id")

final <- keep |>
    filter(Representative) |>
    select(-GencoveID, -ReplicationID, -DiscoveryID, -PMBBID1, -PMBBID2) |>
    left_join(merged_ids, by = "REDCap_record_id")

# remove the duplicates from the ss
dup_remove_ids <- dup$Final_ID
length(dup_remove_ids)
# 40

# merge this back?
ss2 <- ss |> filter(!(Final_ID %in% dup_remove_ids))
dim(ss2)

ss3 <- bind_rows(ss2, final)
dim(ss3)
# 1116

table(ss3$GenoSource2)
# Discovery     Gencove        PMBB Replication
# 61          72         921          62

ss4 <- ss3 |> filter(!is.na(SCORESUM))
sum(is.na(ss4$SCORESUM))
# 24 without PRS

write_xlsx(ss3, here("tgct", "ss", "20260812_tgct_master.xlsx"))

# ========================
# ADJUSTING IDS - original
# ========================
# ids_out <- ss2 %>%
#   mutate(write_id = if_else(GenoSourceFinal == "Replication", paste0(Final_ID, "_", Final_ID), Final_ID)) %>%
#   transmute(FID = write_id, IID = write_id)
#
# write_tsv(ids_out, here("tgct", "data", "tgct_prs_ids.txt"), col_names = FALSE)
# table(ss2$GenoSourceFinal)
# # Discovery     Gencove        PMBB Replication
# # 67          81         922          66
#
# # Freeze 1 Freeze 2 Freeze 3 Freeze 4
# # 25      126      768        3
#
# # write_xlsx(ss2, here("tgct", "ss", "20260710_tgct_master.xlsx"))
#
# # get the right IDs for writing ...
# ss2 <- ss2 %>%
#   mutate(Final_ID2 = if_else(GenoSourceFinal %in% c("Replication", "Discovery"), paste0(Final_ID, "_", Final_ID), Final_ID))
#
# ss3 <- ss2 %>% left_join(prs, by = c("Final_ID2" = "FID"))
#
# write_xlsx(ss3, here("tgct", "ss", "20260730_tgct_master2.xlsx"))





