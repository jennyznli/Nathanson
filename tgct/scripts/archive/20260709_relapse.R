# ========================
# RELAPSE
# ========================
library(here)
source(here("R", "config.R"))
setwd(here("tgct"))

# ========================
# LOAD DATA
# ========================
# fam <- read.table(here("tgct", "data", "pmbb_prs.fam"), header = FALSE)
ss <- read_tsv(here("tgct", "ss", "pennsamples_masterslist_062226.txt"))
dim(ss)
# 1350

cov <- fread(here("PMBB", "4.0", "PMBB-Release-2026-4.0_phenotype_covariates.txt"), header = TRUE)
flags <- fread(here("PMBB", "4.0", "rgcname_pmbbid_metadata_flags_freeze4.csv"), header = TRUE)

# ids <- fam$V2
# length(ids)
# # 915

# ### GET CONSENTED
# table(flags$PMBB_consent)
# # current    ended     null withdraw
# # 63411     3479     4030        5
#
# consent <- flags %>% filter(PMBB_consent == "current")
# consent_ids <- consent$PMBB_ID
# length(consent_ids)
# # 63411

# ========================
# SELECT IDs
# ========================
ss2 <- ss %>%
    left_join(flags %>% select(PMBB_ID, Freeze, PMBB_consent), by = c("PMBBID2" = "PMBB_ID")) %>%
    # filter(Geno.Avail == TRUE | is.na(Geno.Avail)) %>%
    mutate(Final_ID = coalesce(PMBBID2, GencoveID, ReplicationID, DiscoveryID)) %>%
    mutate(GenoSource2 = case_when(
        !is.na(PMBBID2) ~ "PMBB",
        !is.na(GencoveID) ~ "Gencove",
        !is.na(ReplicationID) ~ "Replication",
        !is.na(DiscoveryID) ~ "Discovery"
    )) %>%
    filter(!is.na(GenoSource2))
dim(ss2)
# 1133

ids <- ss2$Final_ID

ids_out <- ss2 %>%
    mutate(write_id = if_else(GenoSource2 == "Replication", paste0(Final_ID, "_", Final_ID), Final_ID)) %>%
    transmute(FID = write_id, IID = write_id)

write_tsv(ids_out, here("tgct", "data", "tgct_prs_ids.txt"), col_names = FALSE)
table(ss2$GenoSource2)
# Discovery     Gencove        PMBB Replication
# 67          81         922          66

# Freeze 1 Freeze 2 Freeze 3 Freeze 4
# 25      126      768        3

write_xlsx(ss2, here("tgct", "ss", "20260710_tgct_master.xlsx"))

### additional...
ss_pmbb <- ss2 %>% filter(GenoSource2 == "PMBB")
dim(ss_pmbb)
# 921

ss_pmbb2 <- ss2 %>% filter(PMBBID2 %in% ids)
dim(ss_pmbb2)
# 914

table(ss2$Freeze)
# Freeze 1 Freeze 2 Freeze 3 Freeze 4
# 25      126      767        3

# ========================
# EXPORT KEEP FILES FOR PLINK2
# ========================
table(ss2$GenoSource2)
# DiscoveryID     GencoveID PMBBID2     ReplicationID
# 65            81           921            66

for (source in unique(ss2$GenoSource2)) {
    df <- ss2 %>% filter(GenoSource2 == source)
    id_col <- df[[paste0(source, "ID")]]
    if (source == "Replication") {
        id_col <- paste0(id_col, "_", id_col)
    }
    out <- tibble(FID = id_col, IID = id_col)
    write_tsv(out, here("tgct", "data", paste0("keep_", source, ".txt")), col_names = FALSE)
}










