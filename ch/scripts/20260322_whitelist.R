# ========================
# Whitelist filtering
# ========================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
library(data.table, quietly = T)

# ========================
# READ IN DATA
# ========================
all_ch <- read.csv(file.path("ch", "data", "ch_seq_vars.csv"), row.names = 1)
dim(all_ch)

gList<-fread(file.path("ch", "data", "whitelist_filter_20230531", "Full_CHIP_gene_list_20260324.txt"))
whitelist.mis<-fread(file.path("ch", "data", "whitelist_filter_20230531", "CHIP_missense_vars_cv_20260324.txt"))
whitelist.splice<-fread(file.path("ch", "data", "whitelist_filter_20230531", "CHIP_splice_vars_cv_20260324.txt"))
whitelist.LoF<-fread(file.path("ch", "data", "whitelist_filter_20230531", "CHIP_nonsense_FS_vars_cv_20260324.txt"))

# ========================
# MATCH TRANSCRIPTS
# ========================
all_ch <- all_ch %>% mutate(across(where(is.character), ~ URLdecode(.x)))

all_ch <- all_ch %>% left_join(gList[, c("Gene", "Accession")], by = "Gene")
all_ch <- all_ch %>%
    dplyr::mutate(
        MANE_Select_Str   = gsub("\\.\\d+$", "", MANE_SELECT),
        Accession_Str     = gsub("\\.\\d+$", "", Accession),
        correct_transcript     = (MANE_Select_Str == Accession_Str)
    )

table(all_ch$correct_transcript)
# FALSE  TRUE
# 5431 24173

# examine transcript mismatch
wrong_transcript <- all_ch %>% filter(!correct_transcript)
table(wrong_transcript$Gene)
# BCORL1  BRCC3  CSF1R   EZH2  GATA2  IKZF2  KDM6A  KMT2A   KRAS    NF1  RUNX1 SETDB1  SRSF2  STAG2   TP53
# 287     33    346    175    195    110    179    770     71    594    217    590    139    862    642

# it seems their RefSeq isn't exactly our MANE...
examine <- all_ch %>%
    filter(!correct_transcript) %>%
    distinct(Gene, MANE_SELECT, Accession) %>%
    print()

write.csv(examine,  file.path("ch", "data", "ch_gene_transcript_mismatches.csv"),   row.names = FALSE)

# ========================
# REMOVE LAST EXON
# ========================
extract_exon_info <- function(exon_field) {
    if (is.na(exon_field) || exon_field == "" || exon_field == ".") {
        return(list(exon_num = NA_real_, exon_total = NA_real_))
    }
    parts <- strsplit(as.character(exon_field), "\\|")[[1]]
    list(exon_num = as.numeric(parts[1]), exon_total = as.numeric(parts[2]))
}

all_ch <- all_ch %>%
    mutate(
        ExonNumber   = sapply(EXON, function(x) extract_exon_info(x)$exon_num),
        ExonTotal    = sapply(EXON, function(x) extract_exon_info(x)$exon_total),
        is_last_exon = !is.na(ExonNumber) & !is.na(ExonTotal) & ExonNumber == ExonTotal
    )

cat("Variants in last exon:", sum(all_ch$is_last_exon, na.rm = TRUE), "\n")
# 6714
all_ch <- all_ch %>% filter(!is_last_exon)

cat("Variants after last exon filter:", nrow(all_ch), "\n")
# 22890

# ========================
# EXTRACT SPECIFIC FIELDS
# ========================
aa_3to1 <- c(
    Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C",
    Gln = "Q", Glu = "E", Gly = "G", His = "H", Ile = "I",
    Leu = "L", Lys = "K", Met = "M", Phe = "F", Pro = "P",
    Ser = "S", Thr = "T", Trp = "W", Tyr = "Y", Val = "V",
    Ter = "*"
)

# Extract protein change from HGVSp (everything after "p.")
extractProteinChange <- function(hgvsp) {
    if (is.na(hgvsp) || hgvsp == "" || hgvsp == ".") return(NA_character_)
    prot_change <- gsub("^p\\.", "", hgvsp)
    prot_change <- gsub("%3D", "=", prot_change)
    return(prot_change)
}

# Replace each 3-letter code with its 1-letter equivalent
convert_3to1 <- function(prot_change) {
    if (is.na(prot_change) || prot_change == "" || prot_change == ".") return(NA_character_)

    result <- prot_change
    for (aa3 in names(aa_3to1)) {
        result <- gsub(aa3, aa_3to1[[aa3]], result, fixed = TRUE)
    }
    return(result)
}

# Extract numeric AA position from 1-letter change for range-based filters
extract_aa_position <- function(prot_change) {
    if (is.na(prot_change) || prot_change == "" || prot_change == ".") return(NA_real_)
    as.numeric(regmatches(prot_change, regexpr("\\d+", prot_change)))
}

all_ch$ProteinChange <- sapply(all_ch$HGVSp, extractProteinChange)
all_ch$ProteinChange_1L <- sapply(all_ch$ProteinChange, convert_3to1)
all_ch$AAPosition <- sapply(all_ch$ProteinChange_1L, extract_aa_position)

table(all_ch$Variant.LoF_level)
# 1     2     3     4
# 980  1917 10397  9596

write.csv(all_ch, file.path("ch", "data", "ch_seq_vars_proc.csv"),   row.names = FALSE)
all_ch2 <- all_ch
all_ch <- all_ch2

# ========================
# INITIALIZE WHITELIST FLAGS
# ========================
all_ch <- all_ch %>%
    dplyr::mutate(
        whitelist    = FALSE,
        wl.mis       = FALSE,
        wl.lof       = FALSE,
        wl.splice    = FALSE,
        wl.exception = FALSE,
        manualreview = FALSE
    )

# ========================
# 1) MISSENSE
# Match on 1-letter ProteinChange_1L directly against whitelist AAChange
# ========================
vmis <- grepl("missense", all_ch$Variant.Consequence)
vCorrectTranscript <- all_ch$correct_transcript

vmis_wl <- paste(all_ch$Gene, all_ch$ProteinChange_1L, sep = "_") %in%
    paste(whitelist.mis$Gene, whitelist.mis$AAChange, sep = "_")

# all_ch$whitelist[vmis & vmis_wl & vCorrectTranscript] <- TRUE
# all_ch$wl.mis[vmis & vmis_wl & vCorrectTranscript] <- TRUE
# cat(sprintf("Missense variants in whitelist: %d\n", sum(vmis & vmis_wl & vCorrectTranscript)))
# # 57 are fine on correct transcript...
#
# mis_manual <- all_ch[vmis & vmis_wl & !vCorrectTranscript, ]
# dim(mis_manual)
# # 15
# write_xlsx(mis_manual, file.path("ch", "data", "ch_mis_check.xlsx"))
# # these are all fine -- different transcript doesn't affect the varinat

all_ch$whitelist[vmis & vmis_wl] <- TRUE
all_ch$wl.mis[vmis & vmis_wl] <- TRUE
cat(sprintf("Missense variants in whitelist: %d\n", sum(vmis & vmis_wl)))
# 72

write_xlsx(all_ch[vmis & vmis_wl, ], file.path("ch", "data", "ch_wl_mis_checked.xlsx"))

# ========================
# 2) LOF / FRAMESHIFT
# ========================
vlof <- grepl("Ter|\\*|X", all_ch$ProteinChange_1L) |
    grepl("stop_gain|stop_lost|start_lost", all_ch$Variant.Consequence) |
    grepl("fs", all_ch$ProteinChange_1L, fixed = TRUE) |
    grepl("frameshift", all_ch$Variant.Consequence)

vLOFgene <- all_ch$Gene %in% whitelist.LoF$Gene

all_ch$whitelist[vlof & vLOFgene] <- TRUE
all_ch$wl.lof[vlof & vLOFgene]    <- TRUE
cat(sprintf("LoF/frameshift variants in whitelist: %d\n", sum(vlof & vLOFgene, na.rm = TRUE)))
# 420

### MANUAL CHECK - INDIVIDUALS W >1 VARIANT
lof <- all_ch[all_ch$wl.lof, ]

# find samples with >1 row within 50 bp of each other
# these could be complex indel events split up
lof_check <- lof %>%
    group_by(Sample.ID, Gene) %>%
    filter(n() > 1) %>%
    arrange(Sample.ID, Gene, Start) %>%
    mutate(
        near_prev = c(FALSE, diff(Start) <= 50)
    ) %>%
    filter(near_prev | lead(near_prev, default = FALSE)) %>%
    ungroup() %>%
    arrange(Sample.ID, Gene, Start)
dim(lof_check)
# 195

### MAKE CSV FOR BRAD IGV ###
lof_check_out <- lof_check %>% group_by(Sample.ID, Gene) %>% summarise(Chr = unique(Chr), min_start = min(Start), max_start = max(Start)) %>% filter(min_start != max_start)

write.csv(lof_check_out, file.path("ch", "data", "ch_lof_examine.csv"), row.names = FALSE)

### FOR NOW TAKE THE LEFTMOST + MANUAL CHECK ###
lof_check_man <- lof_check %>%
    group_by(Sample.ID, Gene) %>%
    mutate(Keep = Start == min(Start)) %>%
    ungroup()

write_xlsx(lof_check_man, file.path("ch", "data", "ch_wl_lof_check.xlsx"))
# manually sorted these into ch_wl_lof_checked.xlsx

# ### ONE VARIANT / INDIVIDUAL
# lof_one <- lof %>% filter(!(Sample.ID %in% lof_check$Sample.ID))
# cat(sprintf("Individuals with potential duplicate missense: %d\n", n_distinct(paste(lof_one$Sample.ID, lof_one$Gene))))
# # 213
#
# # not overlap check - good
# sum(unique(lof_one$Sample.ID) %in% unique(lof_check$Sample.ID))
#
# table(lof_one$Variant.LoF_level)
# # 1   2   3
# # 189  23   2
# # these are prob fine..so keep them as is

# ========================
# 3) SPLICE
# ========================
vSplice <- grepl("splice", all_ch$Variant.Consequence) & (all_ch$Variant.LoF_level %in% c(1, 2))
vSplicegene <- all_ch$Gene %in% whitelist.splice$Gene

vSpliceCorrectTranscript <- all_ch$correct_transcript
cat(sprintf("Splice variants in whitelist: %d\n", sum(vSplice & vSplicegene & vSpliceCorrectTranscript)))
# 99
cat(sprintf("Splice variants in whitelist: %d\n", sum(vSplice & vSplicegene & !vSpliceCorrectTranscript)))
# 138

all_ch$whitelist[vSplice & vSplicegene & vSpliceCorrectTranscript] <- TRUE
all_ch$wl.splice[vSplice & vSplicegene & vSpliceCorrectTranscript] <- TRUE
all_ch$manualreview[(vSplice & vSplicegene) & (!vSpliceCorrectTranscript)] <- TRUE

write_xlsx(all_ch[all_ch$wl.splice, ], file.path("ch", "data", "ch_wl_splice_check.xlsx"))
# manually checked these into ch_wl_splice_checked.xlsx

# ========================
# GENE-SPECIFIC RULES
# ========================
vFS <- grepl("fs", all_ch$ProteinChange_1L, fixed = TRUE) | grepl("frameshift_variant", all_ch$Variant.Consequence)
vexon11 <- all_ch$ExonNumber == 11
vexon12 <- all_ch$ExonNumber == 12
vexon5  <- all_ch$ExonNumber == 5
vexon6  <- all_ch$ExonNumber == 6

# ASXL1: frameshift/nonsense/splice in exon 11-12
asxl1Exception <- (all_ch$Gene == "ASXL1") & (vlof | vFS) & (vexon11 | vexon12)
all_ch$whitelist[asxl1Exception]    <- TRUE
all_ch$wl.lof[asxl1Exception]       <- TRUE
all_ch$wl.exception[asxl1Exception] <- TRUE

asxl1ExceptionSplice <- (all_ch$Gene == "ASXL1") & vSplice & vSpliceCorrectTranscript & (vexon11 | vexon12)
all_ch$whitelist[asxl1ExceptionSplice]    <- TRUE
all_ch$wl.splice[asxl1ExceptionSplice]    <- TRUE
all_ch$wl.exception[asxl1ExceptionSplice] <- TRUE
cat(sprintf("ASXL1 exceptions: %d\n", sum(asxl1Exception | asxl1ExceptionSplice, na.rm = TRUE)))
# 4

# ASXL2: frameshift/nonsense/splice in exon 11-12
asxl2Exception <- (all_ch$Gene == "ASXL2") & (vlof | vFS) & (vexon11 | vexon12)
all_ch$whitelist[asxl2Exception]    <- TRUE
all_ch$wl.lof[asxl2Exception]       <- TRUE
all_ch$wl.exception[asxl2Exception] <- TRUE

asxl2ExceptionSplice <- (all_ch$Gene == "ASXL2") & vSplice & vSpliceCorrectTranscript & (vexon11 | vexon12)
all_ch$whitelist[asxl2ExceptionSplice]    <- TRUE
all_ch$wl.splice[asxl2ExceptionSplice]    <- TRUE
all_ch$wl.exception[asxl2ExceptionSplice] <- TRUE
cat(sprintf("ASXL2 exceptions: %d\n", sum(asxl2Exception | asxl2ExceptionSplice, na.rm = TRUE)))
# 8

# PPM1D: frameshift/nonsense in exon 5-6
ppm1dException <- (all_ch$Gene == "PPM1D") & (vlof | vFS) & (vexon5 | vexon6)
all_ch$whitelist[ppm1dException]    <- TRUE
all_ch$wl.lof[ppm1dException]       <- TRUE
all_ch$wl.exception[ppm1dException] <- TRUE
cat(sprintf("PPM1D exceptions: %d\n", sum(ppm1dException, na.rm = TRUE)))
# 0

# TET2: missense in catalytic domains (p.1104-1481 and p.1843-2002)
vmis <- grepl("missense", all_ch$Variant.Consequence)
TETidx <- which(all_ch$Gene == "TET2" & vmis & !is.na(all_ch$AAPosition))
for (i in TETidx) {
    AApos <- all_ch$AAPosition[i]
    if (!is.na(AApos) && ((AApos >= 1104 & AApos <= 1481) | (AApos >= 1843 & AApos <= 2002))) {
        all_ch$whitelist[i]    <- TRUE
        all_ch$wl.mis[i]       <- TRUE
        all_ch$wl.exception[i] <- TRUE
    }
}
cat(sprintf("TET2 catalytic domain exceptions: %d\n", sum(all_ch$wl.exception & all_ch$Gene == "TET2")))
# 25

# CBL: RING finger missense p.381-421
CBLidx <- which(all_ch$Gene == "CBL" & vmis & !is.na(all_ch$AAPosition))
for (i in CBLidx) {
    AApos <- all_ch$AAPosition[i]
    if (!is.na(AApos) && AApos >= 381 & AApos <= 421) {
        all_ch$whitelist[i]    <- TRUE
        all_ch$wl.mis[i]       <- TRUE
        all_ch$wl.exception[i] <- TRUE
    }
}
cat(sprintf("CBL RING finger exceptions: %d\n", sum(all_ch$wl.exception & all_ch$Gene == "CBL")))
# 1

# CBLB: RING finger missense p.372-412
CBLBidx <- which(all_ch$Gene == "CBLB" & vmis & !is.na(all_ch$AAPosition))
for (i in CBLBidx) {
    AApos <- all_ch$AAPosition[i]
    if (!is.na(AApos) && AApos >= 372 & AApos <= 412) {
        all_ch$whitelist[i]    <- TRUE
        all_ch$wl.mis[i]       <- TRUE
        all_ch$wl.exception[i] <- TRUE
    }
}
cat(sprintf("CBLB RING finger exceptions: %d\n", sum(all_ch$wl.exception & all_ch$Gene == "CBLB")))
# 2

gene_exception <- all_ch %>% filter(wl.exception == TRUE)
table(gene_exception$correct_transcript)
# 40 these are all correct

write_xlsx(all_ch[all_ch$wl.exception, ], file.path("ch", "data", "ch_wl_gene_checked.xlsx"))

# ========================
# MANUAL REVIEW
# ========================
vNonFrameshiftDel <- grepl("inframe_deletion",  all_ch$Variant.Consequence)
vNonFrameshiftIns <- grepl("inframe_insertion", all_ch$Variant.Consequence)
vFrameshiftIndel  <- grepl("frameshift_variant", all_ch$Variant.Consequence)

all_ch$manualreview[(vlof | vFS | vSplice) & (all_ch$Gene == "GATA3")]                    <- TRUE
all_ch$manualreview[vNonFrameshiftDel       & (all_ch$Gene == "CREBBP")]                   <- TRUE
all_ch$manualreview[(vlof | vFS | vSplice) & (all_ch$Gene == "CSF3R")]                    <- TRUE
all_ch$manualreview[vNonFrameshiftDel       & (all_ch$Gene == "DNMT3A")]                   <- TRUE
all_ch$manualreview[vNonFrameshiftDel       & (all_ch$Gene == "EP300")]                    <- TRUE
all_ch$manualreview[(vNonFrameshiftDel | vNonFrameshiftIns) & (all_ch$Gene == "FLT3")]    <- TRUE
all_ch$manualreview[(vNonFrameshiftDel | vNonFrameshiftIns) & (all_ch$Gene == "JAK2")]    <- TRUE
all_ch$manualreview[vNonFrameshiftDel       & (all_ch$Gene == "KDM6A")]                   <- TRUE
all_ch$manualreview[(vNonFrameshiftDel | vNonFrameshiftIns) & (all_ch$Gene == "KIT")]     <- TRUE
all_ch$manualreview[(vNonFrameshiftDel | vNonFrameshiftIns) & (all_ch$Gene == "MPL")]     <- TRUE
all_ch$manualreview[vFrameshiftIndel& (all_ch$Gene == "NPM1")]                    <- TRUE

all_ch$manualreview[vmis & vmis_wl & !vCorrectTranscript] <- TRUE

cat(sprintf("Total flagged for manual review: %d\n", sum(all_ch$manualreview, na.rm = TRUE)))
# 175
write_xlsx(all_ch[all_ch$manualreview, ], file.path("ch", "data", "ch_wl_manual.xlsx"))
# sorted these into ch_wl_manual_checked.xlsx

# ========================
# COMBINE MANUAL REVIEWS
# ========================
mis_checked <- read_excel(file.path("ch", "data", "ch_wl_mis_checked.xlsx"))
lof_checked <- read_excel(file.path("ch", "data", "ch_wl_lof_checked.xlsx"))
spl_checked <- read_excel(file.path("ch", "data", "ch_wl_splice_checked.xlsx"), sheet = "keep")
gene_checked <- read_excel(file.path("ch", "data", "ch_wl_gene_checked.xlsx"))
manual_checked <- read_excel(file.path("ch", "data", "ch_wl_manual_checked.xlsx"))

ch_final <- rbind(lof_checked %>% filter(Keep == 1) %>% select(-Keep, -near_prev),
                   gene_checked,
                   mis_checked,
                   spl_checked,
                   manual_checked %>% filter(Keep == 1) %>% select(-Keep, -Artifact)
                   )

write_xlsx(ch_final, file.path("ch", "data", "ch_wl_final.xlsx"))
dim(ch_final)
# 315

length(unique(ch_final$Sample.ID))
# 273

