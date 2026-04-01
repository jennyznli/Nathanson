# ========================
# Whitelist filtering
# ========================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))

# ========================
# READ IN DATA
# ========================
all_ch <- read.csv(file.path("ch", "data", "ch_seq_vars.csv"), row.names = 1)
# all_ch <- read.csv(file.path("ch", "data", "ch_all_vars.csv"), row.names = 1)

dim(all_ch)
# 29656

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
# 5436 24224

# examine transcript mismatch
wrong_transcript <- all_ch %>% filter(!correct_transcript)
table(wrong_transcript$Gene)
# BCORL1  BRCC3  CSF1R   EZH2  GATA2  IKZF2  KDM6A  KMT2A   KRAS    NF1  RUNX1 SETDB1  SRSF2  STAG2   TP53
# 307     33    387    178    211    152    179    774     93    604    277    593    139    863    646

examine <- all_ch %>%
    filter(!correct_transcript) %>%
    distinct(Gene, MANE_SELECT, Accession) %>%
    print()
examine
# Gene    MANE_SELECT    Accession
# 1    TP53    NM_000546.6 NM_001126112
# 2   GATA2    NM_032638.5 NM_001145661
# 3    EZH2    NM_004456.5 NM_001203247
# 4     NF1 NM_001042492.3    NM_000267
# 5  SETDB1 NM_001366418.1 NM_001145415
# 6   RUNX1    NM_001754.5 NM_001001890
# 7   KMT2A NM_001197104.2    NM_005933
# 8   IKZF2 NM_001387220.1    NM_016260
# 9   CSF1R NM_001288705.3    NM_005211
# 10  SRSF2 NM_001195427.2    NM_003016
# 11  STAG2 NM_001042750.2    NM_006603
# 12   KRAS    NM_004985.5    NM_033360
# 13  KDM6A NM_001291415.2    NM_021140
# 14 BCORL1 NM_001379451.1    NM_021946
# 15  BRCC3 NM_001018055.3    NM_024332

write_xlsx(examine, file.path("ch", "data", "ch_gene_transcript_mismatches.xlsx"))

# it seems their RefSeq isn't exactly our MANE...
# we'll keep using MANE and manually examine mismatching transcripts

# ========================
# MARK LAST EXON
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
# 6726

# ========================
# EXTRACT SPECIFIC FIELDS
# ========================
aa_3to1 <- c(
    Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C",
    Gln = "Q", Glu = "E", Gly = "G", His = "H", Ile = "I",
    Leu = "L", Lys = "K", Met = "M", Phe = "F", Pro = "P",
    Ser = "S", Thr = "T", Trp = "W", Tyr = "Y", Val = "V",
    fsTer = "fs*", Ter = "X"
)

# Extract protein change from HGVSp
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

all_ch$ProteinChange <- sapply(all_ch$HGVSp, extractProteinChange)
all_ch$ProteinChange_1L <- sapply(all_ch$ProteinChange, convert_3to1)
all_ch$Protein.position <- as.numeric(all_ch$Protein.position)
all_ch$cDNA.position <- as.numeric(all_ch$cDNA.position)
all_ch$CDS.position <- as.numeric(all_ch$CDS.position)

write.csv(all_ch, file.path("ch", "data", "ch_seq_vars_proc.csv"),   row.names = FALSE)
all_ch2 <- all_ch
# all_ch <- all_ch2

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

# x <- as.data.frame(paste(whitelist.mis$Gene, whitelist.mis$AAChange, sep = "_"))

vmis_wl <- paste(all_ch$Gene, all_ch$ProteinChange_1L, sep = "_") %in%
    paste(whitelist.mis$Gene, whitelist.mis$AAChange, sep = "_")
# correct_transcript <- all_ch$correct_transcript

all_ch$whitelist[vmis & vmis_wl & correct_transcript] <- TRUE
all_ch$wl.mis[vmis & vmis_wl & correct_transcript] <- TRUE

cat(sprintf("Missense variants in whitelist: %d\n", sum(vmis & vmis_wl)))
# 87

# could consider just getting variants in same position not necessarily same AA changes..

# ========================
# 2) LOF / FRAMESHIFT
# ========================
vlof <- grepl("Ter|\\*|X|fs*", all_ch$ProteinChange_1L) |
    grepl("stop_gain|stop_lost|start_lost", all_ch$Variant.Consequence) |
    grepl("fs", all_ch$ProteinChange_1L, fixed = TRUE) |
    grepl("frameshift", all_ch$Variant.Consequence)

vLOFgene <- all_ch$Gene %in% whitelist.LoF$Gene

all_ch$whitelist[vlof & vLOFgene] <- TRUE
all_ch$wl.lof[vlof & vLOFgene] <- TRUE
cat(sprintf("LoF/frameshift variants: %d\n", sum(vlof & vLOFgene, na.rm = TRUE)))
# 517

# ========================
# 3) SPLICE
# ========================
#  & (all_ch$Variant.LoF_level %in% c(1, 2))
vSplice <- grepl("splice", all_ch$Variant.Consequence)
vSpliceStringent <- grepl("splice", all_ch$Variant.Consequence) & (all_ch$Variant.LoF_level %in% c(1, 2))
vSplicegene <- all_ch$Gene %in% whitelist.splice$Gene

all_ch$whitelist[vSpliceStringent & vSplicegene] <- TRUE
all_ch$wl.splice[vSpliceStringent & vSplicegene] <- TRUE

cat(sprintf("Splice variants in whitelist: %d\n", sum(vSpliceStringent & vSplicegene)))
# 237

# ========================
# GENE-SPECIFIC RULES
# ========================
vFS <- grepl("fs", all_ch$ProteinChange_1L, fixed = TRUE) | grepl("frameshift_variant", all_ch$Variant.Consequence)
# vexon11 <- all_ch$ExonNumber == 11
vexon12 <- all_ch$ExonNumber == 12
vexon13 <- all_ch$ExonNumber == 13
vexon5  <- all_ch$ExonNumber == 5
vexon6  <- all_ch$ExonNumber == 6

# ASXL1: frameshift/nonsense/splice in exon 12-13 NOT 11-12!!!
asxl1Exception <- (all_ch$Gene == "ASXL1") & (vlof | vFS) & (vexon12 | vexon13)
all_ch$whitelist[asxl1Exception]    <- TRUE
all_ch$wl.lof[asxl1Exception]       <- TRUE
all_ch$wl.exception[asxl1Exception] <- TRUE
cat(sprintf("ASXL1 frameshift exceptions: %d\n", sum(asxl1Exception, na.rm = TRUE)))
# 26

asxl1ExceptionSplice <- (all_ch$Gene == "ASXL1") & vSplice & (vexon12 | vexon13)
all_ch$whitelist[asxl1ExceptionSplice]    <- TRUE
all_ch$wl.splice[asxl1ExceptionSplice]    <- TRUE
all_ch$wl.exception[asxl1ExceptionSplice] <- TRUE
cat(sprintf("ASXL1 splice exceptions: %d\n", sum(asxl1ExceptionSplice, na.rm = TRUE)))
# 0

# ASXL2: frameshift/nonsense/splice in in exon 12-13 NOT 11-12!!!
asxl2Exception <- (all_ch$Gene == "ASXL2") & (vlof | vFS) & (vexon12 | vexon13)
all_ch$whitelist[asxl2Exception]    <- TRUE
all_ch$wl.lof[asxl2Exception]       <- TRUE
all_ch$wl.exception[asxl2Exception] <- TRUE
cat(sprintf("ASXL2 frameshift exceptions: %d\n", sum(asxl2Exception, na.rm = TRUE)))
# 11

asxl2ExceptionSplice <- (all_ch$Gene == "ASXL2") & vSplice & (vexon12 | vexon13)
all_ch$whitelist[asxl2ExceptionSplice]    <- TRUE
all_ch$wl.splice[asxl2ExceptionSplice]    <- TRUE
all_ch$wl.exception[asxl2ExceptionSplice] <- TRUE
cat(sprintf("ASXL2 exceptions: %d\n", sum(asxl2ExceptionSplice, na.rm = TRUE)))
# 1

# PPM1D: frameshift/nonsense in exon 5-6
ppm1dException <- (all_ch$Gene == "PPM1D") & (vlof | vFS) & (vexon5 | vexon6)
all_ch$whitelist[ppm1dException]    <- TRUE
all_ch$wl.lof[ppm1dException]       <- TRUE
all_ch$wl.exception[ppm1dException] <- TRUE
cat(sprintf("PPM1D exceptions: %d\n", sum(ppm1dException, na.rm = TRUE)))
# 3

# View(all_ch[ppm1dException, ])
# View(all_ch %>% filter(Gene == "TP53", ExonNumber %in% c(5, 6)))
# View(all_ch %>% filter(Gene == "TP53", Protein.position < 400, Protein.position > 100,
#                        Variant.Consequence == "missense_variant") %>% arrange(Protein.position))
#
# View(all_ch %>% filter(Gene == "SRSF2") %>% arrange(Protein.position))
# View(all_ch %>% filter(Gene == "JAK2", wl.exception) %>% arrange(Protein.position))


# TET2: missense in catalytic domains (p.1104-1481 and p.1843-2002)
vmis <- grepl("missense", all_ch$Variant.Consequence)
TETidx <- which(all_ch$Gene == "TET2" & vmis & !is.na(all_ch$Protein.position))
for (i in TETidx) {
    AApos <- all_ch$Protein.position[i]
    if (!is.na(AApos) && ((AApos >= 1104 & AApos <= 1481) | (AApos >= 1843 & AApos <= 2002))) {
        all_ch$whitelist[i]    <- TRUE
        all_ch$wl.mis[i]       <- TRUE
        all_ch$wl.exception[i] <- TRUE
    }
}
cat(sprintf("TET2 catalytic domain exceptions: %d\n", sum(all_ch$wl.exception & all_ch$Gene == "TET2")))
# 25

# CBL: RING finger missense p.381-421 - correct transcript
# View(all_ch %>% filter(Gene == "CBL") %>% select(correct_transcript))
CBLidx <- which(all_ch$Gene == "CBL" & vmis & !is.na(all_ch$Protein.position))
for (i in CBLidx) {
    AApos <- all_ch$Protein.position[i]
    if (!is.na(AApos) && AApos >= 381 & AApos <= 421) {
        all_ch$whitelist[i]    <- TRUE
        all_ch$wl.mis[i]       <- TRUE
        all_ch$wl.exception[i] <- TRUE
    }
}
cat(sprintf("CBL RING finger exceptions: %d\n", sum(all_ch$wl.exception & all_ch$Gene == "CBL")))
# 1

# CBLB: RING finger missense p.372-412 - correct transcript
# View(all_ch %>% filter(Gene == "CBLB") %>% select(correct_transcript))
CBLBidx <- which(all_ch$Gene == "CBLB" & vmis & !is.na(all_ch$Protein.position))
for (i in CBLBidx) {
    AApos <- all_ch$Protein.position[i]
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
# 88 - these are all correct transcript

# ========================
# SUMMARIZE GENE-SPECIFIC EXCEPTIONS
# ========================
exception_summary <- all_ch %>%
    filter(wl.exception) %>%
    group_by(Gene) %>%
    summarise(
        n_lof    = sum(wl.lof,    na.rm = TRUE),
        n_mis    = sum(wl.mis,    na.rm = TRUE),
        n_splice = sum(wl.splice, na.rm = TRUE),
        n_total  = n(),
        .groups  = "drop"
    ) %>%
    arrange(desc(n_total))

print(exception_summary)
# Gene  n_lof n_mis n_splice n_total
# <chr> <int> <int>    <int>   <int>
#     1 TET2      0    44        0      44
# 2 ASXL1    26     0        0      26
# 3 ASXL2    11     0        1      12
# 4 PPM1D     3     0        0       3
# 5 CBLB      0     2        0       2
# 6 CBL       0     1        0       1

write.csv(exception_summary,
          file.path("ch", "data", "ch_gene_exception_summary.csv"),
          row.names = FALSE)

# Plot
exception_long <- exception_summary %>%
    tidyr::pivot_longer(cols = c(n_lof, n_mis, n_splice),
                        names_to  = "type",
                        values_to = "count") %>%
    mutate(
        type = factor(type,
                      levels = c("n_lof", "n_mis", "n_splice"),
                      labels = c("LoF / frameshift", "Missense", "Splice")),
        Gene = factor(Gene, levels = exception_summary$Gene)
    )

p_exceptions <- ggplot(exception_long, aes(x = Gene, y = count, fill = type)) +
    geom_col(position = "stack", width = 0.65) +
    scale_fill_manual(
        values = c("LoF / frameshift" = "#1D9E75",
                   "Missense"         = "#378ADD",
                   "Splice"           = "#D85A30"),
        name = NULL
    ) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 6),
                       expand = expansion(mult = c(0, 0.1))) +
    labs(title    = "Gene-specific whitelist exceptions",
         subtitle = sprintf("Total: %d exceptions across %d genes",
                            sum(exception_summary$n_total),
                            nrow(exception_summary)),
         x = NULL, y = "Exception count") +
    theme_minimal(base_size = 11) +
    theme(
        plot.title       = element_text(face = "bold", hjust = 0.5),
        plot.subtitle    = element_text(hjust = 0.5, size = 9, color = "grey40"),
        legend.position  = "bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank()
    )

ggsave(file.path("ch", "figures", "qc_gene_exception_counts.pdf"),
       p_exceptions, width = 7, height = 5)

# ========================
# MANUAL REVIEW
# ========================
vNonFrameshiftDel <- grepl("inframe_deletion",  all_ch$Variant.Consequence)
vNonFrameshiftIns <- grepl("inframe_insertion", all_ch$Variant.Consequence)
vFrameshiftIndel  <- grepl("frameshift_variant", all_ch$Variant.Consequence)

# loss of function / splice
all_ch$manualreview[(vlof | vFS | vSplice) & (all_ch$Gene == "GATA3")]                    <- TRUE
all_ch$manualreview[(vlof | vFS | vSplice) & (all_ch$Gene == "CSF3R")]                    <- TRUE

# non frameshift
all_ch$manualreview[vNonFrameshiftDel       & (all_ch$Gene == "CREBBP")]                   <- TRUE
all_ch$manualreview[vNonFrameshiftDel       & (all_ch$Gene == "DNMT3A")]                   <- TRUE
all_ch$manualreview[vNonFrameshiftDel       & (all_ch$Gene == "EP300")]                    <- TRUE
all_ch$manualreview[(vNonFrameshiftDel | vNonFrameshiftIns) & (all_ch$Gene == "FLT3")]    <- TRUE
all_ch$manualreview[(vNonFrameshiftDel | vNonFrameshiftIns) & (all_ch$Gene == "JAK2")]    <- TRUE
all_ch$manualreview[vNonFrameshiftDel       & (all_ch$Gene == "KDM6A")]                   <- TRUE
all_ch$manualreview[(vNonFrameshiftDel | vNonFrameshiftIns) & (all_ch$Gene == "KIT")]     <- TRUE
all_ch$manualreview[(vNonFrameshiftDel | vNonFrameshiftIns) & (all_ch$Gene == "MPL")]     <- TRUE

# frameshifts
all_ch$manualreview[vFrameshiftIndel& (all_ch$Gene == "NPM1")]                    <- TRUE

# missense without correct transcript
all_ch$manualreview[vmis & vmis_wl & !correct_transcript] <- TRUE

cat(sprintf("Gene specific flagged for manual review: %d\n", sum(all_ch$manualreview, na.rm = TRUE)))
# 127

# filter down to any variant of interest
all_ch2 <- all_ch %>% filter(whitelist | manualreview | wl.exception)

# View(all_ch_raw %>% filter(Gene == "JAK2"))
# View(all_ch %>% filter(Gene == "JAK2"))
#
View(all_ch %>% filter(Gene == "SF3B1"))
# View(all_ch2 %>% filter(Gene == "JAK2"))

# ========================
# ASSIGN VARIANT CATEGORY
# ========================
vlof <- grepl("Ter|\\*|X|fs*", all_ch2$ProteinChange_1L) |
    grepl("stop_gain|stop_lost|start_lost", all_ch2$Variant.Consequence) |
    grepl("fs", all_ch2$ProteinChange_1L, fixed = TRUE) |
    grepl("frameshift", all_ch2$Variant.Consequence)
vmis <- grepl("missense", all_ch2$Variant.Consequence)
vSplice <- grepl("splice", all_ch2$Variant.Consequence) & (all_ch2$Variant.LoF_level %in% c(1, 2))
inframe <- grepl("inframe", all_ch2$Variant.Consequence)

all_ch2 <- all_ch2 %>%
    dplyr::mutate(
        variant_category = dplyr::case_when(
            vlof   ~ "lof",
            vmis   ~ "missense",
            vSplice ~ "splice",
            inframe ~ "inframe",
            TRUE    ~ NA_character_
        )
    )

print(table(all_ch2$variant_category, useNA = "always"))
# inframe      lof missense   splice     <NA>
#     95      529      137      232        6

# ========================
# CREATE VARIANT ID
# ========================
trim_splice_hgvsc <- function(hgvsc) {
    if (is.na(hgvsc) || hgvsc == "" || hgvsc == ".") return(NA_character_)
    hgvsc_clean <- sub("^[A-Z]{2}_[0-9]+(?:\\.[0-9]+)?:", "", hgvsc)
    sub("(c\\.\\d+)[+\\-].*", "\\1", hgvsc_clean)
}

clean_hgvsc <- function(hgvsc) {
    if (is.na(hgvsc) || hgvsc == "" || hgvsc == ".") return(NA_character_)
    sub("^[A-Z]{2}_[0-9]+(?:\\.[0-9]+)?:", "", hgvsc)
}

all_ch2 <- all_ch2 %>%
    dplyr::mutate(
        HGVSp_clean = dplyr::case_when(
            is.na(HGVSp) | HGVSp == "" | HGVSp == "." ~ NA_character_,
            TRUE ~ sub("^[A-Z]{2}_[0-9]+(?:\\.[0-9]+)?:", "", HGVSp)
        ),
        exon_label = dplyr::case_when(
            !is.na(ExonNumber) ~ paste0("exon", ExonNumber),
            TRUE               ~ "exonUnk"
        ),
        variant_id = dplyr::case_when(
            variant_category == "non_whitelist" ~ NA_character_,
            variant_category == "splice" ~ paste(
                Gene, MANE_Select_Str, exon_label,
                sapply(HGVSc, trim_splice_hgvsc), sep = ":"
            ),
            TRUE ~ paste(
                Gene, MANE_Select_Str, exon_label,
                sapply(HGVSc, clean_hgvsc),
                ifelse(!is.na(ProteinChange_1L), paste0("p.", ProteinChange_1L), ""),
                sep = ":"
            )
        ),
        variant_id = sub(":$", "", variant_id)
    ) %>%
    dplyr::select(-exon_label, -HGVSp_clean)

# won't match exactly, will have to grep the individual fields probably

# ========================
# SAVE
# ========================
write_xlsx(all_ch2, file.path("ch", "data", "ch_seq_wl_vars.xlsx"))

cat(sprintf("Final variants: %d\n", nrow(all_ch2)))
# 999
cat(sprintf("Final unique individuals: %d\n", length(unique(all_ch2$Sample.ID))))
# 710
cat(sprintf("Final unique genes: %d\n", length(unique(all_ch2$Gene))))
# 57

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

plot_batch_gene_freq(all_ch2, label = "post whitelist")


cov <- read.csv(file.path("ch", "data", "pmbb_brca12_cov_df.csv"), row.names = 1)
all_ids <- read.csv(file.path("ch", "data", "ch_psm_matched4_case_control_ids.csv"))$x
cov <- cov %>% filter(person_id %in% all_ids)
table(cov$Batch)
# 1    2
# 760 2244
dim(cov)
# 3004
controls <- cov %>% filter(BRCA12_Case == 0)
dim(controls)
# 2327

control_vars <- all_ch2 %>% filter(Sample.ID %in% controls$person_id)
plot_batch_gene_freq(control_vars, label = "post whitelist controls")

