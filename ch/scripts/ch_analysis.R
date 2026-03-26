library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
library(data.table, quietly=T)
library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr)
library(stringr)

# ========================
# prelim qc
# ========================
vars <- read.csv(file.path("ch", "data", "f3_ch.csv"))
dim(vars)
# 74230    65

length(unique(vars$Sample.ID))
# 1707
length(unique(vars$Gene))
# 72

found_ids <- sort(unique(vars$Sample.ID))

gList         <- fread(file.path("ch", "data", "whitelist_filter_20230531", "NEJM_2017_genes_01262020.txt"))
whitelist.mis <- fread(file.path("ch", "data", "whitelist_filter_20230531", "CHIP_missense_vars_cv_04102022.txt"))
whitelist.splice <- fread(file.path("ch", "data", "whitelist_filter_20230531", "CHIP_splice_vars_agb_01262020.txt"))
whitelist.LoF <- fread(file.path("ch", "data", "whitelist_filter_20230531", "CHIP_nonsense_FS_vars_agb_01262020.txt"))

# ========================
# Parse F1R2/F2R1 alt counts from "ref,alt" string format
# Do this BEFORE process_ch so the parsed columns survive into all steps
# ========================
parse_alt_count <- function(x) {
    sapply(x, function(v) {
        if (is.na(v) || v == ".") return(NA_real_)
        parts <- strsplit(as.character(v), ",")[[1]]
        if (length(parts) >= 2) as.numeric(parts[2]) else NA_real_
    })
}

vars$F1R2.alt <- parse_alt_count(vars$Sample.F1R2)
vars$F2R1.alt <- parse_alt_count(vars$Sample.F2R1)
vars$TLOD     <- as.numeric(vars$TLOD)

# ========================
# QC FILTER FUNCTION
# ========================
process_ch <- function(data, dataset_label, AF_LOWER, AF_UPPER, DP, AD) {

    datetime_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    log_file <- file.path("ch", "log", paste0(datetime_stamp, "_ch_qc_", dataset_label, ".txt"))
    dir.create(file.path("ch", "log"), showWarnings = FALSE, recursive = TRUE)

    data <- data %>%
        dplyr::mutate(
            pass_protein_coding = FALSE,
            pass_vaf            = FALSE,
            pass_depth          = FALSE,
            qc_pass             = FALSE
        )

    # Step 1: Protein coding filter
    data$pass_protein_coding <- (data$Bio.type == "protein_coding") &
        (data$Sample.AltFrac != ".") &
        (data$Sample.Depth   != ".") &
        (data$Sample.AltDepth != ".")

    # Convert to numeric
    data <- data %>%
        dplyr::mutate(
            Sample.AltFrac  = as.numeric(Sample.AltFrac),
            Sample.Depth    = as.numeric(Sample.Depth),
            Sample.AltDepth = as.numeric(Sample.AltDepth)
        )

    # Remove -1 sentinel values
    data$pass_protein_coding <- data$pass_protein_coding &
        !is.na(data$Sample.AltFrac)  &
        !is.na(data$Sample.Depth)    &
        !is.na(data$Sample.AltDepth) &
        (data$Sample.AltFrac  != -1) &
        (data$Sample.Depth    != -1) &
        (data$Sample.AltDepth != -1)

    # Step 2: VAF filter
    data$pass_vaf <- data$pass_protein_coding &
        (data$Sample.AltFrac >= AF_LOWER) &
        (data$Sample.AltFrac <= AF_UPPER)

    # Step 3: Depth filter
    data$pass_depth <- data$pass_vaf &
        (data$Sample.Depth    >= DP) &
        (data$Sample.AltDepth >= AD)

    data$qc_pass <- data$pass_depth

    # Log
    sink(log_file)
    cat("=== POST-VEP GENE QC LOG ===\n")
    cat("Dataset:", dataset_label, "\n\n")
    cat("=== FILTERING SUMMARY ===\n")
    cat("Initial variants:", nrow(data), "\n")
    cat("Pass protein coding:", sum(data$pass_protein_coding),
        sprintf("(%.1f%%)\n", 100 * sum(data$pass_protein_coding) / nrow(data)))
    cat("Pass VAF filter:", sum(data$pass_vaf),
        sprintf("(%.1f%%)\n", 100 * sum(data$pass_vaf) / nrow(data)))
    cat("Pass depth filter:", sum(data$pass_depth),
        sprintf("(%.1f%%)\n", 100 * sum(data$pass_depth) / nrow(data)))
    cat("\nFinal QC pass:", sum(data$qc_pass),
        sprintf("(%.1f%% of initial)\n\n", 100 * sum(data$qc_pass) / nrow(data)))
    sink()

    cat("Log saved to:", log_file, "\n")
    return(data)
}

vars <- process_ch(vars, "f3_mutect2_ch", 0.02, 0.45, 20, 3)

# ========================
# STEP 1 for QC plots = all variants before whitelist
# STEP 2 for QC plots = after QC pass
# ========================
step1 <- vars
step2 <- vars %>% filter(qc_pass)

# ========================
# DIAGNOSTIC PDF FUNCTION
# ========================
create_diagnostic_pdf <- function(step1, step2, step3, output_file,
                                  title = "CHIP Mutect2 - Filtering QC") {

    metrics <- list(
        list(col="Sample.Depth",    label="Total Depth",    xlab="DP",    color="steelblue"),
        list(col="Sample.AltDepth", label="Alt Depth",      xlab="AD",    color="darkgreen"),
        list(col="Sample.AltFrac",  label="VAF",            xlab="VAF",   color="purple"),
        list(col="TLOD",            label="TLOD",           xlab="TLOD",  color="firebrick"),
        list(col="F1R2.alt",        label="F1R2 Alt Reads", xlab="F1R2",  color="darkorange"),
        list(col="F2R1.alt",        label="F2R1 Alt Reads", xlab="F2R1",  color="darkgoldenrod")
    )

    steps <- list(
        list(data=step1, label="Initial"),
        list(data=step2, label="After QC (DP/VAF)"),
        list(data=step3, label="After Whitelist")
    )

    plots <- list()
    idx   <- 1

    for (s in steps) {
        n_samples <- length(unique(s$data$Sample.ID))
        n_vars    <- nrow(s$data)
        for (m in metrics) {
            if (m$col %in% names(s$data)) {
                plot_data <- s$data %>%
                    dplyr::select(val = !!sym(m$col)) %>%
                    filter(!is.na(val))

                plots[[idx]] <- ggplot(plot_data, aes(x = val)) +
                    geom_histogram(bins=80, fill=m$color, alpha=0.7,
                                   color="white", linewidth=0.1) +
                    labs(
                        title    = paste0(s$label, ": ", m$label),
                        subtitle = paste0("n=", n_vars, " variants | ", n_samples, " samples"),
                        x        = m$xlab,
                        y        = "Count"
                    ) +
                    theme_minimal() +
                    theme(
                        plot.title    = element_text(size=9, face="bold"),
                        plot.subtitle = element_text(size=7, color="gray40"),
                        axis.title    = element_text(size=8),
                        axis.text     = element_text(size=7)
                    )
                idx <- idx + 1
            }
        }
    }

    ncols <- length(metrics)   # 6 columns = one per metric
    nrows <- length(steps)     # 3 rows = one per filtering stage

    pdf(output_file, width=ncols * 3.5, height=nrows * 3.5)
    grid.arrange(
        grobs = plots,
        ncol  = ncols,
        top   = textGrob(title, gp=gpar(fontsize=14, fontface="bold"))
    )
    dev.off()
    cat("QC PDF saved to:", output_file, "\n")
}

# ========================
# SELECT COLUMNS FOR DOWNSTREAM ANALYSIS
# ========================
all_ch <- vars %>% filter(qc_pass) %>%
    select("Sample.ID", "Chr", "Start", "REF", "ALT", "Gene",
           "Variant.Class", "Variant.Consequence", "Bio.type",
           "HGVSc", "HGVSp", "EXON", "INTRON",
           "MANE_SELECT", "MANE_PLUS_CLINICAL",
           "Sample.Zyg", "Sample.Depth", "Sample.AltDepth", "Sample.AltFrac",
           "Sample.F1R2", "Sample.F2R1", "TLOD",
           "F1R2.alt", "F2R1.alt")   # keep parsed columns

# ========================
# MATCH TRANSCRIPTS
# ========================
all_ch <- all_ch %>%
    left_join(gList[, c("Gene", "Accession")], by="Gene")

# Strip version numbers for robust matching
# e.g. MANE_SELECT="NM_015338.3" vs Accession="NM_015338" -> both become "NM_015338"
all_ch <- all_ch %>%
    dplyr::mutate(
        correct_transcript = gsub("\\.\\d+$", "", MANE_SELECT) == gsub("\\.\\d+$", "", Accession)
    )

cat("Transcript mismatches:\n")
# all_ch %>%
#     filter(!correct_transcript) %>%
#     distinct(Gene, MANE_SELECT, Accession) %>%
#     print()
table(all_ch$correct_transcript)
# FALSE  TRUE
# 5734 12825

# ========================
# Extract protein change from HGVSp
# ========================
extractProteinChange <- function(hgvsp) {
    if (is.na(hgvsp) || hgvsp == "" || hgvsp == ".") return(NA)
    prot_change <- gsub("^p\\.", "", hgvsp)
    prot_change <- gsub("%3D", "=", prot_change)
    return(prot_change)
}

all_ch$ProteinChange <- sapply(all_ch$HGVSp, extractProteinChange)

# ========================
# Initialize whitelist flags
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
# 1) MISSENSE variants
# ========================
vmis <- grepl("missense", all_ch$Variant.Consequence)

convert_1to3 <- function(aa_change) {
    if (is.na(aa_change) || aa_change == "" || aa_change == ".") return(NA)

    aa_map <- c(
        "A"="Ala","C"="Cys","D"="Asp","E"="Glu","F"="Phe",
        "G"="Gly","H"="His","I"="Ile","K"="Lys","L"="Leu",
        "M"="Met","N"="Asn","P"="Pro","Q"="Gln","R"="Arg",
        "S"="Ser","T"="Thr","V"="Val","W"="Trp","Y"="Tyr",
        "*"="Ter","="="="
    )

    pattern <- "^([A-Z*])(\\d+)([A-Z*=])$"
    if (grepl(pattern, aa_change)) {
        parts  <- str_match(aa_change, pattern)
        ref_aa <- parts[2]
        position <- parts[3]
        alt_aa <- parts[4]
        ref_3  <- ifelse(ref_aa %in% names(aa_map), aa_map[ref_aa], ref_aa)
        alt_3  <- ifelse(alt_aa %in% names(aa_map), aa_map[alt_aa], alt_aa)
        return(paste0(ref_3, position, alt_3))
    } else {
        return(aa_change)
    }
}

whitelist.mis$AAChange_3letter <- sapply(whitelist.mis$AAChange, convert_1to3)

vmis_wl <- paste(all_ch$Gene, all_ch$ProteinChange, sep="_") %in%
    paste(whitelist.mis$Gene, whitelist.mis$AAChange_3letter, sep="_")

all_ch$whitelist[vmis & vmis_wl] <- TRUE
all_ch$wl.mis[vmis & vmis_wl]    <- TRUE

cat(sprintf("Missense variants in whitelist: %d\n", sum(vmis & vmis_wl)))

# ========================
# 2) LoF and frameshift variants
# ========================
vlof <- grepl("Ter|\\*|X", all_ch$ProteinChange) |
    grepl("stop_gain|stop_lost|start_lost", all_ch$Variant.Consequence) |
    grepl("fs", all_ch$ProteinChange, fixed=TRUE) |
    grepl("frameshift", all_ch$Variant.Consequence)

vLOFgene <- all_ch$Gene %in% whitelist.LoF$Gene

all_ch$whitelist[vlof & vLOFgene] <- TRUE
all_ch$wl.lof[vlof & vLOFgene]    <- TRUE

cat(sprintf("LoF/frameshift variants in whitelist: %d\n", sum(vlof & vLOFgene, na.rm=TRUE)))

# ========================
# 3) SPLICE variants — gate on correct transcript
# ========================
vSplice     <- grepl("splice", all_ch$Variant.Consequence)
vSplicegene <- all_ch$Gene %in% whitelist.splice$Gene
vSpliceCorrectTranscript <- all_ch$correct_transcript

# Only flag splice variants on the correct transcript
all_ch$whitelist[vSplice & vSplicegene & vSpliceCorrectTranscript] <- TRUE
all_ch$wl.splice[vSplice & vSplicegene & vSpliceCorrectTranscript] <- TRUE

# Flag wrong-transcript splice variants for manual review
all_ch$manualreview[(vSplice & vSplicegene) & (!vSpliceCorrectTranscript)] <- TRUE

cat(sprintf("Splice variants in whitelist (correct transcript only): %d\n",
            sum(vSplice & vSplicegene & vSpliceCorrectTranscript)))
cat(sprintf("Splice variants flagged for manual review (wrong transcript): %d\n",
            sum((vSplice & vSplicegene) & (!vSpliceCorrectTranscript))))

cat("\n=== WHITELIST FILTERING SUMMARY ===\n")
cat(sprintf("Variants in ANY whitelist: %d (%.1f%%)\n",
            sum(all_ch$whitelist), 100 * sum(all_ch$whitelist) / nrow(all_ch)))
cat(sprintf("  - Missense whitelist: %d\n", sum(all_ch$wl.mis)))
cat(sprintf("  - LoF whitelist: %d\n",      sum(all_ch$wl.lof)))
cat(sprintf("  - Splice whitelist: %d\n",   sum(all_ch$wl.splice)))
cat(sprintf("  - Flagged for manual review: %d\n", sum(all_ch$manualreview)))

# ========================
# GENE SPECIFIC RULES
# ========================
extract_exon_number <- function(exon_field) {
    if (is.na(exon_field) || exon_field == "" || exon_field == ".") return(NA)
    as.numeric(gsub("\\|.*", "", exon_field))
}
all_ch$ExonNumber <- sapply(all_ch$EXON, extract_exon_number)

extract_aa_position <- function(prot_change) {
    if (is.na(prot_change) || prot_change == "" || prot_change == ".") return(NA)
    as.numeric(gsub("\\D", "", prot_change))
}
all_ch$AAPosition <- sapply(all_ch$ProteinChange, extract_aa_position)

# Redefine variant type vectors with full criteria
vlof <- grepl("Ter|\\*|X", all_ch$ProteinChange) |
    grepl("stop_gain|stop_lost|start_lost", all_ch$Variant.Consequence) |
    grepl("fs", all_ch$ProteinChange, fixed=TRUE) |
    grepl("frameshift", all_ch$Variant.Consequence)
vFS     <- grepl("fs", all_ch$ProteinChange, fixed=TRUE) |
    grepl("frameshift_variant", all_ch$Variant.Consequence)
vSplice <- grepl("splice", all_ch$Variant.Consequence)
vSpliceCorrectTranscript <- all_ch$correct_transcript

# ========================
# 4) Gene-specific exceptions
# ========================

# ASXL1: Frameshift/nonsense/splice in exon 11-12
vexon11 <- all_ch$ExonNumber == 11
vexon12 <- all_ch$ExonNumber == 12

asxl1Exception <- (all_ch$Gene=="ASXL1") & (vlof | vFS) & (vexon11 | vexon12)
all_ch$whitelist[asxl1Exception]    <- TRUE
all_ch$wl.lof[asxl1Exception]       <- TRUE
all_ch$wl.exception[asxl1Exception] <- TRUE

asxl1ExceptionSplice <- (all_ch$Gene=="ASXL1") & vSplice & vSpliceCorrectTranscript & (vexon11 | vexon12)
all_ch$whitelist[asxl1ExceptionSplice]    <- TRUE
all_ch$wl.splice[asxl1ExceptionSplice]    <- TRUE
all_ch$wl.exception[asxl1ExceptionSplice] <- TRUE

cat(sprintf("ASXL1 exceptions: %d\n", sum(asxl1Exception | asxl1ExceptionSplice, na.rm=TRUE)))

# ASXL2: Frameshift/nonsense/splice in exon 11-12
asxl2Exception <- (all_ch$Gene=="ASXL2") & (vlof | vFS) & (vexon11 | vexon12)
all_ch$whitelist[asxl2Exception]    <- TRUE
all_ch$wl.lof[asxl2Exception]       <- TRUE
all_ch$wl.exception[asxl2Exception] <- TRUE

asxl2ExceptionSplice <- (all_ch$Gene=="ASXL2") & vSplice & vSpliceCorrectTranscript & (vexon11 | vexon12)
all_ch$whitelist[asxl2ExceptionSplice]    <- TRUE
all_ch$wl.splice[asxl2ExceptionSplice]    <- TRUE
all_ch$wl.exception[asxl2ExceptionSplice] <- TRUE

cat(sprintf("ASXL2 exceptions: %d\n", sum(asxl2Exception | asxl2ExceptionSplice, na.rm=TRUE)))

# PPM1D: Frameshift/nonsense in exon 5 or 6
vexon5 <- all_ch$ExonNumber == 5
vexon6 <- all_ch$ExonNumber == 6

ppm1dException <- (all_ch$Gene=="PPM1D") & (vlof | vFS) & (vexon5 | vexon6)
all_ch$whitelist[ppm1dException]    <- TRUE
all_ch$wl.lof[ppm1dException]       <- TRUE
all_ch$wl.exception[ppm1dException] <- TRUE

cat(sprintf("PPM1D exceptions: %d\n", sum(ppm1dException, na.rm=TRUE)))

# TET2: Missense in catalytic domains (p.1104-1481 and p.1843-2002)
vmis    <- grepl("missense", all_ch$Variant.Consequence)
TETidx  <- which(all_ch$Gene=="TET2" & vmis & !is.na(all_ch$AAPosition))

for (i in TETidx) {
    AApos <- all_ch$AAPosition[i]
    if (!is.na(AApos) && ((AApos>=1104 & AApos<=1481) | (AApos>=1843 & AApos<=2002))) {
        all_ch$whitelist[i]    <- TRUE
        all_ch$wl.mis[i]       <- TRUE
        all_ch$wl.exception[i] <- TRUE
    }
}
cat(sprintf("TET2 catalytic domain exceptions: %d\n",
            sum(all_ch$wl.exception & all_ch$Gene=="TET2", na.rm=TRUE)))

# CBL: RING finger missense p.381-421
CBLidx <- which(all_ch$Gene=="CBL" & vmis & !is.na(all_ch$AAPosition))
for (i in CBLidx) {
    AApos <- all_ch$AAPosition[i]
    if (!is.na(AApos) && AApos>=381 & AApos<=421) {
        all_ch$whitelist[i]    <- TRUE
        all_ch$wl.mis[i]       <- TRUE
        all_ch$wl.exception[i] <- TRUE
    }
}
cat(sprintf("CBL RING finger exceptions: %d\n",
            sum(all_ch$wl.exception & all_ch$Gene=="CBL", na.rm=TRUE)))

# CBLB: RING finger missense p.372-412
CBLBidx <- which(all_ch$Gene=="CBLB" & vmis & !is.na(all_ch$AAPosition))
for (i in CBLBidx) {
    AApos <- all_ch$AAPosition[i]
    if (!is.na(AApos) && AApos>=372 & AApos<=412) {
        all_ch$whitelist[i]    <- TRUE
        all_ch$wl.mis[i]       <- TRUE
        all_ch$wl.exception[i] <- TRUE
    }
}
cat(sprintf("CBLB RING finger exceptions: %d\n",
            sum(all_ch$wl.exception & all_ch$Gene=="CBLB", na.rm=TRUE)))

# ========================
# 5) Manual review flags
# ========================
vNonFrameshiftDel <- grepl("inframe_deletion",  all_ch$Variant.Consequence)
vNonFrameshiftIns <- grepl("inframe_insertion", all_ch$Variant.Consequence)
vFrameshiftIndel  <- grepl("frameshift_variant", all_ch$Variant.Consequence)

all_ch$manualreview[(vlof | vFS | vSplice) & (all_ch$Gene=="GATA3")]  <- TRUE
all_ch$manualreview[vNonFrameshiftDel       & (all_ch$Gene=="CREBBP")] <- TRUE
all_ch$manualreview[(vlof | vFS | vSplice)  & (all_ch$Gene=="CSF3R")]  <- TRUE
all_ch$manualreview[vNonFrameshiftDel       & (all_ch$Gene=="DNMT3A")] <- TRUE
all_ch$manualreview[vNonFrameshiftDel       & (all_ch$Gene=="EP300")]  <- TRUE
all_ch$manualreview[(vNonFrameshiftDel | vNonFrameshiftIns) & (all_ch$Gene=="FLT3")]  <- TRUE
all_ch$manualreview[(vNonFrameshiftDel | vNonFrameshiftIns) & (all_ch$Gene=="JAK2")]  <- TRUE
all_ch$manualreview[vNonFrameshiftDel       & (all_ch$Gene=="KDM6A")] <- TRUE
all_ch$manualreview[(vNonFrameshiftDel | vNonFrameshiftIns) & (all_ch$Gene=="KIT")]   <- TRUE
all_ch$manualreview[(vNonFrameshiftDel | vNonFrameshiftIns) & (all_ch$Gene=="MPL")]   <- TRUE
all_ch$manualreview[vFrameshiftIndel        & (all_ch$Gene=="NPM1")]  <- TRUE

cat("\n=== MANUAL REVIEW FLAGS ===\n")
cat(sprintf("Total flagged for manual review: %d\n", sum(all_ch$manualreview, na.rm=TRUE)))

# ========================
# Final summary
# ========================
cat("\n=== FINAL WHITELIST + EXCEPTION SUMMARY ===\n")
cat(sprintf("Total variants: %d\n", nrow(all_ch)))
cat(sprintf("Variants in whitelist (including exceptions): %d (%.1f%%)\n",
            sum(all_ch$whitelist), 100 * sum(all_ch$whitelist) / nrow(all_ch)))
cat(sprintf("  - Missense whitelist: %d\n",  sum(all_ch$wl.mis)))
cat(sprintf("  - LoF whitelist: %d\n",       sum(all_ch$wl.lof)))
cat(sprintf("  - Splice whitelist: %d\n",    sum(all_ch$wl.splice)))
cat(sprintf("  - Exception rules: %d\n",     sum(all_ch$wl.exception)))
cat(sprintf("  - Flagged for manual review: %d\n", sum(all_ch$manualreview)))

# ========================
# STEP 3 for QC plots = after whitelist
# ========================
step3 <- all_ch %>% filter(whitelist)

# ========================
# GENERATE QC PDF
# ========================
datetime_stamp <- format(Sys.time(), "%Y%m%d")
dir.create(file.path("ch", "output"), showWarnings=FALSE, recursive=TRUE)

create_diagnostic_pdf(
    step1, step2, step3,
    output_file = file.path("ch", "output", paste0(datetime_stamp, "_f3_mutect2_qc.pdf")),
    title       = "CHIP Mutect2 - Filtering QC"
)

# ========================
# Write output files
# ========================
check_vars <- data.frame(
    total_num_variants    = nrow(all_ch),
    total_num_whitelist   = sum(all_ch$whitelist),
    total_num_manualreview = sum(all_ch$manualreview)
)

dim(all_ch)
dim(all_ch[all_ch$whitelist, ])
dim(all_ch[all_ch$manualreview, ])

write.csv(check_vars,
          file.path("ch", "output", paste0(datetime_stamp, "_f3_mutect2_varcount.csv")),
          row.names=FALSE)

write.csv(all_ch,
          file.path("ch", "output", paste0(datetime_stamp, "_f3_mutect2_allvariants.csv")),
          row.names=FALSE)

write.csv(all_ch[all_ch$whitelist, ],
          file.path("ch", "output", paste0(datetime_stamp, "_f3_mutect2_wl.csv")),
          row.names=FALSE)

write.csv(all_ch[all_ch$manualreview, ],
          file.path("ch", "output", paste0(datetime_stamp, "_f3_mutect2_manualreview.csv")),
          row.names=FALSE)

# ========================
# PRELIM PLOTS
# ========================
# ========================
# WHITELIST BREAKDOWN FIGURE
# ========================
library(ggplot2)
library(dplyr)
library(scales)

# Build summary data
wl_breakdown <- data.frame(
    Category = c(
        "Splice",
        "LoF / Frameshift",
        "Missense",
        "Gene-specific\nexceptions",
        "Manual review"
    ),
    N = c(
        sum(all_ch$wl.splice & !all_ch$wl.exception),
        sum(all_ch$wl.lof & !all_ch$wl.exception),
        sum(all_ch$wl.mis & !all_ch$wl.exception),
        sum(all_ch$wl.exception),
        sum(all_ch$manualreview)
    ),
    Group = c("Whitelisted", "Whitelisted", "Whitelisted", "Whitelisted", "Manual review")
)

wl_breakdown$Category <- factor(wl_breakdown$Category,
                                levels = wl_breakdown$Category[order(wl_breakdown$N)])

colors <- c(
    "Splice"                   = "#1D9E75",
    "LoF / Frameshift"         = "#378ADD",
    "Missense"                 = "#D4537E",
    "Gene-specific\nexceptions"= "#BA7517",
    "Manual review"            = "#888780"
)

p1 <- ggplot(wl_breakdown, aes(x = Category, y = N, fill = Category)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = paste0(N, "\n(", round(100 * N / sum(N), 1), "%)")),
              hjust = -0.1, size = 3.2, color = "gray30") +
    coord_flip() +
    scale_fill_manual(values = colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25)),
                       labels = comma) +
    labs(
        title    = "CHIP variant whitelist breakdown",
        subtitle = paste0("Total whitelisted: ", sum(all_ch$whitelist),
                          " | Manual review: ", sum(all_ch$manualreview),
                          " | QC-pass input: ", nrow(all_ch)),
        x        = NULL,
        y        = "Number of variants"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        legend.position  = "none",
        plot.title       = element_text(face = "bold", size = 13),
        plot.subtitle    = element_text(size = 10, color = "gray50"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text.y      = element_text(size = 11)
    )

# ========================
# GENE-SPECIFIC EXCEPTION BREAKDOWN
# ========================
gene_exceptions <- all_ch %>%
    filter(wl.exception) %>%
    count(Gene, name = "N") %>%
    arrange(desc(N)) %>%
    mutate(Gene = factor(Gene, levels = Gene[order(N)]))

p2 <- ggplot(gene_exceptions, aes(x = Gene, y = N)) +
    geom_col(width = 0.6, fill = "#BA7517") +
    geom_text(aes(label = N), hjust = -0.2, size = 3.2, color = "gray30") +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(
        title = "Gene-specific exceptions",
        x     = NULL,
        y     = "Number of variants"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        plot.title         = element_text(face = "bold", size = 12),
        axis.text.y        = element_text(size = 11)
    )

# ========================
# COMBINE AND SAVE
# ========================
library(patchwork)

final_plot <- p1 / p2 + plot_layout(heights = c(2, 1))

datetime_stamp <- format(Sys.time(), "%Y%m%d")
ggsave(
    file.path("ch", "output", paste0(datetime_stamp, "_f3_mutect2_whitelist_breakdown.pdf")),
    final_plot,
    width  = 8,
    height = 9
)

ggsave(
    file.path("ch", "output", paste0(datetime_stamp, "_f3_mutect2_whitelist_breakdown.png")),
    final_plot,
    width  = 8,
    height = 9,
    dpi    = 300
)

cat("Whitelist breakdown figure saved.\n")

