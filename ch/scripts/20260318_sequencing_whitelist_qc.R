library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
library(data.table, quietly=T)
library(ggplot2)
library(patchwork)

# ========================
# READ IN WHITELISTS
# ========================
gList<-fread(file.path("ch", "data", "whitelist_filter_20230531", "NEJM_2017_genes_01262020.txt"))
whitelist.mis<-fread(file.path("ch", "data", "whitelist_filter_20230531", "CHIP_missense_vars_cv_04102022.txt"))
whitelist.splice<-fread(file.path("ch", "data", "whitelist_filter_20230531", "CHIP_splice_vars_agb_01262020.txt"))
whitelist.LoF<-fread(file.path("ch", "data", "whitelist_filter_20230531", "CHIP_nonsense_FS_vars_agb_01262020.txt"))

# ========================
# READ IN SOMATIC CALLS
# ========================
vars3 <- read.csv(file.path("ch", "data", "f3_ch.csv"))
dim(vars3)
# 98654    65
length(unique(vars3$Sample.ID))
# 2257
length(unique(vars3$Gene))
# 72

vars2 <- read.csv(file.path("ch", "data", "f2_ch.csv"))
dim(vars2)
# 34067    65
length(unique(vars2$Sample.ID))
# 759
length(unique(vars2$Gene))
# 72

# combine into one mega df with Batch col with 1 if it came from f2 and 2 if it came from f3
vars2$Batch <- 1
vars3$Batch <- 2

vars <- rbind(vars2, vars3)
write.csv(vars, file.path("ch", "data", "ch_vars.csv"))

cat("Unique samples:", length(unique(vars$Sample.ID)), "\n")
# 3003
cat("Unique genes:", length(unique(vars$Gene)), "\n")
# 72

genes_found <- unique(vars$Gene)
length(genes_found)
gList$Gene[!(gList$Gene %in% genes_found)]
# "MLL"  "MLL2"
# didn't find these...

# ========================
# PARSE F1R2/F2R1
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
vars$TLOD <- as.numeric(vars$TLOD)

# histogram of vars alt
p1 <- ggplot(vars, aes(x = F1R2.alt)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.8) +
    labs(title = "F1R2 Alt Read Count", x = "Alt Reads (F1R2)", y = "Count") +
    theme_minimal()

p2 <- ggplot(vars, aes(x = F2R1.alt)) +
    geom_histogram(bins = 50, fill = "darkgreen", color = "white", alpha = 0.8) +
    labs(title = "F2R1 Alt Read Count", x = "Alt Reads (F2R1)", y = "Count") +
    theme_minimal()

combined <- p1 + p2 + plot_annotation(title = "Strand Orientation Alt Read Distributions")

ggsave(file.path("ch", "figures", "ch_f1r2_f2r1_alt_hist.png"), plot = combined, width = 8, height = 5, dpi = 300)

# ========================
# SEQUENCING QC
# ========================
process_ch <- function(data, dataset_label, AF_LOWER, AF_UPPER, DP, AD) {
    # Generate datetime stamp and log file
    datetime_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    log_file <- file.path("ch", "log", paste0(datetime_stamp, "_ch_qc_", dataset_label, ".txt"))
    summary_file <- file.path("ch", "log", paste0(datetime_stamp, "_ch_qc_", dataset_label, "_summary.csv"))
    dir.create(file.path("ch", "log"), showWarnings = FALSE, recursive = TRUE)

    # Initialize flags
    data <- data %>%
        dplyr::mutate(
            pass_protein_coding = FALSE,
            pass_vaf = FALSE,
            pass_depth = FALSE,
            pass_orientation = FALSE,
            qc_pass = FALSE
        )

    # Step 1: Protein coding filter
    data$pass_protein_coding <- (data$Bio.type == "protein_coding") &
        (data$Sample.AltFrac != ".") &
        (data$Sample.Depth != ".") &
        (data$Sample.AltDepth != ".")

    # Convert to numeric
    data <- data %>%
        dplyr::mutate(
            Sample.AltFrac = as.numeric(Sample.AltFrac),
            Sample.Depth = as.numeric(Sample.Depth),
            Sample.AltDepth = as.numeric(Sample.AltDepth)
        )

    # Remove -1 values from passing
    data$pass_protein_coding <- data$pass_protein_coding &
        !is.na(data$Sample.AltFrac) &
        !is.na(data$Sample.Depth) &
        !is.na(data$Sample.AltDepth) &
        (data$Sample.AltFrac != -1) &
        (data$Sample.Depth != -1) &
        (data$Sample.AltDepth != -1)

    # Step 2: VAF filter
    data$pass_vaf <- data$pass_protein_coding &
        (data$Sample.AltFrac >= AF_LOWER) &
        (data$Sample.AltFrac <= AF_UPPER)

    # Step 3: Depth filter
    data$pass_depth <- data$pass_vaf &
        (data$Sample.Depth >= DP) &
        (data$Sample.AltDepth >= AD)

    # Step 4: Strand orientation filter
    data$pass_orientation <- data$pass_depth &
        !is.na(data$F1R2.alt) &
        !is.na(data$F2R1.alt) &
        (data$F1R2.alt >= 1) &
        (data$F2R1.alt >= 1)

    # Final QC pass
    data$qc_pass <- data$pass_orientation

    # Build summary table
    n_initial <- nrow(data)
    summary_df <- data.frame(
        Step = c(
            "Initial",
            "Pass protein coding",
            "Pass VAF filter",
            "Pass depth filter",
            "Pass orientation filter"
        ),
        Variants_Remaining = c(
            n_initial,
            sum(data$pass_protein_coding),
            sum(data$pass_vaf),
            sum(data$pass_depth),
            sum(data$pass_orientation)
        )
    )
    summary_df$Percent_of_Original <- round(
        100 * summary_df$Variants_Remaining / n_initial, 1
    )

    # Log summary
    sink(log_file)
    cat("=== POST-VEP GENE QC LOG ===\n")
    cat("Dataset:", dataset_label, "\n\n")
    cat("=== FILTERING SUMMARY ===\n")
    print(summary_df, row.names = FALSE)
    cat("\n")
    sink()
    cat("Log saved to:", log_file, "\n")

    # Save summary CSV
    write.csv(summary_df, summary_file, row.names = FALSE)
    cat("Summary CSV saved to:", summary_file, "\n")

    return(data)
}

vars <- process_ch(vars, "ch", 0.02, 0.45, 20, 3)

### PLOT DIAGNOSTICS
create_diagnostic_pdf <- function(step1, step2, step3, gene_name, output_file) {
    plots <- list()

    metrics <- list(
        list(col = "Sample.Depth", name = "Total Depth", xlab = "Total Depth (DP)",
             color = "steelblue"),
        list(col = "Sample.AltDepth", name = "Alt Depth", xlab = "Alt Depth (AD)",
             color = "darkgreen"),
        list(col = "Sample.AltFrac", name = "VAF", xlab = "Variant Allele Fraction",
             color = "purple"),
        list(col = "F1R2.alt", name = "F1R2 Alt", xlab = "Alt Reads (F1R2)",
             color = "steelblue"),
        list(col = "F2R1.alt", name = "F2R1 Alt", xlab = "Alt Reads (F2R1)",
             color = "darkgreen")
    )

    steps <- list(
        list(data = step1, label = "Initial"),
        list(data = step2, label = "After Depth/VAF"),
        list(data = step3, label = "Final")
    )

    idx <- 1
    for (s in 1:3) {
        step_data <- steps[[s]]$data
        step_label <- steps[[s]]$label
        n <- length(unique(step_data$Sample.ID))

        for (m in metrics) {
            if (m$col %in% names(step_data)) {
                plots[[idx]] <- ggplot(step_data, aes(x = .data[[m$col]])) +
                    geom_histogram(bins = 100, fill = m$color, alpha = 0.7, color = "white") +
                    labs(title = paste0(step_label, ": ", m$name),
                         subtitle = paste("n =", n, "individuals"),
                         x = m$xlab, y = "Count") +
                    theme_minimal() +
                    theme(plot.title = element_text(size = 10),
                          plot.subtitle = element_text(size = 8))
                idx <- idx + 1
            }
        }
    }

    pdf(output_file, width = 20, height = 12)
    grid.arrange(
        grobs = plots,
        ncol = 5,
        top = textGrob(paste(gene_name, "Filtering Pipeline"),
                       gp = gpar(fontsize = 16, fontface = "bold"))
    )
    dev.off()
}

step1 <- vars                            # all variants, pre-filtering
step2 <- vars %>% filter(pass_depth)     # after VAF + depth
step3 <- vars %>% filter(pass_orientation) # final / orientation pass

create_diagnostic_pdf(
    step1 = step1,
    step2 = step2,
    step3 = step3,
    gene_name = "CHIP",
    output_file = file.path("ch", "figures", "ch_filtering_diagnostic.pdf")
)

all_ch <- vars %>% filter(qc_pass)

dim(all_ch)
# 18312    26

write.csv(all_ch, file.path("ch", "data", "ch_seq_vars.csv"))

# ========================
# MATCH TRANSCRIPTS
# ========================
# convert HTML to character coding
all_ch <- all_ch %>% mutate(across(where(is.character), ~ URLdecode(.x)))

# ensure that none are in last exon?

# join gene list accessions
all_ch <- all_ch %>%
    left_join(gList[, c("Gene", "Accession")], by = "Gene")

all_ch <- all_ch %>%
    dplyr::mutate(
        MANE_SELECT_stripped = gsub("\\.\\d+$", "", MANE_SELECT),
        Accession_stripped   = gsub("\\.\\d+$", "", Accession),
        correct_transcript   = MANE_SELECT_stripped == Accession_stripped
    )

# Check mismatches
all_ch %>%
    filter(!correct_transcript) %>%
    distinct(Gene, MANE_SELECT, Accession) %>%
    print()

table(all_ch$correct_transcript)

# FALSE  TRUE
# 4116 14196

all_ch <- all_ch %>% filter(correct_transcript == TRUE)
dim(all_ch)
# 14196    26

# ========================
# Extract protein change from HGVSp
# ========================
# everything after "p."
extractProteinChange <- function(hgvsp) {
    if (is.na(hgvsp) || hgvsp == "" || hgvsp == ".") {
        return(NA)
    }
    # Remove "p." prefix and decode URL encoding
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
        whitelist = FALSE,
        wl.mis = FALSE,
        wl.lof = FALSE,
        wl.splice = FALSE,
        wl.exception = FALSE,
        manualreview = FALSE
    )

# ========================
# PREP COLUMNS
# ========================
### EXONS ###
extract_exon_info <- function(exon_field) {
    if (is.na(exon_field) || exon_field == "" || exon_field == ".") {
        return(list(exon_num = NA_real_, exon_total = NA_real_))
    }
    parts <- strsplit(as.character(exon_field), "\\|")[[1]]
    list(
        exon_num   = as.numeric(parts[1]),
        exon_total = as.numeric(parts[2])
    )
}

all_ch <- all_ch %>%
    mutate(
        ExonNumber = sapply(EXON, function(x) extract_exon_info(x)$exon_num),
        ExonTotal  = sapply(EXON, function(x) extract_exon_info(x)$exon_total),
        is_last_exon = !is.na(ExonNumber) & !is.na(ExonTotal) & ExonNumber == ExonTotal
    )

# sanity check
cat("Variants in last exon:", sum(all_ch$is_last_exon, na.rm = TRUE), "\n")
# 4535
# print(all_ch %>% filter(is_last_exon) %>% count(Gene, ExonTotal) %>% arrange(Gene))

# filter
all_ch <- all_ch %>% filter(!is_last_exon)
cat("Variants after last exon filter:", nrow(all_ch), "\n")
# 9661

### EXONS ###
# Extract amino acid position from ProteinChange for range-based filters
extract_aa_position <- function(prot_change) {
    if (is.na(prot_change) || prot_change == "" || prot_change == ".") {
        return(NA)
    }
    # Extract digits from protein change (e.g., "Arg573Gln" -> 573)
    position <- as.numeric(gsub("\\D", "", prot_change))
    return(position)
}

all_ch$AAPosition <- sapply(all_ch$ProteinChange, extract_aa_position)

# ========================
# 1) Handle MISSENSE variants
# ========================
vmis <- grepl("missense", all_ch$Variant.Consequence)

convert_1to3 <- function(aa_change) {
    if (is.na(aa_change) || aa_change == "" || aa_change == ".") {
        return(NA)
    }

    aa_map <- c(
        "A" = "Ala", "C" = "Cys", "D" = "Asp", "E" = "Glu", "F" = "Phe",
        "G" = "Gly", "H" = "His", "I" = "Ile", "K" = "Lys", "L" = "Leu",
        "M" = "Met", "N" = "Asn", "P" = "Pro", "Q" = "Gln", "R" = "Arg",
        "S" = "Ser", "T" = "Thr", "V" = "Val", "W" = "Trp", "Y" = "Tyr",
        "*" = "Ter", "=" = "="
    )

    # Pattern: one letter + digits + one letter/symbol
    pattern <- "^([A-Z*])(\\d+)([A-Z*=])$"

    if (grepl(pattern, aa_change)) {
        parts <- str_match(aa_change, pattern)
        ref_aa <- parts[2]
        position <- parts[3]
        alt_aa <- parts[4]

        # Convert
        ref_3 <- ifelse(ref_aa %in% names(aa_map), aa_map[ref_aa], ref_aa)
        alt_3 <- ifelse(alt_aa %in% names(aa_map), aa_map[alt_aa], alt_aa)

        return(paste0(ref_3, position, alt_3))
    } else {
        return(aa_change)
    }
}

# convert whitelist to 3 letter codes
whitelist.mis$AAChange_3letter <- sapply(whitelist.mis$AAChange, convert_1to3)

vmis_wl <- paste(all_ch$Gene, all_ch$ProteinChange, sep = "_") %in%
    paste(whitelist.mis$Gene, whitelist.mis$AAChange_3letter, sep = "_")

all_ch$whitelist[vmis & vmis_wl] <- TRUE
all_ch$wl.mis[vmis & vmis_wl] <- TRUE

cat(sprintf("Missense variants in whitelist: %d\n", sum(vmis & vmis_wl)))
# 55

# examine
all_ch %>% filter(wl.mis == TRUE)

# ========================
# 2) Handle LoF and frameshift variants
# ========================
# Stop gain/loss/start lost OR frameshift
vlof <- grepl("Ter|\\*|X", all_ch$ProteinChange) |
    grepl("stop_gain|stop_lost|start_lost", all_ch$Variant.Consequence) |
    grepl("fs", all_ch$ProteinChange, fixed = TRUE) |
    grepl("frameshift", all_ch$Variant.Consequence)

vLOFgene <- all_ch$Gene %in% whitelist.LoF$Gene

all_ch$whitelist[vlof & vLOFgene] <- TRUE
all_ch$wl.lof[vlof & vLOFgene] <- TRUE

cat(sprintf("LoF/frameshift variants in whitelist: %d\n", sum(vlof & vLOFgene, na.rm = TRUE)))

all_ch %>% filter(wl.lof == TRUE)
# 250

# ========================
# 3) Handle SPLICE variants
# ========================
vSplice <- grepl("splice", all_ch$Variant.Consequence)

vSplicegene <- all_ch$Gene %in% whitelist.splice$Gene

# check if correct transcript
vSpliceCorrectTranscript <- all_ch$correct_transcript

# Flag splice variants
all_ch$whitelist[vSplice & vSplicegene & vSpliceCorrectTranscript] <- TRUE
all_ch$wl.splice[vSplice & vSplicegene & vSpliceCorrectTranscript] <- TRUE

# # Flag for manual review if splice gene but wrong transcript (shouldn't happen after canonical filter)
all_ch$manualreview[(vSplice & vSplicegene) & (!vSpliceCorrectTranscript)] <- TRUE

cat(sprintf("Splice variants in whitelist: %d\n",
            sum(vSplice & vSplicegene)))
# 2236

# ========================
# GENE SPECIFIC RULES
# ========================
### Define LoF and frameshift variants ###
vlof <- grepl("Ter|\\*|X", all_ch$ProteinChange) |
    grepl("stop_gain|stop_lost|start_lost", all_ch$Variant.Consequence) |
    grepl("fs", all_ch$ProteinChange, fixed = TRUE) |
    grepl("frameshift", all_ch$Variant.Consequence)
vFS <- grepl("fs", all_ch$ProteinChange, fixed = TRUE) |
    grepl("frameshift_variant", all_ch$Variant.Consequence)
vSplice <- grepl("splice", all_ch$Variant.Consequence)
vSpliceCorrectTranscript <- all_ch$correct_transcript

#### ASXL1 ###
# Frameshift/nonsense/splice-site in exon 11-12
vexon11 <- all_ch$ExonNumber == 11
vexon12 <- all_ch$ExonNumber == 12

asxl1Exception <- (all_ch$Gene == "ASXL1") & (vlof | vFS) & (vexon11 | vexon12)

all_ch$whitelist[asxl1Exception] <- TRUE
all_ch$wl.lof[asxl1Exception] <- TRUE
all_ch$wl.exception[asxl1Exception] <- TRUE

# Splicing
asxl1ExceptionSplice <- (all_ch$Gene == "ASXL1") &
    vSplice & vSpliceCorrectTranscript &
    (vexon11 | vexon12)

all_ch$whitelist[asxl1ExceptionSplice] <- TRUE
all_ch$wl.splice[asxl1ExceptionSplice] <- TRUE
all_ch$wl.exception[asxl1ExceptionSplice] <- TRUE

cat(sprintf("ASXL1 exceptions: %d\n", sum(asxl1Exception | asxl1ExceptionSplice, na.rm = TRUE)))
# 4

### ASXL2: Frameshift/nonsense/splice-site in exon 11-12
# LOF
asxl2Exception <- (all_ch$Gene == "ASXL2") & (vlof | vFS) & (vexon11 | vexon12)
all_ch$whitelist[asxl2Exception] <- TRUE
all_ch$wl.lof[asxl2Exception] <- TRUE
all_ch$wl.exception[asxl2Exception] <- TRUE

# Splicing
asxl2ExceptionSplice <- (all_ch$Gene == "ASXL2") &
    vSplice & vSpliceCorrectTranscript &
    (vexon11 | vexon12)

all_ch$whitelist[asxl2ExceptionSplice] <- TRUE
all_ch$wl.splice[asxl2ExceptionSplice] <- TRUE
all_ch$wl.exception[asxl2ExceptionSplice] <- TRUE

cat(sprintf("ASXL2 exceptions: %d\n", sum(asxl2Exception | asxl2ExceptionSplice, na.rm = TRUE)))
# 8

#### PPM1D: Frameshift/nonsense in exon 5 or 6
vexon5 <- all_ch$ExonNumber == 5
vexon6 <- all_ch$ExonNumber == 6

ppm1dException <- (all_ch$Gene == "PPM1D") & (vlof | vFS) & (vexon5 | vexon6)
all_ch$whitelist[ppm1dException] <- TRUE
all_ch$wl.lof[ppm1dException] <- TRUE
all_ch$wl.exception[ppm1dException] <- TRUE

cat(sprintf("PPM1D exceptions: %d\n", sum(ppm1dException, na.rm = TRUE)))
# 0

### TET2: Missense mutations in catalytic domains (p.1104-1481 and 1843-2002)
vmis <- grepl("missense", all_ch$Variant.Consequence)
TETidx <- which(all_ch$Gene == "TET2" &  vmis & !is.na(all_ch$AAPosition))

for(i in TETidx) {
    AApos <- all_ch$AAPosition[i]
    if(!is.na(AApos) && ((AApos >= 1104 & AApos <= 1481) | (AApos >= 1843 & AApos <= 2002))) {
        all_ch$whitelist[i] <- TRUE
        all_ch$wl.mis[i] <- TRUE
        all_ch$wl.exception[i] <- TRUE
    }
}

cat(sprintf("TET2 catalytic domain exceptions: %d\n", sum(all_ch$wl.exception & all_ch$Gene == "TET2", na.rm = TRUE)))
# 22

# CBL: RING finger missense p.381-421
CBLidx <- which(all_ch$Gene == "CBL" & vmis & !is.na(all_ch$AAPosition))
for(i in CBLidx) {
    AApos <- all_ch$AAPosition[i]
    if(!is.na(AApos) && AApos >= 381 & AApos <= 421) {
        all_ch$whitelist[i] <- TRUE
        all_ch$wl.mis[i] <- TRUE
        all_ch$wl.exception[i] <- TRUE
    }
}

cat(sprintf("CBL RING finger exceptions: %d\n", sum(all_ch$wl.exception & all_ch$Gene == "CBL", na.rm = TRUE)))
# 1

# CBLB: RING finger missense p.372-412
CBLBidx <- which(all_ch$Gene == "CBLB" & vmis & !is.na(all_ch$AAPosition))

for(i in CBLBidx) {
    AApos <- all_ch$AAPosition[i]
    if(!is.na(AApos) && AApos >= 372 & AApos <= 412) {
        all_ch$whitelist[i] <- TRUE
        all_ch$wl.mis[i] <- TRUE
        all_ch$wl.exception[i] <- TRUE
    }
}

cat(sprintf("CBLB RING finger exceptions: %d\n",
            sum(all_ch$wl.exception & all_ch$Gene == "CBLB", na.rm = TRUE)))
# 1

# ========================
# 5) Flag remaining exceptions for manual review
# ========================
# Non-frameshift indels (inframe insertions/deletions)
vNonFrameshiftDel <- grepl("inframe_deletion", all_ch$Variant.Consequence)
vNonFrameshiftIns <- grepl("inframe_insertion", all_ch$Variant.Consequence)

# GATA3: Frameshift/nonsense/splice-site ZNF domain
all_ch$manualreview[(vlof | vFS | vSplice) & (all_ch$Gene == "GATA3")] <- TRUE

# CREBBP: S1680del
all_ch$manualreview[vNonFrameshiftDel & (all_ch$Gene == "CREBBP")] <- TRUE

# CSF3R: truncating c.741-791
all_ch$manualreview[(vlof | vFS | vSplice) & (all_ch$Gene == "CSF3R")] <- TRUE

# DNMT3A: F732del, F752del
all_ch$manualreview[vNonFrameshiftDel & (all_ch$Gene == "DNMT3A")] <- TRUE

# EP300: VF1148_1149del
all_ch$manualreview[vNonFrameshiftDel & (all_ch$Gene == "EP300")] <- TRUE

# FLT3: FY590-591GD, del835
all_ch$manualreview[(vNonFrameshiftDel | vNonFrameshiftIns) & (all_ch$Gene == "FLT3")] <- TRUE

# JAK2: del/ins537-539L, del/ins538-539L, del/ins540-543MK, del/ins540-544MK, del542-543, del543-544, ins11546-547
all_ch$manualreview[(vNonFrameshiftDel | vNonFrameshiftIns) & (all_ch$Gene == "JAK2")] <- TRUE

# KDM6A: del419
all_ch$manualreview[vNonFrameshiftDel & (all_ch$Gene == "KDM6A")] <- TRUE

# KIT: ins503, del560, del579, del551-559
all_ch$manualreview[(vNonFrameshiftDel | vNonFrameshiftIns) & (all_ch$Gene == "KIT")] <- TRUE

# MPL: del513 W515-518KT
all_ch$manualreview[(vNonFrameshiftDel | vNonFrameshiftIns) & (all_ch$Gene == "MPL")] <- TRUE

# NPM1: Frameshift p.W288fs (insertion at c.859_860, 860_861, 862_863, 863_864)
vFrameshiftIndel <- grepl("frameshift_variant", all_ch$Variant.Consequence)
all_ch$manualreview[vFrameshiftIndel & (all_ch$Gene == "NPM1")] <- TRUE

cat("\n=== MANUAL REVIEW FLAGS ===\n")
cat(sprintf("Total flagged for manual review: %d\n", sum(all_ch$manualreview, na.rm = TRUE)))
# 18

# ========================
# Final summary and output
# ========================
cat("\n=== FINAL WHITELIST + EXCEPTION SUMMARY ===\n")
cat(sprintf("Total variants: %d\n", nrow(all_ch)))
cat(sprintf("Variants in whitelist (including exceptions): %d (%.1f%%)\n",
            sum(all_ch$whitelist),
            100 * sum(all_ch$whitelist) / nrow(all_ch)))
cat(sprintf("  - Missense whitelist: %d\n", sum(all_ch$wl.mis)))
cat(sprintf("  - LoF whitelist: %d\n", sum(all_ch$wl.lof)))
cat(sprintf("  - Splice whitelist: %d\n", sum(all_ch$wl.splice)))
cat(sprintf("  - Exception rules: %d\n", sum(all_ch$wl.exception)))
cat(sprintf("  - Flagged for manual review: %d\n", sum(all_ch$manualreview)))

# Variants in whitelist (including exceptions): 2558 (26.5%)
# > cat(sprintf("  - Missense whitelist: %d\n", sum(all_ch$wl.mis)))
# - Missense whitelist: 79
# > cat(sprintf("  - LoF whitelist: %d\n", sum(all_ch$wl.lof)))
# - LoF whitelist: 250
# > cat(sprintf("  - Splice whitelist: %d\n", sum(all_ch$wl.splice)))
# - Splice whitelist: 2236
# > cat(sprintf("  - Exception rules: %d\n", sum(all_ch$wl.exception)))
# - Exception rules: 36
# > cat(sprintf("  - Flagged for manual review: %d\n", sum(all_ch$manualreview)))
# - Flagged for manual review: 18

# ========================
# Write output files
# ========================
check_vars <- data.frame(
    total_num_variants = nrow(all_ch),
    total_num_whitelist = sum(all_ch$whitelist),
    total_num_manualreview = sum(all_ch$manualreview)
)
print(check_vars)

write.csv(check_vars,
          file.path("ch", "data", "ch_varcount.csv"),
          row.names = FALSE)

write.csv(all_ch,
          file.path("ch", "data", "ch_allvars.csv"),
          row.names = FALSE)

write.csv(all_ch[all_ch$whitelist, ],
          file.path("ch", "data", "ch_wl.csv"),
          row.names = FALSE)

write.csv(all_ch[all_ch$manualreview, ],
          file.path("ch", "data", "ch_manual.csv"),
          row.names = FALSE)


# ========================
# WHITELIST BREAKDOWN FIGURE
# ========================
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

p1 <- ggplot(wl_breakdown, aes(x = Category, y = N, fill = Category)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = paste0(N, "\n(", round(100 * N / sum(N), 1), "%)")),
              vjust = -0.4, size = 3.2, color = "gray30") +
    # scale_fill_manual(values = colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25)),
                       labels = comma) +
    labs(
        title    = "CHIP variant whitelist breakdown",
        # subtitle = paste0("Total whitelisted: ", sum(all_ch$whitelist),
        #                   " | Manual review: ", sum(all_ch$manualreview),
        #                   " | QC-pass input: ", nrow(all_ch)),
        x        = NULL,
        y        = "Number of variants"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        legend.position    = "none",
        plot.title         = element_text(face = "bold", size = 13),
        plot.subtitle      = element_text(size = 10, color = "gray50"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text.x        = element_text(size = 11)
    )

ggsave(file.path("ch", "figures", "ch_wl_categories.pdf"), p1, width = 6, height = 5)

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
    geom_text(aes(label = N), vjust = -0.4, size = 3.2, color = "gray30") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(
        title = "Gene-specific exceptions",
        x     = NULL,
        y     = "Number of variants"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        plot.title         = element_text(face = "bold", size = 12),
        axis.text.x        = element_text(size = 11, angle = 45, hjust = 1)
    )

ggsave(file.path("ch", "figures", "ch_wl_gene_exceptions.pdf"), p2, width = 4, height = 4)

