
# ========================
# ch_whitelist_filter.R
# General CH variant whitelist filtering pipeline
# ========================

library(tidyverse)
library(data.table)

# Load whitelists (adjust paths as needed)
load_whitelists <- function() {
    list(
        gList = fread(file.path("ch", "data", "whitelist_filter_20230531", "NEJM_2017_genes_01262020.txt")),
        whitelist.mis = fread(file.path("ch", "data", "whitelist_filter_20230531", "CHIP_missense_vars_cv_04102022.txt")),
        whitelist.splice = fread(file.path("ch", "data", "whitelist_filter_20230531", "CHIP_splice_vars_agb_01262020.txt")),
        whitelist.LoF = fread(file.path("ch", "data", "whitelist_filter_20230531", "CHIP_nonsense_FS_vars_agb_01262020.txt"))
    )
}

# Helper functions
extractProteinChange <- function(hgvsp) {
    if (is.na(hgvsp) || hgvsp == "" || hgvsp == ".") return(NA)
    prot_change <- gsub("^p\\.", "", hgvsp)
    prot_change <- gsub("%3D", "=", prot_change)
    return(prot_change)
}

extract_exon_number <- function(exon_field) {
    if (is.na(exon_field) || exon_field == "" || exon_field == ".") return(NA)
    exon_num <- as.numeric(gsub("\\|.*", "", exon_field))
    return(exon_num)
}

extract_aa_position <- function(prot_change) {
    if (is.na(prot_change) || prot_change == "" || prot_change == ".") return(NA)
    position <- as.numeric(gsub("\\D", "", prot_change))
    return(position)
}

convert_1to3 <- function(aa_change) {
    if (is.na(aa_change) || aa_change == "" || aa_change == ".") return(NA)

    aa_map <- c(
        "A" = "Ala", "C" = "Cys", "D" = "Asp", "E" = "Glu", "F" = "Phe",
        "G" = "Gly", "H" = "His", "I" = "Ile", "K" = "Lys", "L" = "Leu",
        "M" = "Met", "N" = "Asn", "P" = "Pro", "Q" = "Gln", "R" = "Arg",
        "S" = "Ser", "T" = "Thr", "V" = "Val", "W" = "Trp", "Y" = "Tyr",
        "*" = "Ter", "=" = "="
    )

    pattern <- "^([A-Z*])(\\d+)([A-Z*=])$"

    if (grepl(pattern, aa_change)) {
        parts <- str_match(aa_change, pattern)
        ref_aa <- parts[2]
        position <- parts[3]
        alt_aa <- parts[4]

        ref_3 <- ifelse(ref_aa %in% names(aa_map), aa_map[ref_aa], ref_aa)
        alt_3 <- ifelse(alt_aa %in% names(aa_map), aa_map[alt_aa], alt_aa)

        return(paste0(ref_3, position, alt_3))
    } else {
        return(aa_change)
    }
}

# ========================
# MAIN FILTERING FUNCTION
# ========================
filter_ch_variants <- function(data, dataset_label, AF_LOWER = 0.02, AF_UPPER = 0.45,
                               DP = 15, AD = 2, output_dir = "ch/output") {

    # Generate datetime stamp and log file
    datetime_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    log_file <- file.path("ch", "log", paste0(datetime_stamp, "_", dataset_label, "_filtering.txt"))
    dir.create(file.path("ch", "log"), showWarnings = FALSE, recursive = TRUE)
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    # Start logging
    log_conn <- file(log_file, open = "wt")
    sink(log_conn, type = "output")
    sink(log_conn, type = "message")

    cat("============================================\n")
    cat("CH VARIANT WHITELIST FILTERING PIPELINE\n")
    cat("============================================\n")
    cat("Dataset:", dataset_label, "\n")
    cat("Timestamp:", datetime_stamp, "\n")
    cat("Parameters:\n")
    cat("  AF range:", AF_LOWER, "-", AF_UPPER, "\n")
    cat("  Min depth:", DP, "\n")
    cat("  Min alt depth:", AD, "\n\n")

    # Load whitelists
    cat("Loading whitelists...\n")
    wl <- load_whitelists()

    # ========================
    # QC FILTERING
    # ========================
    cat("\n=== QUALITY CONTROL FILTERING ===\n")
    cat("Initial variants:", nrow(data), "\n\n")

    # Initialize flags
    data <- data %>%
        mutate(
            pass_protein_coding = FALSE,
            pass_vaf = FALSE,
            pass_depth = FALSE,
            qc_pass = FALSE
        )

    # Step 1: Protein coding filter
    data$pass_protein_coding <- (data$Bio.type == "protein_coding") &
        (data$Sample.AltFrac != ".") &
        (data$Sample.Depth != ".") &
        (data$Sample.AltDepth != ".")

    # Convert to numeric
    data <- data %>%
        mutate(
            Sample.AltFrac = as.numeric(Sample.AltFrac),
            Sample.Depth = as.numeric(Sample.Depth),
            Sample.AltDepth = as.numeric(Sample.AltDepth)
        )

    # Remove -1 values
    data$pass_protein_coding <- data$pass_protein_coding &
        !is.na(data$Sample.AltFrac) &
        !is.na(data$Sample.Depth) &
        !is.na(data$Sample.AltDepth) &
        (data$Sample.AltFrac != -1) &
        (data$Sample.Depth != -1) &
        (data$Sample.AltDepth != -1)

    cat("Step 1 - Protein coding filter:", sum(data$pass_protein_coding),
        sprintf("(%.1f%%)\n", 100 * sum(data$pass_protein_coding) / nrow(data)))

    # Step 2: VAF filter
    data$pass_vaf <- data$pass_protein_coding &
        (data$Sample.AltFrac >= AF_LOWER) &
        (data$Sample.AltFrac <= AF_UPPER)

    cat("Step 2 - VAF filter:", sum(data$pass_vaf),
        sprintf("(%.1f%%)\n", 100 * sum(data$pass_vaf) / nrow(data)))

    # Step 3: Depth filter
    data$pass_depth <- data$pass_vaf &
        (data$Sample.Depth >= DP) &
        (data$Sample.AltDepth >= AD)

    cat("Step 3 - Depth filter:", sum(data$pass_depth),
        sprintf("(%.1f%%)\n", 100 * sum(data$pass_depth) / nrow(data)))

    data$qc_pass <- data$pass_depth

    cat("\nFinal QC pass:", sum(data$qc_pass),
        sprintf("(%.1f%% of initial)\n\n", 100 * sum(data$qc_pass) / nrow(data)))

    # Work with passing variants
    all_ch <- data %>% filter(qc_pass)

    # ========================
    # WHITELIST FILTERING
    # ========================
    cat("=== WHITELIST FILTERING ===\n")

    # Extract protein change
    all_ch$ProteinChange <- sapply(all_ch$HGVSp, extractProteinChange)
    all_ch$ExonNumber <- sapply(all_ch$EXON, extract_exon_number)
    all_ch$AAPosition <- sapply(all_ch$ProteinChange, extract_aa_position)

    # Initialize whitelist flags
    all_ch <- all_ch %>%
        mutate(
            whitelist = FALSE,
            wl.mis = FALSE,
            wl.lof = FALSE,
            wl.splice = FALSE,
            wl.exception = FALSE,
            manualreview = FALSE
        )

    # 1) MISSENSE variants
    vmis <- grepl("missense", all_ch$Variant.Consequence)
    wl$whitelist.mis$AAChange_3letter <- sapply(wl$whitelist.mis$AAChange, convert_1to3)

    vmis_wl <- paste(all_ch$Gene, all_ch$ProteinChange, sep = "_") %in%
        paste(wl$whitelist.mis$Gene, wl$whitelist.mis$AAChange_3letter, sep = "_")

    all_ch$whitelist[vmis & vmis_wl] <- TRUE
    all_ch$wl.mis[vmis & vmis_wl] <- TRUE

    cat("Missense whitelist variants:", sum(vmis & vmis_wl), "\n")

    # 2) LOF/FRAMESHIFT variants
    vlof <- grepl("Ter|\\*|X", all_ch$ProteinChange) |
        grepl("stop_gain|stop_lost|start_lost", all_ch$Variant.Consequence) |
        grepl("fs", all_ch$ProteinChange, fixed = TRUE) |
        grepl("frameshift", all_ch$Variant.Consequence)

    vLOFgene <- all_ch$Gene %in% wl$whitelist.LoF$Gene

    all_ch$whitelist[vlof & vLOFgene] <- TRUE
    all_ch$wl.lof[vlof & vLOFgene] <- TRUE

    cat("LoF/frameshift whitelist variants:", sum(vlof & vLOFgene, na.rm = TRUE), "\n")

    # 3) SPLICE variants
    vSplice <- grepl("splice", all_ch$Variant.Consequence)
    vSplicegene <- all_ch$Gene %in% wl$whitelist.splice$Gene

    all_ch$whitelist[vSplice & vSplicegene] <- TRUE
    all_ch$wl.splice[vSplice & vSplicegene] <- TRUE

    cat("Splice whitelist variants:", sum(vSplice & vSplicegene), "\n\n")

    # ========================
    # GENE-SPECIFIC EXCEPTIONS
    # ========================
    cat("=== GENE-SPECIFIC EXCEPTION RULES ===\n")

    vFS <- grepl("fs", all_ch$ProteinChange, fixed = TRUE) |
        grepl("frameshift_variant", all_ch$Variant.Consequence)
    vSpliceCorrectTranscript <- TRUE  # Temporary

    # ASXL1
    vexon11 <- all_ch$ExonNumber == 11
    vexon12 <- all_ch$ExonNumber == 12

    asxl1Exception <- (all_ch$Gene == "ASXL1") & (vlof | vFS) & (vexon11 | vexon12)
    asxl1ExceptionSplice <- (all_ch$Gene == "ASXL1") & vSplice & vSpliceCorrectTranscript & (vexon11 | vexon12)

    all_ch$whitelist[asxl1Exception | asxl1ExceptionSplice] <- TRUE
    all_ch$wl.lof[asxl1Exception] <- TRUE
    all_ch$wl.splice[asxl1ExceptionSplice] <- TRUE
    all_ch$wl.exception[asxl1Exception | asxl1ExceptionSplice] <- TRUE

    cat("ASXL1 exceptions:", sum(asxl1Exception | asxl1ExceptionSplice, na.rm = TRUE), "\n")

    # ASXL2
    asxl2Exception <- (all_ch$Gene == "ASXL2") & (vlof | vFS) & (vexon11 | vexon12)
    asxl2ExceptionSplice <- (all_ch$Gene == "ASXL2") & vSplice & vSpliceCorrectTranscript & (vexon11 | vexon12)

    all_ch$whitelist[asxl2Exception | asxl2ExceptionSplice] <- TRUE
    all_ch$wl.lof[asxl2Exception] <- TRUE
    all_ch$wl.splice[asxl2ExceptionSplice] <- TRUE
    all_ch$wl.exception[asxl2Exception | asxl2ExceptionSplice] <- TRUE

    cat("ASXL2 exceptions:", sum(asxl2Exception | asxl2ExceptionSplice, na.rm = TRUE), "\n")

    # PPM1D
    vexon5 <- all_ch$ExonNumber == 5
    vexon6 <- all_ch$ExonNumber == 6

    ppm1dException <- (all_ch$Gene == "PPM1D") & (vlof | vFS) & (vexon5 | vexon6)

    all_ch$whitelist[ppm1dException] <- TRUE
    all_ch$wl.lof[ppm1dException] <- TRUE
    all_ch$wl.exception[ppm1dException] <- TRUE

    cat("PPM1D exceptions:", sum(ppm1dException, na.rm = TRUE), "\n")

    # TET2
    TETidx <- which(all_ch$Gene == "TET2" & vmis & !is.na(all_ch$AAPosition))

    for(i in TETidx) {
        AApos <- all_ch$AAPosition[i]
        if(!is.na(AApos) && ((AApos >= 1104 & AApos <= 1481) | (AApos >= 1843 & AApos <= 2002))) {
            all_ch$whitelist[i] <- TRUE
            all_ch$wl.mis[i] <- TRUE
            all_ch$wl.exception[i] <- TRUE
        }
    }

    cat("TET2 catalytic domain exceptions:", sum(all_ch$wl.exception & all_ch$Gene == "TET2", na.rm = TRUE), "\n")

    # CBL
    CBLidx <- which(all_ch$Gene == "CBL" & vmis & !is.na(all_ch$AAPosition))

    for(i in CBLidx) {
        AApos <- all_ch$AAPosition[i]
        if(!is.na(AApos) && AApos >= 381 & AApos <= 421) {
            all_ch$whitelist[i] <- TRUE
            all_ch$wl.mis[i] <- TRUE
            all_ch$wl.exception[i] <- TRUE
        }
    }

    cat("CBL RING finger exceptions:", sum(all_ch$wl.exception & all_ch$Gene == "CBL", na.rm = TRUE), "\n")

    # CBLB
    CBLBidx <- which(all_ch$Gene == "CBLB" & vmis & !is.na(all_ch$AAPosition))

    for(i in CBLBidx) {
        AApos <- all_ch$AAPosition[i]
        if(!is.na(AApos) && AApos >= 372 & AApos <= 412) {
            all_ch$whitelist[i] <- TRUE
            all_ch$wl.mis[i] <- TRUE
            all_ch$wl.exception[i] <- TRUE
        }
    }

    cat("CBLB RING finger exceptions:", sum(all_ch$wl.exception & all_ch$Gene == "CBLB", na.rm = TRUE), "\n\n")

    # ========================
    # MANUAL REVIEW FLAGS
    # ========================
    cat("=== MANUAL REVIEW FLAGS ===\n")

    vNonFrameshiftDel <- grepl("inframe_deletion", all_ch$Variant.Consequence)
    vNonFrameshiftIns <- grepl("inframe_insertion", all_ch$Variant.Consequence)

    all_ch$manualreview[(vlof | vFS | vSplice) & (all_ch$Gene == "GATA3")] <- TRUE
    all_ch$manualreview[vNonFrameshiftDel & (all_ch$Gene == "CREBBP")] <- TRUE
    all_ch$manualreview[(vlof | vFS | vSplice) & (all_ch$Gene == "CSF3R")] <- TRUE
    all_ch$manualreview[vNonFrameshiftDel & (all_ch$Gene == "DNMT3A")] <- TRUE
    all_ch$manualreview[vNonFrameshiftDel & (all_ch$Gene == "EP300")] <- TRUE
    all_ch$manualreview[(vNonFrameshiftDel | vNonFrameshiftIns) & (all_ch$Gene == "FLT3")] <- TRUE
    all_ch$manualreview[(vNonFrameshiftDel | vNonFrameshiftIns) & (all_ch$Gene == "JAK2")] <- TRUE
    all_ch$manualreview[vNonFrameshiftDel & (all_ch$Gene == "KDM6A")] <- TRUE
    all_ch$manualreview[(vNonFrameshiftDel | vNonFrameshiftIns) & (all_ch$Gene == "KIT")] <- TRUE
    all_ch$manualreview[(vNonFrameshiftDel | vNonFrameshiftIns) & (all_ch$Gene == "MPL")] <- TRUE
    all_ch$manualreview[grepl("frameshift_variant", all_ch$Variant.Consequence) & (all_ch$Gene == "NPM1")] <- TRUE

    cat("Total flagged for manual review:", sum(all_ch$manualreview, na.rm = TRUE), "\n\n")

    # ========================
    # FINAL SUMMARY
    # ========================
    cat("=== FINAL SUMMARY ===\n")
    cat("Total variants passing QC:", nrow(all_ch), "\n")
    cat("Variants in whitelist:", sum(all_ch$whitelist),
        sprintf("(%.1f%%)\n", 100 * sum(all_ch$whitelist) / nrow(all_ch)))
    cat("  - Missense:", sum(all_ch$wl.mis), "\n")
    cat("  - LoF:", sum(all_ch$wl.lof), "\n")
    cat("  - Splice:", sum(all_ch$wl.splice), "\n")
    cat("  - Exception rules:", sum(all_ch$wl.exception), "\n")
    cat("  - Manual review:", sum(all_ch$manualreview), "\n\n")

    # Close logging
    sink(type = "message")
    sink(type = "output")
    close(log_conn)

    # ========================
    # WRITE OUTPUT FILES
    # ========================
    check_vars <- data.frame(
        total_num_variants = nrow(all_ch),
        total_num_whitelist = sum(all_ch$whitelist),
        total_num_manualreview = sum(all_ch$manualreview)
    )

    write.csv(check_vars,
              file.path(output_dir, paste0(datetime_stamp, "_", dataset_label, "_varcount.csv")),
              row.names = FALSE)

    write.csv(all_ch,
              file.path(output_dir, paste0(datetime_stamp, "_", dataset_label, "_allvariants.csv")),
              row.names = FALSE)

    write.csv(all_ch[all_ch$whitelist, ],
              file.path(output_dir, paste0(datetime_stamp, "_", dataset_label, "_wl.csv")),
              row.names = FALSE)

    write.csv(all_ch[all_ch$manualreview, ],
              file.path(output_dir, paste0(datetime_stamp, "_", dataset_label, "_manualreview.csv")),
              row.names = FALSE)

    cat("Output files saved to:", output_dir, "\n")
    cat("Log file saved to:", log_file, "\n")

    return(list(
        filtered_data = all_ch,
        all_data = data,
        log_file = log_file
    ))
}
