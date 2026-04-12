# ========================
# U2AF1 ANALYSIS
# ========================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
library(data.table, quietly = T)

library(here)
library(readr)
library(purrr)

# ========================
# READ IN
# ========================
pileup_dir <- here("ch", "data", "freeze2_u2af1_vars_pileup")

COL_NAMES <- c("Chr", "Start", "REF", "ALT", "AA_Change", "Total_Depth", "Alt_Count")

parse_one <- function(filepath) {
    person_id <- basename(filepath) |>
        sub("_pileup\\.tsv$",  "", x = _) |>
        sub("\\.pileup\\.tsv$", "", x = _)

    read_tsv(
        filepath,
        col_names = COL_NAMES,
        col_types = cols(
            Chr     = col_character(),
            Start       = col_integer(),
            REF       = col_character(),
            ALT       = col_character(),
            AA_Change = col_character(),
            Total_Depth     = col_integer(),
            Alt_Count = col_integer()
        ),
        comment   = ""
    ) |>
        mutate(person_id = person_id, .before = everything())
}


# Collect all matching files
files <- list.files(
    path       = pileup_dir,
    pattern    = "\\.(pileup\\.tsv|pileup\\.tsv)$|_pileup\\.tsv$",
    full.names = TRUE
)

combined <- map(files, parse_one, .progress = TRUE) |>
    list_rbind()

write_tsv(combined, file.path("ch", "data", "ch_f2_u2af1_pileup.tsv"))

# ========================
# ANALYSIS
# ========================
length(unique(combined$person_id))

# look at which of them were split ??
split <- combined %>% filter(Alt_Count != 0)
length(unique(split$person_id))

# 114, but most of these are 1-2....
# ifl these are fine...


