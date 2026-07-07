library(readxl)
library(here)

master <- read_excel(here("ch", "data", "brca_carriers_ch_freq_w_sampnum_20250716.xlsx"), sheet = "Data_from_master_table")
tumor <- read_excel(here("ch", "data", "brca_carriers_ch_freq_w_sampnum_20250716.xlsx"), sheet = "Tumor")
treatment <- read_excel(here("ch", "data", "brca_carriers_ch_freq_w_sampnum_20250716.xlsx"), sheet = "Treatment")
recurrence <- read_excel(here("ch", "data", "brca_carriers_ch_freq_w_sampnum_20250716.xlsx"), sheet = "Recurrence")
pathology <- read_excel(here("ch", "data", "brca_carriers_ch_freq_w_sampnum_20250716.xlsx"), sheet = "Pathology")

length(unique(master$globalid))
# 4336 brca carriers

# samples w/ DNA
master <- master[!is.na(master$DNA), ]
length(unique(master$globalid))
# 2816 samples w/ DNA brca carriers

# those w/ cancer
cancer <- master[(master$globalid %in% tumor$globalid),] #2376 w/ cancer total
length(unique(cancer$globalid))
dim(cancer)
# 1560   33

# recurrence
rec <- master[!(master$globalid %in% recurrence$globalid),]
dim(rec) #2441   33

# recurrence
treat <- master[!(master$globalid %in% treatment$globalid),]
dim(treat) #2441   33

# recurrence
treat <- master[!(master$globalid %in% treatment$globalid),]
dim(rec) #2441   33

