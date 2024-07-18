# Update Sample Sheet and Process Data for GRIMAGE2
#### Data Prep ####
# Go to project Dir
setwd("~/project-ophoff/BP-DNAm")

# load packages
library(pacman)
p_load(dplyr,tidyr,readxl,data.table,lubridate)

# Read in Bipolar 2023 Sample Sheet.csv as data frame and remove Pool_ID col
BPDNAm_SS <- read.csv("~/project-ophoff/BP-DNAm/Bipolar 2023 Sample Sheet.csv")
BPDNAm_SS <- BPDNAm_SS %>% select(-Pool_ID)

# Read in Complete BIG Data.xlsx as data frame
bp_master <- read_excel("~/project-ophoff/BP-DNAm/Complete BIG Data.xlsx")

## Update Sample Sheet ##
# Rename 'Sample_id' in bp_master to 'Sample_Name' for matching and remove duplicate entries
bp_master <- bp_master %>% 
  rename(Sample_Name = Sample_id) %>% 
  group_by(Sample_Name) %>%
  filter(`Date of sample collection` == max(`Date of sample collection`)) %>%
  slice(1) %>%
  ungroup()

# Identify columns in bp_master that are not in BPDNAm_SS
new_cols <- setdiff(colnames(bp_master), colnames(BPDNAm_SS))

# Merge BPDNAm_SS with only the new columns from bp_master
BPDNAm_SS_updated <- BPDNAm_SS %>%
  left_join(bp_master %>% select(Sample_Name, all_of(new_cols)), by = "Sample_Name") %>%
  mutate(
    Age_Months = interval(`Date of birth`, `Date of sample collection`) %>% 
      time_length(unit = "months") %>% 
      floor(),
    Age_Years = interval(`Date of birth`, `Date of sample collection`) %>% 
      time_length(unit = "years") %>% 
      floor()
  ) %>%
  select(Sample_Name, all_of(new_cols), Age_Years, Age_Months, everything())

# summarize NAs in BPDNAm_SS_updated
na_summary <- BPDNAm_SS_updated %>%
  summarise(across(everything(), ~sum(is.na(.)) / n())) %>%
  pivot_longer(everything(), names_to = "Column", values_to = "NA_Proportion") %>%
  mutate(
    NA_Count = round(NA_Proportion * nrow(bp_master))
  ) %>%
  filter(NA_Count > 0) %>%
  arrange(desc(NA_Proportion))

# Export na_summary
fwrite(na_summary, "BPDNAm_na_summary.csv")

# Identify columns to keep (NA proportion < 0.2)
cols_to_keep <- na_summary %>%
  filter(NA_Proportion < 0.2) %>%
  pull(Column) %>%
  setdiff(c("Serum", "Plasma", "Type of sample"))

# Create a subset of BPDNAm_SS_updated with cols_to_keep and NA rows
BPDNAm_SS_NAs <- BPDNAm_SS_updated %>%
  select(Sample_Name, all_of(cols_to_keep)) %>%
  filter(if_any(everything(), is.na))

# Export BPDNAm_SS_NAs
fwrite(BPDNAm_SS_NAs, "BPDNAm_NA_Samples.csv")

# Create a subset of BPDNAm_SS_updated with cols_to_keep and no NA rows
BPDNAm_SS_noNAs <- BPDNAm_SS_updated %>%
  select(Sample_Name, all_of(cols_to_keep)) %>%
  filter(if_all(everything(), ~!is.na(.))) %>%
  select(Sample_Name:`Date of sample collection`, Age_Years, Age_Months, everything())

# Export BPDNAm_SS_noNAs
fwrite(BPDNAm_SS_noNAs, "BPDNAm_Samples_noNAs.csv")

#### Process Data for GRIMAGE2 ####

