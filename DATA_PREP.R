#### Update Sample Sheet and Process Data for GRIMAGE2 ####

#### Data Prep ####
# load packages
library(pacman)
p_load(dplyr,tidyr,readxl)

# Read in Bipolar 2023 Sample Sheet.csv as data frame and remove Pool_ID col
BPDNAm_SS <- read.csv("~/project-ophoff/BP-DNAm/Bipolar 2023 Sample Sheet.csv")
BPDNAm_SS <- BPDNAm_SS %>% select(-Pool_ID)

# Read in Complete BIG Data.xlsx as data frame
bp_master <- read_excel("Complete BIG Data.xlsx")


#### Update Sample Sheet ####
# Rename 'Sample_id' in bp_master to 'Sample_Name' for matching
bp_master <- bp_master %>% rename(Sample_Name = Sample_id)

# Identify columns in bp_master that are not in BPDNAm_SS
new_cols <- setdiff(colnames(bp_master), colnames(BPDNAm_SS))

# Merge BPDNAm_SS with only the new columns from bp_master
BPDNAm_SS_updated <- BPDNAm_SS %>%
  left_join(bp_master %>% select(Sample_Name, all_of(new_cols)), by = "Sample_Name")

# summarize NAs in BPDNAm_SS_updated
na_summary <- BPDNAm_SS_updated %>%
  summarise(across(everything(), ~sum(is.na(.)) / n())) %>%
  pivot_longer(everything(), names_to = "Column", values_to = "NA_Proportion") %>%
  mutate(
    NA_Count = round(NA_Proportion * nrow(bp_master))
  ) %>%
  filter(NA_Count > 0) %>%
  arrange(desc(NA_Proportion))

print(na_summary)

# Identify columns to keep (NA proportion < 0.2)
cols_to_keep <- na_summary %>%
  filter(NA_Proportion < 0.2) %>%
  pull(Column)

# Create a subset of BPDNAm_SS_updated with cols_to_keep and NA rows
BPDNAm_SS_NAs <- BPDNAm_SS_updated %>%
  select(Sample_Name, all_of(cols_to_keep)) %>%
  filter(if_any(everything(), is.na))

# Create a subset of BPDNAm_SS_updated with cols_to_keep and no NA rows
BPDNAm_SS_updated_noNAs <- BPDNAm_SS_updated %>%
  select(Sample_Name, all_of(cols_to_keep)) %>%
  filter(if_all(everything(), ~!is.na(.)))

# Print summaries
cat("Dimensions of BPDNAm_SS_NAs:", dim(BPDNAm_SS_NAs), "\n")
cat("Dimensions of BPDNAm_SS_updated_noNAs:", dim(BPDNAm_SS_updated_noNAs), "\n")

# Verify that the sum of rows in both dataframes equals the total rows in the original
total_rows <- nrow(BPDNAm_SS_updated)
cat("Total rows in original dataframe:", total_rows, "\n")
cat("Sum of rows in NA and no-NA dataframes:", 
    nrow(BPDNAm_SS_NAs) + nrow(BPDNAm_SS_updated_noNAs), "\n")

#### Process Data for GRIMAGE2 ####



