# Update Sample Sheet and Process Data for GRIMAGE2
#### Data Prep ####
# Go to project dir or skip and run in current project dir
setwd("~/project-ophoff/BP-DNAm")

# load packages
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, tidyr, stringr, readr, readxl, data.table, lubridate, tibble)

# load external functions
source('BPDNAm_external_functions.R')

# Read in Bipolar 2023 Sample Sheet.csv as data frame and remove Pool_ID col
BPDNAm_SS <- read.csv("Bipolar 2023 Sample Sheet.csv") %>%
  select(-Pool_ID)

# Read in 2000_sample_covariates.csv as data frame, remove Pool_ID col and rename 'Sample_id' col
BPDNAm_ext <- read.csv("From_Roel/2000_sample_covariates.csv") %>%
  select(-RIN) %>% 
  rename(Sample_Name = Sample_id,
         Age_Years = Sample.Age
         ) %>%
  mutate(Age_Years = as.numeric(gsub("[^0-9.]", "", Age_Years)),
         Gender = recode(Gender, "Female" = "F", "Male" = "M")) # Remove non-numeric characters and convert to numeric

# Read in highcov_technical_covariates.txt as data frame
BPDNAm_cov <- read.table("From_Roel/highcov_technical_covariates.txt", sep = "\t", header = TRUE, row.names = 1) %>%
  t() %>%
  as_tibble(rownames = "Sample_Name") %>%
  select(Sample_Name, !starts_with("/u/project/")) %>%
  rename(Gender = sex,
         Diagnosis = diagnosis,
         Age_Years = age) %>%
  mutate(Age_Years = as.numeric(gsub("[^0-9.]", "", Age_Years))) # Remove non-numeric characters and convert to numeric

# Read in Complete BIG Data.xlsx as data frame
bp_master <- read_excel("Complete BIG Data.xlsx")

# Rename 'Sample_id' in bp_master to 'Sample_Name' for matching and remove duplicate entries
bp_master <- bp_master %>% 
  rename(Sample_Name = Sample_id) %>% 
  group_by(Sample_Name) %>%
  filter(`Date of sample collection` == max(`Date of sample collection`)) %>%
  slice(1) %>%
  ungroup()

# Create BPDNAm_SS missing entries df for samples not present in bp_master
missing_samples <- BPDNAm_SS %>%
  anti_join(bp_master, by = "Sample_Name")

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
  select(Sample_Name, all_of(new_cols), Age_Years, Age_Months)

# Columns to check for missing values
columns_to_check <- c("Age_Years", "Gender", "Diagnosis")

# Merging and filling more missing values from BPDNAm_cov and BPDNAm_ext
BPDNAm_SS_updated <- BPDNAm_SS_updated %>%
  left_join(select(BPDNAm_cov, Sample_Name, Gender, Diagnosis, Age_Years), by = "Sample_Name", suffix = c("", ".cov")) %>%
  left_join(select(BPDNAm_ext, Sample_Name, Gender, Diagnosis, Age_Years), by = "Sample_Name", suffix = c("", ".ext")) %>%
  mutate(
    Gender = coalesce(Gender.cov, Gender.ext, Gender),
    Diagnosis = coalesce(Diagnosis.cov, Diagnosis.ext, Diagnosis),
    Age_Years = coalesce(Age_Years.cov, Age_Years.ext, Age_Years),
    Age_Months = if_else(
      is.na(interval(`Date of birth`, `Date of sample collection`) %>% 
              time_length(unit = "months") %>% 
              floor()),
      Age_Years * 12,
      interval(`Date of birth`, `Date of sample collection`) %>% 
        time_length(unit = "months") %>% 
        floor()
    )
  ) %>%
  select(-ends_with(".cov"), -ends_with(".ext"))

# summarize Sample_Names with _Rep pair
BPDNAm_SS_REP <- BPDNAm_SS_updated %>%
  mutate(base_name = str_remove(Sample_Name, "_REP$")) %>%
  group_by(base_name) %>%
  filter(n() == 2 & sum(str_detect(Sample_Name, "_REP$")) == 1) %>%
  ungroup() %>%
  select(-base_name) %>%
  arrange(Sample_Name) %>%
  left_join(BPDNAm_SS %>% select(Sample_Name, Basename, Sample_Plate, Sample_Well, Sentrix_ID, Sentrix_Position, Sample_Group), 
            by = "Sample_Name")

# Export BPDNAm_SS_REP with a filename including the count of "_REP" samples
BPDNAm_SS_REP %>%
  {
    num_rep_samples <- sum(str_detect(.$Sample_Name, "_REP$"))
    filename <- sprintf("BPDNAm_SS_REP_Pairs_%d.csv", num_rep_samples)
    fwrite(., filename)
  }

# summarize NAs in BPDNAm_SS_updated
NA_summary <- BPDNAm_SS_updated %>%
  summarise(across(everything(), ~sum(is.na(.)) / n())) %>%
  pivot_longer(everything(), names_to = "Column", values_to = "NA_Proportion") %>%
  mutate(
    NA_Count = round(NA_Proportion * nrow(bp_master))
  ) %>%
  filter(NA_Count > 0) %>%
  arrange(desc(NA_Proportion))

# Export NA_summary
fwrite(NA_summary, "BPDNAm_NA_Summary.csv")

# Identify columns to keep (NA proportion < 0.2)
cols_to_keep <- NA_summary %>%
  filter(NA_Proportion < 0.2) %>%
  pull(Column) %>%
  setdiff(c("Serum", "Plasma", "Type of sample", "Time of inclusion", "AZU NR"))

# Create a subset of BPDNAm_SS_updated with cols_to_keep and NA rows joining in cols from missing_data
BPDNAm_SS_NAs <- BPDNAm_SS_updated %>%
  select(Sample_Name, all_of(cols_to_keep), -c("Fam_code", "Date of birth", "Date of sample collection", "Age_Months")) %>%
  filter(if_any(everything(), is.na)) %>%
  left_join(missing_samples, by = "Sample_Name") %>%
  select(-c("Sample_Group", "Basename"))

# Export BPDNAm_SS_NAs
BPDNAm_SS_NAs %>%
  {
    num_samples <- nrow(.)
    filename <- sprintf("BPDNAm_NA_Samples_%d.csv", num_samples)
    fwrite(., filename)
  }

# Export comma separated list of Sample_Names in BPDNAm_SS_NAs
BPDNAm_SS_NAs %>%
  pull(Sample_Name) %>%
  {
    sample_names <- .
    filename <- sprintf("BPDNAm_NA_Sample_Names_%d.txt", length(sample_names))
    paste(sample_names, collapse = ",") %>%
      write_file(filename)
  }

# Create a subset of BPDNAm_SS_updated with cols_to_keep and no NA rows
BPDNAm_SS_noNAs <- BPDNAm_SS_updated %>%
  select(Sample_Name, all_of(cols_to_keep), -c("Fam_code", "Date of birth", "Date of sample collection")) %>%
  anti_join(BPDNAm_SS_NAs, by = "Sample_Name")

# Export BPDNAm_SS_noNAs
BPDNAm_SS_noNAs %>%
  {
    num_samples <- nrow(.)
    filename <- sprintf("BPDNAm_noNA_Samples_%d.csv", num_samples)
    fwrite(., filename)
  }

#### Process Data for GRIMAGE2 ####
# load in required packages
pacman::p_load(qs, ggplot2, reshape2, minfi, GenomicRanges, SummarizedExperiment)

# load in mSetSqFlt.qs S4 object of GenomicRatioSet class
qload(mSetSqFlt.qs, nthreads = 36)

# summarize contents of mSetSqFlt S4 object of GenomicRatioSet class
slot_names <- slotNames(mSetSqFlt)
print(slot_names)
print_slots(mSetSqFlt)

# Perform quality control checks
# Extract M-values
M_values <- assays(mSetSqFlt)$M

# Convert M-values to long format for ggplot2
M_values_long <- melt(M_values)

# Plot density
ggplot(M_values_long, aes(x = value, color = Var2)) +
  geom_density() +
  labs(title = "Density Plot of M-values", x = "M-values", y = "Density") +
  theme_minimal()