# Update Sample Sheet and Process Data for GRIMAGE2
#### Data Prep ####
# Go to project dir or skip and run in current project dir
setwd("~/project-ophoff/BP-DNAm")

# load packages
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, tidyr, stringr, readr, readxl, data.table, lubridate, tibble)

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
pacman::p_load(qs, ggplot2, plotly, RColorBrewer, reshape2, minfi, GenomicRanges, SummarizedExperiment)

# load external functions
source('BPDNAm_external_functions.R')

# save mSetSqFlt.qs S4 object of GenomicRatioSet class
# load("BipolarMethylationData.RData")
# qsavem(mSetSqFlt, file = "BipolarMethylationData.qs", nthreads = 36, preset = "uncompressed")

# load in mSetSqFlt.qs S4 object of GenomicRatioSet class
qload("BipolarMethylationData.qs", nthreads = 36)

# summarize contents of mSetSqFlt
slot_names <- slotNames(mSetSqFlt)
print(slot_names)
print_slots(mSetSqFlt)

# Visualize the distribution of the beta values for each sample after normalization
# Extract Beta Values from GenomicRatioSet
# beta_values <- getBeta(mSetSqFlt)

# Extract Sample Groups from GenomicRatioSet
# sample_groups <- as.factor(pData(mSetSqFlt)$Sample_Group)

# Plot Density of Beta Values
# densityPlot <- densityPlot(beta_values, sampGroups = sample_groups, main = "Normalized", legend = FALSE)
# densityPlot <- recordPlot()
# load density plot data
qread("Density_Data.qs", nthreads=36)$densityPlot
print(densityPlot)

# Reshape beta_values to long format
# beta_values_long <- melt(beta_values, varnames = c("CpG_Site", "Sample"), value.name = "Beta_Value")
qread(beta_values_long.qs, nthreads = 36)
sample_groups <- qread("Density_Data.qs", nthreads=36)$sample_groups

# Merge with sample_groups
beta_values_long <- beta_values_long %>%
  mutate(SampleGroup = sample_groups[match(beta_values_long$Sample, colnames(beta_values))])

# Plot using ggplot2 with facets
ggplot(beta_long, aes(x = BetaValue, color = SampleGroup)) +
  geom_density() +
  facet_wrap(~ SampleGroup, scales = "free_y") +
  theme_minimal() +
  labs(title = "Density Plot of Beta Values by Sample Group", x = "Beta Value", y = "Density")

#### Remaining Analysis (Work in Progress) ####
# Source: https://shorturl.at/hKHuc

I. R Workflow (Hoffman2)

# --- R SECTION ---

# 1. Load Libraries and Set Working Directory
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19) 
library(dplyr)
library(data.table)
library(dnaMethyAge) 
library(meffil)       
library(methylClock)  
# ... Add other benchmarking libraries as needed (e.g., BioAge) ...

setwd("~/project-ophoff/BP-DNAm")  # Replace with your actual project directory

# 2. Load and Preprocess Methylation Data
qload("BipolarMethylationData.qs", nthreads = 36) # Assuming .qs format
sample_annotation <- fread("BPDNAm_noNA_Samples_2351.csv") 

# Quality Control (adapt based on your data):
# detP <- detectionP(mSetSqFlt)
# failed <- detP > 0.01
# mSetSqFlt <- mSetSqFlt[rowSums(failed) == 0, ] 
# ... other QC steps ...

# Normalization (if not already done):
# mSetSqFlt <- preprocessQuantile(mSetSqFlt)  

# Blood Cell Composition Estimation
cellCounts <- estimateCellCounts(mSetSqFlt, compositeCellType = "Blood", 
                                 referencePlatform = "IlluminaHumanMethylation450k") 

# Extract beta values 
beta_values <- getBeta(mSetSqFlt)

# 3. Sample and CpG Verification 
# ------------------------------

# Get sample names from methylation data
meth_sample_names <- colnames(beta_values) 

# Check for samples in annotation NOT in methylation data
missing_in_meth <- setdiff(sample_annotation$Sample_Name, meth_sample_names)
if (length(missing_in_meth) > 0) {
  warning("Samples in annotation not found in methylation data: ", 
          paste(missing_in_meth, collapse = ", "))
}

# Check for samples in methylation data NOT in annotation
missing_in_annot <- setdiff(meth_sample_names, sample_annotation$Sample_Name)
if (length(missing_in_annot) > 0) {
  warning("Samples in methylation data not found in annotation: ", 
          paste(missing_in_annot, collapse = ", "))
}

# Load GrimAge2 CpG list 
grimage2_cpgs <- fread("DNAmGrimAge2_1030CpGs.csv")

# Get CpG names from methylation data 
data_cpgs <- rownames(mSetSqFlt)

# Check for missing CpGs
missing_cpgs <- setdiff(grimage2_cpgs$var, data_cpgs)
if (length(missing_cpgs) > 0) {
  stop("GrimAge2 CpGs missing in your methylation data: ", 
       paste(missing_cpgs, collapse = ", "))
}

# 4. Calculate GrimAge2 (Using Provided Source Code)
# ---------------------------------------------------
# Load GrimAge2 source code data
grimage2 <- readRDS("DNAmGrimAge2_final.rds") # Replace with your file if needed
cpgs <- grimage2[[1]]
glmnet.final1 <- grimage2[[2]]
gold <- grimage2[[3]]

# Define functions from source code
F_scale <- function(INPUT0, Y.pred0.name, Y.pred.name, gold) {
  out.para <- subset(gold, var == 'COX')
  out.para.age <- subset(gold, var == 'Age')
  m.age <- out.para.age$mean
  sd.age <- out.para.age$sd
  Y0 <- INPUT0[, Y.pred0.name]
  Y <- (Y0 - out.para$mean) / out.para$sd
  INPUT0[, Y.pred.name] <- as.numeric((Y * sd.age) + m.age)
  return(INPUT0)
}

# Step 1: Generate DNAm Protein Variables 
Ys <- unique(cpgs$Y.pred) 
for (k in 1:length(Ys)) {
  cpgs1 <- subset(cpgs, Y.pred == Ys[k])
  Xs <- subset(dat.meth, select = cpgs1$var)
  Y.pred <- as.numeric(as.matrix(Xs) %*% cpgs1$beta)
  dat.meth[, Ys[k]] <- Y.pred
  attr(dat.meth[, Ys[k]], "dimnames") <- NULL
}

# Step 2: Generate DNAmGrimAge2 and AgeAccelGrim2
vars <- c('SampleID', 'Age', 'Female', Ys)
output.all <- subset(dat.meth, select = vars)

# Calculate raw GrimAge2 ('COX')
output.all$COX <- as.numeric(as.matrix(subset(output.all, select = glmnet.final1$var)) %*% glmnet.final1$beta)
output.all <- F_scale(output.all, 'COX', 'DNAmGrimAge2', gold) 

# Calculate AgeAccelGrim2
output.all$DNAmtemp <- output.all$DNAmGrimAge2
output.all$AgeAccelGrim2 <- residuals(lm(DNAmtemp ~ Age, data = output.all, na.action = na.exclude))
output.all$DNAmtemp <- output.all$COX <- NULL

# Rename columns
old.name <- c('DNAmadm', 'DNAmCystatin_C', 'DNAmGDF_15', 'DNAmleptin', 
              'DNAmpai_1', 'DNAmTIMP_1', 'DNAmlog.CRP', 'DNAmlog.A1C')
new.name <- c('DNAmADM', 'DNAmCystatinC', 'DNAmGDF15', 'DNAmLeptin', 
              'DNAmPAI1', 'DNAmTIMP1', 'DNAmlogCRP', 'DNAmlogA1C')
for (k in 1:length(old.name)) {
  id <- which(names(output.all) == old.name[k])
  names(output.all)[id] <- new.name[k]
}

# 5. Run dnaMethyAge Clocks
# --------------------------
dna_methy_age_results <- data.frame(SampleID = sample_annotation$Sample_Name)
dnam_clocks <- c("HannumG2013", "HorvathS2013", "YangZ2016", "ZhangY2017",
                 "HorvathS2018", "LevineM2018", "McEwenL2019", "ZhangQ2019", 
                 "LuA2019", "epiTOC2", "ShirebyG2020", "BernabeuE2023c", 
                 "LuA2023p2", "LuA2023p3") 

for (clock in dnam_clocks) {
  dnam_age <- DNAmAge(beta_values, clock = clock)
  dna_methy_age_results[[clock]] <- dnam_age$Age
  dna_methy_age_results[[paste0("AgeAccel", clock)]] <- dnam_age$AgeAccel
}

# 6. Run DunedinPoAm and DunedinPACE
# -------------------------------------
dunedin_poam <- DunedinPoAm(beta_values)
dunedin_pace <- DunedinPACE(beta_values)

# 7. Run PC-Clocks (if needed)
# ----------------------------
# ... (Code for PC-Clocks similar to dnaMethyAge loop - adapt data format) ... 

# 8. Combine R-based Clock Results and Prepare for Python
# ---------------------------------------------------------
r_clock_results <- data.frame(SampleID = sample_annotation$Sample_Name,
                             DunedinPoAm = dunedin_poam$PoAmAge,
                             DunedinPACE = dunedin_pace$PACE) %>%
                    left_join(dna_methy_age_results, by = "SampleID") 
# ... (Add left_join for PC-Clocks results if calculated) ... 

# Save R-based clock results for Python import
fwrite(r_clock_results, "R_Clock_Results.csv", sep = ",", row.names = F, quote = F)

# --- End of R Section (Part 1) ---
content_copy
Use code with caution.
R
II. Python Workflow (Hoffman2)

# --- PYTHON SECTION ---
import pandas as pd
import pyaging 

# Load methylation data (adapt file path and format)
meth_data = pd.read_csv("mymetharray.csv", index_col="SampleID") # Example for beta values

# Load R-based clock results
r_clock_results = pd.read_csv("R_Clock_Results.csv", index_col="SampleID")

# Calculate ALL pyaging clocks
pyaging_results = pyaging.calculate_all_clocks(meth_data, platform='450k') 

# Select relevant clocks and create a data frame
python_clock_results = pd.DataFrame({
    'SampleID': meth_data.index, 
    'PhenoAge': pyaging_results['phenoage'],
    'AgeAccelPheno': pyaging_results['phenoage_acceleration'],
    'DunedinPoAm_pyaging': pyaging_results['dunedin_poam'], 
    'AgeAccelDunedinPoAm_pyaging': pyaging_results['dunedin_poam_acceleration'],
    'DunedinPACE_pyaging': pyaging_results['dunedin_pace'],
    'AgeAccelDunedinPACE_pyaging': pyaging_results['dunedin_pace_acceleration'],
    'MethylCipherAge': pyaging_results['methyl_cipher'],
    'AgeAccelMethylCipher': pyaging_results['methyl_cipher_acceleration']
    # ... Add other pyaging clock names and age accelerations ...
    })

# Combine with R-based clocks 
python_clock_results = python_clock_results.join(r_clock_results, on="SampleID")

# Save Python-based clock results 
python_clock_results.to_csv("Python_Clock_Results.csv", index=False)

# --- End of Python Section ---
content_copy
Use code with caution.
Python
III. R Workflow (Hoffman2 - Continued)

# --- R SECTION (Part 2) ---

# 9. Import Python Clock Results
python_clock_results <- fread("Python_Clock_Results.csv")

# 10. Combine All Clock Results and Sample Annotation
master_data <- left_join(sample_annotation, output.all, by = "SampleID") %>%
              left_join(python_clock_results, by = "SampleID")

# 11. Comprehensive Benchmarking
clocks_to_benchmark <- names(master_data)[grepl("Age$|GrimAge2$", names(master_data))] # Select all clock variables

benchmark_results <- methylClock::benchmarkClocks(master_data, 
                                                   clocks = clocks_to_benchmark, 
                                                   ageCol = "Age", 
                                                   # ... Other methylClock parameters ...
                                                   )
# ... (Add benchmarking with other tools as desired) ...

print(benchmark_results)
fwrite(benchmark_results, "Clock_Benchmark_Results.csv", sep = ",", row.names = F, quote = F) 

# 12. Create a Composite Clock (using PCR as an example)
# -------------------------------------------------------
# Select top-performing clocks based on benchmarking (using correlation as example)
top_clocks <- benchmark_results %>%
                slice_max(n = 3, order_by = Correlation) %>% # Get top 3 by correlation
                pull(Clock) 
# ... Or select manually based on multiple metrics and relevance to bipolar disorder ...

# Prepare data for PCR 
pcr_data <- master_data[, c(top_clocks, "Age")] 
pcr_data <- na.omit(pcr_data)

# Perform principal component regression
pcr_model <- pcr(Age ~ ., data = pcr_data, validation = "CV") 

# Extract composite clock predictions
master_data$CompositeClock <- predict(pcr_model, newdata = master_data)

# Calculate composite clock age acceleration 
master_data$AgeAccelComposite <- residuals(lm(CompositeClock ~ Age, 
                                                data = master_data, na.action = na.exclude))

# 13. Analyze and Explore Results
# --------------------------------
# - Explore correlations between clocks
# - Analyze associations with bipolar disorder features 
# - ...

# Save master data frame with ALL clock results 
fwrite(master_data, "Bipolar_Epigenetic_Clock_Data.csv", sep = ",", 
       row.names = F, quote = F)

# --- End of R Workflow ---
content_copy
Use code with caution.
R
Execution:

Run the R script (Part 1) on Hoffman2. This will calculate GrimAge2, dnaMethyAge clocks, Dunedin clocks, and prepare data for Python.
Run the Python script on Hoffman2 (request GPU resources). This will calculate all pyaging clocks, including PhenoAge, DunedinPoAm, DunedinPACE, and methylCIPHER, and combine them with the R-based clock results.
Once the Python script is finished, continue running the R script (Part 2) on Hoffman2. This will perform comprehensive benchmarking, create the composite clock, and allow you to explore the integrated results.
This workflow is more efficient, using pyaging's meta-tool capabilities to simplify the Python steps. By thoroughly benchmarking and selectively combining clocks, you will gain powerful insights into biological aging in your bipolar disorder cohort!
