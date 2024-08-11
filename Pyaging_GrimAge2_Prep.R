# Clean environment
rm(list=ls())
options(stringsAsFactors = F)

# Load packages
pacman::p_load(dplyr, tidyverse, data.table, reshape2, purrr, qs, ENmix)

# optional log memory use during workflow
# Rprofmem("memory_profile.log")

# Set working directory and source external functions
setwd("~/project-ophoff/Tools/DNAmGrimAgeGitHub")
source('~/project-ophoff/BP-DNAm/BPDNAm_external_functions.R')

# Load the data
beta_values <- qread("/u/scratch/c/cobeaman/Density_Data.qs", nthreads = 36)$beta_values
# grimage2_cpgs <- (readRDS("input/DNAmGrimAge2_final.Rds")[[1]])$var
# grimage2_cpgs <- grimage2_cpgs[-c(1, 2)]
sample_annotation <- fread("~/project-ophoff/BP-DNAm/BPDNAm_SS_updated_2464.csv")

# Clean CpG names
beta_values_c <- rm.cgsuffix(beta_values)
rownames(beta_values_c) <- gsub("_.*$", "", rownames(beta_values_c))
# Subset beta_values to include only GrimAge2 CpGs
# beta_values_subset <- beta_values[rownames(beta_values) %in% unique(grimage2_cpgs), ]

# Handle duplicate CpG names
duplicate_cpgs <- duplicated(rownames(beta_values_c))
if(any(duplicate_cpgs)) {
  cat("Warning: Found", sum(beta_values_c), "duplicate CpG names. Using the first occurrence of each.\n")
  beta_values_c <- beta_values_c[!duplicate_cpgs, ]
}

# Convert to data.table
beta_values_dt <- as.data.table(t(beta_values_c), keep.rownames = "Basename")
# head(colnames(beta_values_dt))

# Convert beta_values_dt to a data frame and perform a left join with selected columns from sample_annotation
temp_data <- beta_values_dt %>%
  as.data.frame() %>%
  left_join(sample_annotation %>%
              select(Basename, Sample_Name, Gender, Diagnosis, Age_Years) %>%
              mutate(Female = case_when(
                Gender == "F" ~ 1,
                Gender == "M" ~ 0,
                TRUE ~ NA_real_
              )),
            by = "Basename")

# Rename the columns using colnames() to avoid issues with non-standard column names
colnames(temp_data)[colnames(temp_data) == "Sample_Name"] <- "SampleID"
colnames(temp_data)[colnames(temp_data) == "Age_Years"] <- "Age"

# Reorder the columns, bringing SampleID, Age, and Female to the front
beta_values_final <- temp_data %>%
  select(SampleID, Age, Female, Diagnosis, everything()) %>%
  as.data.table()

# ID missing grimage2 CpGs in beta_values_final
grimage2_cpgs <- grimage2_cpgs[grepl("^cg", grimage2_cpgs)]
cleaned_grimage2_cpgs <- unique(grimage2_cpgs)
missing_cpgs <- setdiff(cleaned_grimage2_cpgs, names(beta_values_final))
print(missing_cpgs)

# Assuming grimage2_cpgs is already defined
beta_values_subset <- beta_values_final %>%
  select(SampleID, Age, Female, Diagnosis, any_of(cleaned_grimage2_cpgs)) %>%
  as.data.table()

# optional log memory use during workflow
# Rprofmem(NULL)
# summaryRprof("memory_profile.log")

# save to qs object
# qsavem(beta_values_final, beta_values_subset, file="beta_values.qs", preset = "fast", nthreads = 36)

# Save beta_values_final and beta_values_subset as CSVs
# save_with_info(beta_values_final, "mymetharray_final")
# save_with_info(beta_values_subset, "mymetharray_subset")
