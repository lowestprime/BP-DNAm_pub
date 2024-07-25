#### Useful snippets ####

# We can list the libraries that are actually loaded doing
(.packages())

# list all dataframes in Global Environment
names(which(unlist(eapply(.GlobalEnv,is.data.frame))))

# Unload all currently loaded packages using pacman
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

# Check which packages are currently loaded
pacman::p_loaded()

# free unused R memory
gc()

# save large files fast
qsave(beta_values_long,file="beta_values_long.qs",nthreads=36,preset="uncompressed")

# creating ".RData" in current working directory
save.image()

# count and list the number of each unique in dataframe column
summary(factor(merged_data_final$mri_info_deviceserialnumber))
summary(factor(merged_data_final$batch))
summary(factor(c(merged_data_final$mri_info_deviceserialnumber, merged_data_final$batch)))
combined_summary_df <- as.data.frame(as.table(summary(factor(c(merged_data_final$mri_info_deviceserialnumber, merged_data_final$batch)))))

# take me back home
setwd('~/')

# take me back to proj dir
setwd('/u/project/ophoff/cobeaman/BP-DNAm')

# take me to lib dir
setwd('~/R/APPTAINER/h2-rstudio_4.4.0')

# list library paths
.libPaths()

# List Package Functions
ls("package:qs")

# Set the R_LIBS_USER environment variable to specific directory
.libPaths("/u/home/c/cobeaman/R/APPTAINER/h2-rstudio_4.4.0")

# force install of packages to specified lib path with BiocManager::install()
BiocManager::install(c(), lib = "/u/home/c/cobeaman/R/APPTAINER/h2-rstudio_4.4.0")

# force install of packages to specified lib path with install.packages()
install.packages(c(), lib = "/u/home/c/cobeaman/R/APPTAINER/h2-rstudio_4.4.0")

# list column names
paste(names(merged_data_final), collapse = ", ")

# List all subdirectories
list.dirs(path = '~/project-lhernand/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data', full.names = T, recursive = T)

# sum different values in df
sum(merged_data_no_na$ethnicity %in% c("AFR", "AMR", "EUR"))

# get col names within a dataframe matching a specific naming scheme/pattern using grepl
colnames(smri.R5.1.all)[grepl("^smri_vol_scs", colnames(smri.R5.1.all))]

# get summary of specific cols basic data characteristics within a dataframe that match a naming scheme/pattern using grepl
summary(merged_data_normalized)[grepl("^smri_vol_scs", colnames(merged_data_normalized))]

# launch gptstudio
gptstudio:::gptstudio_chat()

# Remove files from R environement
remove()
rm()
# Remove multiple files by pattern
rm(list = ls(pattern = "_df"))

# Restart RStudio session
.rs.restartR()
