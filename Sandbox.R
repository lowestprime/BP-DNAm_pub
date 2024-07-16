#### Useful snippets ####

# We can list the libraries that are actually loaded doing
(.packages())

# Unload all currently loaded packages using pacman
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

# Check which packages are currently loaded
pacman::p_loaded()

# count and list the number of each unique in dataframe column
summary(factor(merged_data_final$mri_info_deviceserialnumber))
summary(factor(merged_data_final$batch))
summary(factor(c(merged_data_final$mri_info_deviceserialnumber, merged_data_final$batch)))
combined_summary_df <- as.data.frame(as.table(summary(factor(c(merged_data_final$mri_info_deviceserialnumber, merged_data_final$batch)))))

# take me back home
setwd('~/')

# take me back to proj dir
setwd('~/project-lhernand/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis')

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
