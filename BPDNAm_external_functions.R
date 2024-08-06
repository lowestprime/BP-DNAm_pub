# Function to get total memory from /proc/meminfo or free
get_memory_info <- function() {
  meminfo <- readLines("/proc/meminfo")
  
  # Extract specific memory information
  mem_total_kb <- as.numeric(gsub("[^0-9]", "", grep("MemTotal", meminfo, value = TRUE)))
  mem_available_kb <- as.numeric(gsub("[^0-9]", "", grep("MemAvailable", meminfo, value = TRUE)))
  
  # Convert to gigabytes
  mem_total_gb <- mem_total_kb / 1024 / 1024
  mem_available_gb <- mem_available_kb / 1024 / 1024
  
  # Create a list with all memory information
  mem_info <- list(
    total = mem_total_gb,
    available = mem_available_gb
  )
  
  # Print detected memory information
  message("Total memory (GB): ", mem_info$total)
  message("Available memory (GB): ", mem_info$available)
}

# Custom core detection function
detect_custom_cores <- function() {
  num_cores <- tryCatch({
    availableCores(
      methods = c("SGE", "system", "mc.cores", "nproc"),
      which = "min",
      omit = 0L
    )
  }, error = function(e) {
    as.numeric(system("nproc", intern = TRUE))
  })
  
  # Print detected resources
  message("Number of cores available: ", num_cores)
}

# Function to recursively print slots of an S4 object
print_slots <- function(obj) {
  slot_names <- slotNames(obj)
  for (slot_name in slot_names) {
    cat("Slot name:", slot_name, "\n")
    slot_content <- slot(obj, slot_name)
    str(slot_content)
    cat("\n")
  }
}

# Function to process each table (searching by Sample_Name)
process_table_by_name <- function(table, table_name, missing_samples) {
  table %>%
    filter(Sample_Name %in% missing_samples) %>%
    select(any_of(columns_to_include)) %>%
    mutate(Source_Table = table_name, Search_Method = "Sample_Name")
}

# Function to process each table (searching by Sentrix_ID)
process_table_by_id <- function(table, table_name, missing_sentrix_ids) {
  if ("Sentrix_ID" %in% names(table)) {
    table %>%
      filter(Sentrix_ID %in% missing_sentrix_ids) %>%
      select(any_of(columns_to_include)) %>%
      mutate(Source_Table = table_name, Search_Method = "Sentrix_ID")
  } else {
    tibble()
  }
}

# Custom function to get the first non-NA value
first_non_na <- function(x) {
  non_na <- na.omit(x)
  if (length(non_na) > 0) non_na[1] else NA
}

# External function to calculate Y.pred for a single Y
calculate_Y_pred <- function(Y, cpgs, beta_values, sample_annotation) {
  cpgs1 <- cpgs[Y.pred == Y]
  cols <- match(cpgs1$var, colnames(beta_values))
  cols <- cols[!is.na(cols)]
  X <- colMeans(beta_values[, cols])
  X <- c(1, sample_annotation$Age_Years[1], X)  # Add intercept and age
  as.numeric(X %*% cpgs1$beta)
}

# External function to scale predictions
F_scale <- function(INPUT0, Y.pred0.name, Y.pred.name, gold) {
  out.para <- gold[var == 'COX']
  out.para.age <- gold[var == 'Age']
  m.age <- out.para.age$mean
  sd.age <- out.para.age$sd
  Y0 <- INPUT0[[Y.pred0.name]]
  Y <- (Y0 - out.para$mean) / out.para$sd
  (Y * sd.age) + m.age
}

# Function to load and prepare inputs
prepare_inputs <- function() {
  sample_annotation <- fread("~/project-ophoff/BP-DNAm/BPDNAm_SS_updated_2464.csv")
  Density_Data <- qread("/u/scratch/c/cobeaman/Density_Data.qs", nthreads = 36)
  common_samples <- intersect(Density_Data$sample_groups, sample_annotation$Sample_Name)
  sample_annotation <- sample_annotation %>%
    filter(Sample_Name %in% common_samples) %>%
    arrange(match(Sample_Name, common_samples))
  setkey(sample_annotation, Basename)
  
  beta_values_file <- "beta_values.bk"
  beta_values_desc <- "beta_values.desc"
  beta_values <- filebacked.big.matrix(nrow = nrow(Density_Data$beta_values), 
                                       ncol = ncol(Density_Data$beta_values),
                                       type = "double",
                                       backingfile = beta_values_file,
                                       descriptorfile = beta_values_desc)
  beta_values[,] <- Density_Data$beta_values
  list(sample_annotation = sample_annotation, beta_values = beta_values, 
       beta_values_file = beta_values_file, beta_values_desc = beta_values_desc)
}

# Function to load GrimAge2 source code data
load_grimage2_data <- function() {
  grimage2 <- readRDS("input/DNAmGrimAge2_final.Rds")
  cpgs <- as.data.table(grimage2[[1]])
  setkey(cpgs, Y.pred)
  list(cpgs = cpgs, glmnet_final1 = as.data.table(grimage2[[2]]), gold = as.data.table(grimage2[[3]]))
}

# Function to calculate DNAm protein variables using parallel processing
calculate_protein_variables <- function(Ys, cpgs, beta_values, sample_annotation) {
  cores <- detectCores() - 1
  registerDoParallel(cores)
  results <- foreach(Y = Ys, .combine = cbind, .packages = c("bigmemory", "data.table", "foreach", "doParallel")) %dopar% {
    calculate_Y_pred(Y, cpgs, beta_values, sample_annotation)
  }
  stopImplicitCluster()
  setnames(as.data.table(results), Ys)
}

# Function to scale predictions and calculate age acceleration
scale_predictions <- function(results, sample_annotation, gold) {
  output <- data.table(Sample_Name = sample_annotation$Sample_Name, Age = sample_annotation$Age_Years)
  output[, DNAmGrimAge2 := F_scale(as.data.table(results), 'COX', 'DNAmGrimAge2', gold)]
  output[, AgeAccelGrim2 := DNAmGrimAge2 - Age]
}

# Function to save output and clean up
save_output <- function(output, beta_values_file, beta_values_desc) {
  fwrite(output, "output/myDNAmGrimAge2.csv")
  print(summary(output))
  unlink(c(beta_values_file, beta_values_desc))
}

# Function to generate DNAmGrimAge2 and AgeAccelGrim2
generate_GrimAge2_AgeAccelGrim2 <- function(sample_annotation, Ys, glmnet_final1, gold) {
  vars <- c('Sample_Name', 'Age_Years', 'Gender', Ys)
  output_all <- sample_annotation[, ..vars]
  setnames(output_all, "Age_Years", "Age")
  output_all[, `:=`(Female = as.integer(Gender == "F"), Gender = NULL)]
  
  model_vars <- intersect(names(output_all), glmnet_final1$var)
  output_all[, COX := as.numeric(as.matrix(.SD) %*% glmnet_final1[var %in% model_vars, beta]), .SDcols = model_vars]
  
  output_all <- F_scale(output_all, 'COX', 'DNAmGrimAge2', gold)
  output_all[, AgeAccelGrim2 := resid(lm(DNAmGrimAge2 ~ Age, data = .SD)), .SDcols = c("DNAmGrimAge2", "Age")]
  output_all[, COX := NULL]
  
  old_names <- c('DNAmadm', 'DNAmCystatin_C', 'DNAmGDF_15', 'DNAmleptin', 
                 'DNAmpai_1', 'DNAmTIMP_1', 'DNAmlog.CRP', 'DNAmlog.A1C')
  new_names <- c('DNAmADM', 'DNAmCystatinC', 'DNAmGDF15', 'DNAmLeptin', 
                 'DNAmPAI1', 'DNAmTIMP1', 'DNAmlogCRP', 'DNAmlogA1C')
  setnames(output_all, old_names, new_names, skip_absent = TRUE)
  
  fwrite(output_all, "output/myDNAmGrimAge2.csv")
}
