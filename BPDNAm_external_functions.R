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
  print(paste("Starting calculation for Y:", Y))
  
  cpgs1 <- cpgs[Y.pred == Y]
  print(paste("cpgs1 dimensions:", dim(cpgs1)))
  print(paste("First few entries in cpgs1$var:", head(cpgs1$var, 10)))
  
  cols <- match(cpgs1$var, colnames(beta_values))
  print(paste("Initial matched cols length:", length(cols)))
  cols <- cols[!is.na(cols)]
  print(paste("Filtered matched cols length:", length(cols)))
  
  # Debugging information before column means calculation
  print(paste("Y:", Y))
  print(paste("Number of matched columns:", length(cols)))
  if (length(cols) > 0) {
    print(paste("Matched columns (first 10):", paste(head(cols, 10), collapse = ", ")))
  }
  
  # Ensure cols are not empty
  if (length(cols) == 0) {
    stop("No columns matched in beta_values for Y:", Y)
  }
  
  # Use biganalytics for efficient column means
  print("Calculating column means for beta_values...")
  X <- colMeans(beta_values[, cols])
  print(paste("X length after colMeans:", length(X)))
  
  print("Adding intercept and age to X...")
  X <- c(1, sample_annotation$Age_Years[1], X)  # Add intercept and age
  print(paste("X length after adding intercept and age:", length(X)))
  
  # Debugging information after column means calculation
  print(paste("Length of X:", length(X)))
  print(paste("Length of cpgs1$beta:", length(cpgs1$beta)))
  
  if (length(X) != length(cpgs1$beta)) {
    stop("Length mismatch: Length of X and cpgs1$beta do not match")
  }
  
  print("Performing matrix multiplication...")
  Y.pred <- as.numeric(X %*% cpgs1$beta)
  print(paste("Finished calculation for Y:", Y))
  return(Y.pred)
}

# Function to calculate DNAm protein variables using parallel processing
calculate_protein_variables <- function(Ys, cpgs, beta_values, sample_annotation) {
  cores <- detectCores() - 1
  registerDoParallel(cores)
  
  results <- foreach(Y = Ys, .combine = cbind, .packages = c("bigmemory", "data.table", "foreach", "doParallel")) %dopar% {
    tryCatch({
      print(paste("Processing Y:", Y))
      Y.pred <- calculate_Y_pred(Y, cpgs, beta_values, sample_annotation)
      print(paste("Processed Y:", Y))
      return(Y.pred)
    }, error = function(e) {
      msg <- paste("Error in calculating Y.pred for Y:", Y, " - ", e$message)
      print(msg)
      return(NULL)  # Return NULL in case of error to avoid stopping the entire foreach loop
    })
  }
  
  stopImplicitCluster()
  
  # Check if results contain NULL values due to errors
  if (any(sapply(results, is.null))) {
    stop("Errors occurred in foreach loop. Check logs for details.")
  }
  
  results <- as.data.table(results)
  if (ncol(results) == 0) {
    stop("No results were generated. All tasks may have failed.")
  }
  setnames(results, Ys)
  results
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
  # Set bigmemory option to allow changing dimnames
  options(bigmemory.allow.dimnames=TRUE)
  
  # Load sample annotations and density data
  sample_annotation <- fread("~/project-ophoff/BP-DNAm/BPDNAm_SS_updated_2464.csv")
  Density_Data <- qread("/u/scratch/c/cobeaman/Density_Data.qs", nthreads = 36)
  
  # Align data
  common_samples <- intersect(Density_Data$sample_groups, sample_annotation$Sample_Name)
  sample_annotation <- sample_annotation %>%
    filter(Sample_Name %in% common_samples) %>%
    arrange(match(Sample_Name, common_samples))
  setkey(sample_annotation, Basename)
  
  # Convert beta values to matrix and assign column names
  beta_values_matrix <- as.matrix(Density_Data$beta_values)
  colnames(beta_values_matrix) <- sample_annotation$Sample_Name
  
  # Verify column names assignment
  if (is.null(colnames(beta_values_matrix))) {
    stop("Column names not assigned to beta_values_matrix")
  }
  print("Column names assigned to beta_values_matrix:")
  print(head(colnames(beta_values_matrix), 10))
  
  # Create a big.matrix backed by file
  beta_values_file <- "beta_values.bk"
  beta_values_desc <- "beta_values.desc"
  backingpath <- "/u/scratch/c/cobeaman/"
  beta_values <- filebacked.big.matrix(
    nrow = nrow(beta_values_matrix), 
    ncol = ncol(beta_values_matrix),
    type = "double",
    backingfile = beta_values_file,
    backingpath = backingpath,
    descriptorfile = beta_values_desc
  )
  
  # Fill the big.matrix with data from the original matrix
  beta_values[,] <- beta_values_matrix
  
  # Ensure the big.matrix retains column names
  colnames(beta_values) <- colnames(beta_values_matrix)
  if (is.null(colnames(beta_values))) {
    stop("Column names not assigned to big.matrix beta_values")
  }
  print("Column names assigned to big.matrix beta_values:")
  print(head(colnames(beta_values), 10))
  
  list(
    sample_annotation = sample_annotation,
    beta_values = beta_values, 
    beta_values_file = file.path(backingpath, beta_values_file), 
    beta_values_desc = file.path(backingpath, beta_values_desc)
  )
}

# Add cleanup function to delete the backing files
cleanup_temp_files <- function(beta_values_file, beta_values_desc) {
  file.remove(beta_values_file, beta_values_desc)
}

# Function to load GrimAge2 source code data
load_grimage2_data <- function() {
  grimage2 <- readRDS("input/DNAmGrimAge2_final.Rds")
  cpgs <- as.data.table(grimage2[[1]])
  setkey(cpgs, Y.pred)
  list(cpgs = cpgs, glmnet_final1 = as.data.table(grimage2[[2]]), gold = as.data.table(grimage2[[3]]))
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
