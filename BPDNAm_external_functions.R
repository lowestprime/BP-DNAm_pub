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

# Function to calculate Y.pred for a single Y
calculate_Y_pred <- function(Y, cpgs_subset, beta_values) {
  Xs <- t(beta_values[cpgs_subset$var, , drop = FALSE])
  Y.pred <- as.numeric(Xs %*% cpgs_subset$beta)
  return(Y.pred)
}

# Define functions from GrimAge2 source code
F_scale <- function(INPUT0, Y.pred0.name, Y.pred.name, gold) {
  out.para <- gold[var == 'COX']
  out.para.age <- gold[var == 'Age']
  m.age <- out.para.age$mean
  sd.age <- out.para.age$sd
  Y0 <- INPUT0[[Y.pred0.name]]
  Y <- (Y0 - out.para$mean) / out.para$sd
  INPUT0[, (Y.pred.name) := as.numeric((Y * sd.age) + m.age)]
  return(INPUT0)
}