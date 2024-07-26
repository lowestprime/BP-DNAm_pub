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
