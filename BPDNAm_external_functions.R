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

