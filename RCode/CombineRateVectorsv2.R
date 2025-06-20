library(stringr)

# Set your directory
ratevec_dir <- "~/Documents/SignalNoiseBiasPaper/Archiving/Percomorpha/MathematicaInputs"  # Replace with your folder path

# List all .ratevec files
ratevec_files <- list.files(ratevec_dir, pattern = "\\.ratevec$", full.names = TRUE)

# Function to extract the vector from a file
extract_ratevec <- function(filepath) {
  contents <- readLines(filepath, warn = FALSE)
  line <- paste(contents, collapse = "")  # In case it's split over lines
  # Extract everything inside the braces
  vec_string <- str_extract(line, "\\{.*\\}")
  vec_string <- str_remove_all(vec_string, "[\\{\\}]")  # remove { and }
  vec_nums <- as.numeric(str_split(vec_string, ",")[[1]])
  return(vec_nums)
}

# Combine all values into one giant vector
all_values <- unlist(lapply(ratevec_files, extract_ratevec))

# Format for output
formatted_vector <- paste0("ratevectorp = {", paste(format(all_values, scientific = FALSE), collapse = ","), "};")

# Write to a file
writeLines(formatted_vector, file.path(ratevec_dir, "combined_AllLoci.ratevec"))