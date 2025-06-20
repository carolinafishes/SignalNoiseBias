# Load necessary libraries
library(dplyr)
library(stringr)
library(readr)
library(Biostrings)

# Step 1: Read the rate file
rate_data <- read.table("~/Documents/SignalNoiseBiasPaper/Archiving/Percomorpha/IQtreeRates/sate-gblocks-clean-min-90-taxa.phylip.mlrate", header = TRUE, sep = "\t")

# Step 2: Read and parse the charset file
charset_lines <- readLines("~/Documents/SignalNoiseBiasPaper/Archiving/Percomorpha/IQtreeRates/sate-gblocks-clean-min-90-taxa.charsets")
charset_lines <- charset_lines[grepl("^charset", charset_lines)]

# Extract UCE name and start-end positions
charset_df <- tibble(raw = charset_lines) %>%
  mutate(
    name = str_extract(raw, "'[^']+'") %>% str_replace_all("'", "") %>% str_replace("\\.nexus$", ""),
    range = str_extract(raw, "=\\s*\\d+-\\d+") %>% str_remove("^=\\s*"),
    start = as.integer(str_extract(range, "^\\d+")),
    end = as.integer(str_extract(range, "\\d+$"))
  ) %>%
  select(name, start, end)

# Step 3: Extract and write rate vectors
for (i in seq_len(nrow(charset_df))) {
  locus_name <- charset_df$name[i]
  start_pos <- charset_df$start[i]
  end_pos <- charset_df$end[i]
  
  # Extract rates for this locus
  rates <- rate_data$Rate[start_pos:end_pos]
  
  # Format as ratevectorp = {....};
  rate_line <- paste0("ratevectorp = {", paste0(sprintf("%.5f", rates), collapse = ","), "};")
  
  # Write to file
  writeLines(rate_line, paste0(locus_name, ".ratevec"))
}


####Get UCE loci from concatenated alignment
setwd("~/Documents/SignalNoiseBiasPaper/Archiving/Percomorpha/IndividualLoci_Fasta")
# Step 1: Read charset file
charset_lines <- readLines("~/Documents/SignalNoiseBiasPaper/Archiving/Percomorpha/IQtreeRates/sate-gblocks-clean-min-90-taxa.charsets")
charset_lines <- charset_lines[grepl("^charset", charset_lines)]

charset_df <- tibble(raw = charset_lines) %>%
  mutate(
    name = str_extract(raw, "'[^']+'") %>% str_replace_all("'", "") %>% str_replace("\\.nexus$", ""),
    range = str_extract(raw, "=\\s*\\d+-\\d+") %>% str_remove("^=\\s*"),
    start = as.integer(str_extract(range, "^\\d+")),
    end = as.integer(str_extract(range, "\\d+$"))
  ) %>%
  select(name, start, end)

# Step 2: Load full alignment fasta (aligned sequences must match positions in charset)
full_alignment <- readBStringSet("~/Documents/SignalNoiseBiasPaper/Archiving/Percomorpha/IQtreeRates/All_Fish.fasta", format = "fasta")

# Step 3: Extract and write each UCE locus
for (i in seq_len(nrow(charset_df))) {
  locus_name <- charset_df$name[i]
  start_pos <- charset_df$start[i]
  end_pos <- charset_df$end[i]
  
  # Extract each sequence substring
  locus_seqs <- subseq(full_alignment, start = start_pos, end = end_pos)
  names(locus_seqs) <- names(full_alignment)  # keep original headers
  
  # Write to individual fasta file
  writeXStringSet(locus_seqs, filepath = paste0(locus_name, ".fasta"), format = "fasta")
}

# Function: TRUE if sequence is empty (only ?, N, n, -, or .)
is_empty <- function(seq) {
  cleaned <- gsub("[-Nn?\\.]", "", as.character(seq))
  return(nchar(cleaned) == 0)
}

# Get all FASTA files in current directory
fasta_files <- list.files(pattern = "\\.fasta$")

for (f in fasta_files) {
  cat("Checking", f, "\n")
  
  # Read alignment using BStringSet to preserve all characters
  aln <- readBStringSet(f, format = "fasta")
  
  # Identify non-empty sequences
  non_empty <- aln[!sapply(aln, is_empty)]
  
  # Overwrite with cleaned version if anything remains
  if (length(non_empty) > 0) {
    writeXStringSet(non_empty, filepath = f, format = "fasta")
  } else {
    warning(paste("File", f, "has no remaining sequences after cleaning."))
  }
}