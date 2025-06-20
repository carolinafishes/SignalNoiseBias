library(Biostrings)
library(tidyverse)

# Define the input FASTA file and taxa of interest
fasta_file <- "~/Documents/SignalNoiseBiasPaper/Archiving/Percomorpha/IQtreeRates/All_Fish.fasta"
taxa <- c(
  "kurtus_gulliveri",
  "apogon_lateralis",
  "odontobutis_obscura",
  "mogurnda_adspersa"
)

# Define suffixes: p, 1, 2, 3
# Assign suffixes 1–4
suffixes <- setNames(as.character(1:4), taxa)


seqs <- readDNAStringSet(fasta_file)
names(seqs) <- str_replace(names(seqs), " .*", "")  # clean headers

# Function to get base frequencies from one DNAString
get_base_freqs <- function(dna_seq) {
  freqs <- alphabetFrequency(dna_seq, baseOnly = TRUE)
  total <- sum(freqs[c("A", "C", "G", "T")])
  freqs[c("A", "C", "G", "T")] / total
}

# Calculate frequencies for the 4 specified taxa
filtered_seqs <- seqs[names(seqs) %in% taxa]

freq_list <- map2(names(filtered_seqs), suffixes[names(filtered_seqs)], function(taxon, suffix) {
  freqs <- get_base_freqs(filtered_seqs[[taxon]])
  tibble(
    label = paste0("piA", suffix), value = freqs["A"]
  ) %>%
    add_row(label = paste0("piC", suffix), value = freqs["C"]) %>%
    add_row(label = paste0("piG", suffix), value = freqs["G"]) %>%
    add_row(label = paste0("piT", suffix), value = freqs["T"])
})

# Combine into one table
taxon_freqs <- bind_rows(freq_list)

# Compute average of the 4 taxa (suffix = 0)
avg_freqs <- taxon_freqs %>%
  mutate(base = str_sub(label, 3, 3)) %>%
  group_by(base) %>%
  summarise(mean_value = mean(value)) %>%
  mutate(label = paste0("pi", base, "0")) %>%
  select(label, value = mean_value)

# Compute base frequencies for the entire alignment (suffix = p)
all_seq_freqs <- get_base_freqs(DNAString(paste0(seqs, collapse = "")))
all_freqs <- tibble(
  label = paste0("piA", "p"), value = all_seq_freqs["A"]
) %>%
  add_row(label = paste0("piC", "p"), value = all_seq_freqs["C"]) %>%
  add_row(label = paste0("piG", "p"), value = all_seq_freqs["G"]) %>%
  add_row(label = paste0("piT", "p"), value = all_seq_freqs["T"])

# Combine all into final output
final_output <- bind_rows(taxon_freqs, avg_freqs, all_freqs)

# Format as required
formatted_lines <- final_output %>%
  mutate(line = sprintf("%s = %.5f;", label, value)) %>%
  pull(line)

# Print final formatted output
cat(paste(formatted_lines, collapse = " "), "\n")