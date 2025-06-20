# Load libraries
library(tidyverse)

# Define folder with CSVs
data_folder <- "~/Documents/SignalNoiseBiasPaper/Archiving/Birds/Results"  # Replace with your actual folder path

# Get list of relevant CSV files (exclude any with "Check" in the name)
file_list <- list.files(data_folder, pattern = "\\.csv$", full.names = TRUE)
file_list <- file_list[!grepl("Check", file_list, ignore.case = TRUE)]

#figure out which files have the noise+bias problem
  # Initialize list to track files of interest
files_bias_only <- c()

  # Initialize list to track ratios 
ratios <- data.frame(nrow=length(file_list), ncol=5)

# Loop through files
for (i in 1:length(file_list)) {
  file<-file_list[i]
  df <- read.csv(file)
  
  if (nrow(df) == 0) next
  
  last_row <- tail(df, 1)
  mu_signal <- as.numeric(last_row$MuValueSignal)
  mu_noise <- as.numeric(last_row$MuValueNoise)
  mu_bias <- as.numeric(last_row$MuValueBias)
  mu_ratio <- as.numeric(last_row$MuValueNoise)/as.numeric(last_row$MuValueSignal)

  
  # Evaluate the condition
  condition1 <- mu_bias > mu_signal
  
  # Track files where condition3 is TRUE but condition1 is FALSE
  if (condition1 ) {
    files_bias_only <- c(files_bias_only, basename(file))
  }
  ratios[i,1]<-file
    ratios[i,2]<-mu_signal
      ratios[i,3]<-mu_noise
        ratios[i,4]<-mu_bias
          ratios[i,5]<-mu_ratio
}

# View result
print("Files where (Noise + Bias > Signal) but (Noise <= Signal):")
print(files_bias_only)

colnames(ratios)<-c("locus","signal","noise","bias","Noise/Signal")

#plot data trends
ggplot(ratios, aes(x = `Noise/Signal`)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = quantile(ratios$`Noise/Signal`, probs = c(0.25, 0.5, 0.75)),
             linetype = "dashed", color = "red") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Density of Noise/Signal Ratios",
    x = "Noise / Signal",
    y = "Density"
  )
  
  ggplot(ratios, aes(x = `Noise/Signal`)) +
  geom_histogram(binwidth = 0.5, fill = "steelblue", color = "white", alpha = 0.8) +
  geom_vline(xintercept = quantile(ratios$`Noise/Signal`, probs = c(0.25, 0.5, 0.75)), 
             linetype = "dashed", color = "red") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Distribution of Noise/Signal Ratios",
    x = "Noise / Signal",
    y = "Count"
  )
  
# plot loci:
ratios_sorted <- ratios %>%
  mutate(locus_id = basename(locus)) %>%
  arrange(`Noise/Signal`) %>%
  mutate(locus_id = factor(locus_id, levels = locus_id))  # preserve sort order

# Barplot
ggplot(ratios_sorted, aes(x = locus_id, y = `Noise/Signal`)) +
  geom_col(fill = "darkgreen") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6)
  ) +
  labs(
    title = "Noise/Signal Ratio by Locus (Sorted)",
    x = "Locus",
    y = "Noise / Signal"
  )
  
 #to get sorted vector for mathematica assumong you already worked with the other file
 # Set your directory
ratevec_dir <- "~/Documents/SignalNoiseBiasPaper/Archiving/Birds/MathematicaInputs"  # Replace with your folder path
ratevec_files<-as.vector(ratios_sorted$locus)
ratevec_files_short <- str_replace(ratevec_files, "^/Users/adornbur", "~")
ratevec_files_short <- str_replace(ratevec_files_short, "Results/", "MathematicaInputs/")
ratevec_files_short <- str_replace(ratevec_files_short, ".csv", ".fasta.ratevec")
# Combine all values into one giant vector
all_values <- unlist(lapply(ratevec_files_short, extract_ratevec))

# Format for output
formatted_vector <- paste0("ratevectorp = {", paste(format(all_values, scientific = FALSE), collapse = ","), "};")

# Write to a file
writeLines(formatted_vector, file.path(ratevec_dir, "combined_OrderedLowtoHigh.ratevec"))