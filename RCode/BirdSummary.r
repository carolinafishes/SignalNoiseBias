# Load libraries
library(tidyverse)

# Define folder with CSVs
data_folder <- "~/Documents/SignalNoiseBiasPaper/Archiving/Birds/Results"  # Replace with your actual folder path

# Get list of relevant CSV files (exclude any with "Check" in the name)
file_list <- list.files(data_folder, pattern = "\\.csv$", full.names = TRUE)
file_list <- file_list[!grepl("Check", file_list, ignore.case = TRUE)]

# Initialize counters
tally_noise_gt_signal <- 0
tally_bias_gt_signal <- 0
tally_noiseplusbias_gt_signal <- 0

# Loop through files
for (file in file_list) {
  df <- read.csv(file)
  
  if (nrow(df) == 0) next
  
  last_row <- tail(df, 1)
  mu_signal <- last_row$MuValueSignal
  mu_noise <- last_row$MuValueNoise
  mu_bias <- last_row$MuValueBias
  
  if (mu_noise > mu_signal) tally_noise_gt_signal <- tally_noise_gt_signal + 1
  if (mu_bias > mu_signal) tally_bias_gt_signal <- tally_bias_gt_signal + 1
  if ((mu_noise + mu_bias) > mu_signal) tally_noiseplusbias_gt_signal <- tally_noiseplusbias_gt_signal + 1
}

# Create summary dataframe
summary_df <- tibble(
  Condition = c("Noise > Signal", "Bias > Signal", "Noise + Bias > Signal"),
  Count = c(tally_noise_gt_signal, tally_bias_gt_signal, tally_noiseplusbias_gt_signal)
)
#convert to percent
summary_df[,2] <- summary_df[,2]/259*100

# Lollipop plot
ggplot(summary_df, aes(x = reorder(Condition, Count), y = Count)) +
  geom_col(fill = "steelblue", width = 0.6) +
  geom_text(aes(label = Count), hjust = -0.2, size = 5) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Percent of Loci where Signal is Overwhelmed by Noise or Bias",
    x = NULL,
    y = "Percent loci"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 11)
  ) +
  expand_limits(y = max(summary_df$Count) * 1.1)