library(tidyverse)

# Define the folder with the parameter files
param_folder <- "~/Documents/SignalNoiseBiasPaper/Archiving/Percomorpha/MathematicaInputs"  # Replace with your actual path

# List all relevant files
param_files <- list.files(
  param_folder, 
  pattern = "fasta\\.params\\.txt$", 
  full.names = TRUE
)

# Function to parse a single file
parse_param_file <- function(filepath) {
  lines <- readLines(filepath)
  
  # Remove semicolons and empty lines
  lines <- lines[lines != ""]
  lines <- gsub(";", "", lines)
  
  # Split on "="
  key_values <- str_split_fixed(lines, "=", 2)
  
  # Create named vector
  values <- as.numeric(key_values[, 2])
  names(values) <- key_values[, 1]
  
  # Convert to tibble with one row
  as_tibble_row(values)
}

# Apply to all files and collect
param_data <- map_dfr(param_files, parse_param_file, .id = "file")

# Optional: use base filename only (remove path)
param_data <- param_data %>%
  mutate(file = basename(file)) %>%
  relocate(file)

# View the result
print(param_data)

pi_long <- param_data %>%
  pivot_longer(-file, names_to = "parameter", values_to = "value") %>%
  filter(str_detect(parameter, "^pi[ATCG][p1234]$")) %>%
  mutate(
    base = str_extract(parameter, "[ATCG]"),
    group = str_extract(parameter, "(?<=pi[ATCG])\\w")  # extract group after base
  ) %>%
  mutate(group = factor(group, levels = c("p", "1", "2", "3", "4")),
         base = factor(base, levels = c("A", "T", "C", "G")))

# Step 2: Summarize (mean and SD)
pi_summary <- pi_long %>%
  group_by(group, base) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE),
    .groups = "drop"
  )

# Step 3: Plot with facets — one for each group
ggplot(pi_summary, aes(x = base, y = mean_value, fill = base)) +
  geom_col(width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value),
    width = 0.2
  ) +
  facet_wrap(~ group, nrow = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Base Composition by pi Group (Mean ± SD)",
    x = "Base",
    y = "Mean pi Value"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(size = 12)
  )
