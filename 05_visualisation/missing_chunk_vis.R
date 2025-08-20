# Install and load required packages
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("patchwork")) install.packages("patchwork")

library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Create the missing data comparison table
missing_data_10k <- data.frame(
  Dataset = c("FX_160", "C9_128", "C9_256", "FX_640", "C9_512", "C9_1024"),
  Chunk_Size = rep(10000, 6),
  Missing_1_Chunk = c(3.00, 31.00, 77.00, 54.00, 22.00, 5.00),
  Missing_2plus_Chunks = c(1.00, 3.00, 4.00, 9.00, 60.00, 88.00),
  Complete_Reads = c(96.00, 66.00, 19.00, 37.00, 18.00, 7.00)
)

missing_data_custom <- data.frame(
  Dataset = c("FX_160", "C9_128", "C9_256", "FX_640", "C9_512", "C9_1024"),
  Chunk_Size = c(35000, 35000, 35000, 55000, 55000, 75000),
  Missing_1_Chunk = c(0.80, 0.70, 2.10, 2.60, 0.00, 1.40),
  Missing_2plus_Chunks = c(0.00, 0.10, 0.10, 0.20, 0.00, 0.00),
  Complete_Reads = c(99.20, 99.20, 97.80, 99.20, 100.00, 97.30)
)

# Create combined comparison table for display
comparison_table <- data.frame(
  Dataset = rep(c("FX_160", "C9_128", "C9_256", "FX_640", "C9_512", "C9_1024"), 2),
  Chunk_Size = c(rep("10,000", 6), c("35,000", "35,000", "35,000", "55,000", "55,000", "75,000")),
  Total_Missing = c(
    missing_data_10k$Missing_1_Chunk + missing_data_10k$Missing_2plus_Chunks,
    missing_data_custom$Missing_1_Chunk + missing_data_custom$Missing_2plus_Chunks
  ),
  Complete_Reads = c(missing_data_10k$Complete_Reads, missing_data_custom$Complete_Reads)
) %>%
  mutate(
    Strategy = ifelse(Chunk_Size == "10,000", "Standard 10K", "Custom"),
    Dataset = factor(Dataset, levels = c("FX_160", "C9_128", "C9_256", "FX_640", "C9_512", "C9_1024"))
  )

# Print the comparison table
cat("Table: Comparison of Missing Data Between Standard and Custom Chunk Sizes\n")
cat("Dataset\tChunk Size\tTotal Missing (%)\tComplete Reads (%)\n")
cat("-------\t----------\t-----------------\t------------------\n")
for(i in 1:nrow(comparison_table)) {
  cat(sprintf("%s\t%s\t\t%.1f\t\t\t%.1f\n", 
              comparison_table$Dataset[i], 
              comparison_table$Chunk_Size[i],
              comparison_table$Total_Missing[i],
              comparison_table$Complete_Reads[i]))
}

# Prepare data for visualization
sequence_order <- c("FX_160", "C9_128", "C9_256", "FX_640", "C9_512", "C9_1024")

prepare_missing_data <- function(data, condition) {
  data %>%
    mutate(Dataset = factor(Dataset, levels = sequence_order)) %>%
    pivot_longer(cols = c("Missing_1_Chunk", "Missing_2plus_Chunks", "Complete_Reads"), 
                 names_to = "Category", 
                 values_to = "Percentage") %>%
    mutate(
      Condition = condition,
      Category = case_when(
        Category == "Complete_Reads" ~ "Complete Reads",
        Category == "Missing_1_Chunk" ~ "Missing 1 Chunk", 
        Category == "Missing_2plus_Chunks" ~ "Missing 2+ Chunks"
      ),
      Category = factor(Category, levels = c("Complete Reads", "Missing 1 Chunk", "Missing 2+ Chunks"))
    )
}

# Prepare both datasets
data_10k_missing <- prepare_missing_data(missing_data_10k, "Standard 10K Chunks")
data_custom_missing <- prepare_missing_data(missing_data_custom, "Custom Chunk Sizes")

# Combine data
combined_missing_data <- bind_rows(data_10k_missing, data_custom_missing) %>%
  mutate(Condition = factor(Condition, levels = c("Standard 10K Chunks", "Custom Chunk Sizes")))

# Create the visualization
p_missing <- ggplot(combined_missing_data, aes(x = Dataset, y = Percentage, fill = Category)) +
  geom_col(position = "stack", alpha = 0.9, width = 0.7) +
  facet_wrap(~Condition, ncol = 2) +
  scale_fill_manual(values = c("Complete Reads" = "#2E8B57",      # Sea green
                               "Missing 1 Chunk" = "#DAA520",     # Goldenrod  
                               "Missing 2+ Chunks" = "#DC143C"),  # Crimson
                    name = "Data Status") +
  labs(
    title = "Impact of Chunk Size Strategy on Training Data Completeness",
    subtitle = "Custom chunk sizes dramatically reduce missing data across all datasets",
    x = "Dataset (ordered by increasing sequence length)",
    y = "Percentage of Training Data (%)",
    caption = "Missing chunks likely correspond to repeat expansion regions removed by Bonito's internal filtering"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30", margin = margin(b = 20)),
    plot.caption = element_text(size = 9, color = "gray50", hjust = 0.5, margin = margin(t = 15)),
    strip.text = element_text(size = 12, face = "bold", color = "gray20", 
                              margin = margin(t = 8, b = 8)),
    strip.background = element_rect(fill = "gray95", color = "gray70"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", size = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.8, "cm"),
    plot.margin = margin(20, 20, 20, 20),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), 
                     breaks = seq(0, 100, 20),
                     labels = function(x) paste0(x, "%"))

# Save the plot
ggsave("missing_data_comparison.png", plot = p_missing, 
       width = 14, height = 8, dpi = 300, bg = "white")

print("Plot saved as 'missing_data_comparison.png'")
print("Comparison table printed above")

# Display the plot
print(p_missing)