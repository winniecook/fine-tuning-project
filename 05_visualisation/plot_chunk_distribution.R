# Install and load required packages
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("viridis")) install.packages("viridis")
if (!require("patchwork")) install.packages("patchwork")

library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(patchwork)

# Create the data
chunk_10k <- data.frame(
  Dataset = c("FX_160", "C9_128", "C9_256", "FX_640", "C9_512", "C9_1024"),
  Chunk_Size = rep(10000, 6),
  Reads_1_Chunk = c(4.00, 0.00, 0.00, 0.00, 0.00, 0.00),
  Reads_2_Chunks = c(77.77, 81.00, 4.00, 3.00, 0.00, 0.00),
  Reads_3_Chunks = c(14.00, 15.00, 86.00, 83.00, 4.00, 12.00),
  Reads_4_Chunks = c(5.00, 4.00, 10.00, 14.00, 96.00, 88.00)
)

chunk_custom <- data.frame(
  Dataset = c("FX_160", "C9_128", "C9_256", "FX_640", "C9_512", "C9_1024"),
  Chunk_Size = c(35000, 35000, 35000, 55000, 55000, 75000),
  Reads_1_Chunk = c(98.80, 98.70, 95.90, 95.70, 100.00, 90.00),
  Reads_2_Chunks = c(1.20, 1.20, 3.90, 4.20, 0.00, 10.00),
  Reads_3_Chunks = c(0.00, 0.10, 0.10, 0.00, 0.00, 0.00),
  Reads_4_Chunks = c(0.00, 0.00, 0.00, 0.00, 0.00, 0.00)
)

# Add sequence length order for proper ordering
sequence_order <- c("FX_160", "C9_128", "C9_256", "FX_640", "C9_512", "C9_1024")

# Transform data for plotting
prepare_data <- function(data, condition) {
  data %>%
    mutate(Dataset = factor(Dataset, levels = sequence_order)) %>%
    pivot_longer(cols = starts_with("Reads_"), 
                 names_to = "Chunk_Category", 
                 values_to = "Percentage") %>%
    mutate(
      Condition = condition,
      Chunk_Category = case_when(
        Chunk_Category == "Reads_1_Chunk" ~ "1 Chunk",
        Chunk_Category == "Reads_2_Chunks" ~ "2 Chunks",
        Chunk_Category == "Reads_3_Chunks" ~ "3 Chunks",
        Chunk_Category == "Reads_4_Chunks" ~ "4+ Chunks"
      ),
      Chunk_Category = factor(Chunk_Category, levels = c("1 Chunk", "2 Chunks", "3 Chunks", "4+ Chunks"))
    )
}

# Prepare both datasets
data_10k <- prepare_data(chunk_10k, "Fixed 10K Chunks")
data_custom <- prepare_data(chunk_custom, "Custom Chunk Sizes")

# Combine data
combined_data <- bind_rows(data_10k, data_custom) %>%
  mutate(Condition = factor(Condition, levels = c("Fixed 10K Chunks", "Custom Chunk Sizes")))

# Create the visualization
p <- ggplot(combined_data, aes(x = Dataset, y = Percentage, fill = Chunk_Category)) +
  geom_col(position = "stack", alpha = 0.9, width = 0.7) +
  facet_wrap(~Condition, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = c("1 Chunk" = "#2E8B57",      # Sea green
                                "2 Chunks" = "#4682B4",    # Steel blue  
                                "3 Chunks" = "#DAA520",    # Goldenrod
                                "4+ Chunks" = "#DC143C"),  # Crimson
                     name = "Read Distribution") +
  labs(
    title = "Impact of Chunk Size Strategy on Read Fragmentation",
    subtitle = "Fixed 10K chunks create more fragmentation as sequence length increases",
    x = "Dataset (ordered by increasing sequence length)",
    y = "Percentage of Reads (%)",
    caption = "Data shows percentage of reads requiring different numbers of chunks for complete coverage"
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
ggsave("chunking_comparison.png", plot = p, 
       width = 14, height = 8, dpi = 300, bg = "white")

print("Plot saved as 'chunking_comparison.png'")

# Display the plot
print(p)