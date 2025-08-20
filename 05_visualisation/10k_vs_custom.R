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

# Create accuracy data
accuracy_data <- data.frame(
  Dataset = rep(c("FX_160", "C9_128", "C9_256", "C9_512", "FX_640", "C9_1024"), 3),
  Model = rep(c("Baseline", "Fine-tuned 10K", "Fine-tuned Custom"), each = 6),
  Positive = c(
    # Baseline
    90.62, 91.00, 89.25, 84.75, 28.05, 83.38,
    # Fine-tuned 10K
    96.25, 96.50, 89.60, 84.75, 28.05, 82.01,
    # Fine-tuned Custom
    98.75, 98.70, 94.40, 92.77, 64.20, 78.71
  ),
  Negative = c(
    # Baseline
    74.37, 75.00, 75.00, 72.46, 25.03, 66.62,
    # Fine-tuned 10K
    90.23, 86.01, 74.95, 72.39, 25.09, 67.01,
    # Fine-tuned Custom
    99.23, 98.86, 99.00, 90.00, 64.20, 91.50
  )
)

# Create purity data
purity_data <- data.frame(
  Dataset = rep(c("FX_160", "C9_128", "C9_256", "C9_512", "FX_640", "C9_1024"), 3),
  Model = rep(c("Baseline", "Fine-tuned 10K", "Fine-tuned Custom"), each = 6),
  Positive = c(
    # Baseline
    98.24, 41.00, 33.23, 40.09, 96.00, 28.92,
    # Fine-tuned 10K
    99.95, 97.83, 89.00, 94.00, 99.00, 98.00,
    # Fine-tuned Custom
    100.00, 99.99, 100.00, 100.00, 100.00, 99.90
  ),
  Negative = c(
    # Baseline
    97.64, 41.00, 79.11, 25.28, 95.97, 17.92,
    # Fine-tuned 10K
    98.00, 98.00, 97.04, 98.14, 97.99, 89.45,
    # Fine-tuned Custom
    99.79, 99.86, 99.92, 99.41, 100.00, 99.83
  )
)

# Set factor levels for proper ordering
dataset_order <- c("FX_160", "C9_128", "C9_256", "C9_512", "FX_640", "C9_1024")
model_order <- c("Baseline", "Fine-tuned 10K", "Fine-tuned Custom")

accuracy_data$Dataset <- factor(accuracy_data$Dataset, levels = dataset_order)
accuracy_data$Model <- factor(accuracy_data$Model, levels = model_order)
purity_data$Dataset <- factor(purity_data$Dataset, levels = dataset_order)
purity_data$Model <- factor(purity_data$Model, levels = model_order)

# Reshape data for plotting
accuracy_long <- accuracy_data %>%
  pivot_longer(cols = c("Positive", "Negative"), names_to = "Class", values_to = "Accuracy")

purity_long <- purity_data %>%
  pivot_longer(cols = c("Positive", "Negative"), names_to = "Class", values_to = "Purity")

# Create color palette
model_colors <- c("Baseline" = "#8B4513", 
                  "Fine-tuned 10K" = "#4682B4", 
                  "Fine-tuned Custom" = "#2E8B57")

# Create accuracy plot
p_accuracy <- ggplot(accuracy_long, aes(x = Dataset, y = Accuracy, fill = Model)) +
  geom_col(position = position_dodge(width = 0.8), alpha = 0.9, width = 0.7) +
  facet_wrap(~Class, ncol = 2, labeller = labeller(Class = c("Positive" = "Forward Strand (+)", 
                                                             "Negative" = "Reverse Strand (-)"))) +
  scale_fill_manual(values = model_colors, name = "Model") +
  labs(
    title = "Repeat Length Accuracy: Impact of Chunk Size on Model Performance",
    subtitle = "Custom chunk sizes enable accuracy improvements across all datasets",
    x = "Dataset (ordered by increasing sequence length)",
    y = "Repeat Length Accuracy (%)",
    caption = "Results shown for forward and reverse strand sequences"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30", margin = margin(b = 15)),
    plot.caption = element_text(size = 9, color = "gray50", hjust = 0.5, margin = margin(t = 10)),
    strip.text = element_text(size = 11, face = "bold", color = "gray20"),
    strip.background = element_rect(fill = "gray95", color = "gray70"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", size = 0.3),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.6, "cm"),
    plot.margin = margin(15, 15, 15, 15),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), 
                     breaks = seq(0, 100, 20),
                     labels = function(x) paste0(x, "%")) +
  guides(fill = guide_legend(nrow = 1))

# Create purity plot
p_purity <- ggplot(purity_long, aes(x = Dataset, y = Purity, fill = Model)) +
  geom_col(position = position_dodge(width = 0.8), alpha = 0.9, width = 0.7) +
  facet_wrap(~Class, ncol = 2, labeller = labeller(Class = c("Positive" = "Forward Strand (+)", 
                                                             "Negative" = "Reverse Strand (-)"))) +
  scale_fill_manual(values = model_colors, name = "Model") +
  labs(
    title = "Prediction Purity: Impact of Chunk Size on Model Performance",
    subtitle = "Both fine-tuning approaches achieve substantial purity improvements",
    x = "Dataset (ordered by increasing sequence length)",
    y = "Prediction Purity (%)",
    caption = "Higher purity indicates more confident, reliable predictions"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30", margin = margin(b = 15)),
    plot.caption = element_text(size = 9, color = "gray50", hjust = 0.5, margin = margin(t = 10)),
    strip.text = element_text(size = 11, face = "bold", color = "gray20"),
    strip.background = element_rect(fill = "gray95", color = "gray70"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", size = 0.3),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.6, "cm"),
    plot.margin = margin(15, 15, 15, 15),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), 
                     breaks = seq(0, 100, 20),
                     labels = function(x) paste0(x, "%")) +
  guides(fill = guide_legend(nrow = 1))

# Save plots
ggsave("accuracy_comparison.png", plot = p_accuracy, 
       width = 12, height = 7, dpi = 300, bg = "white")

ggsave("purity_comparison.png", plot = p_purity, 
       width = 12, height = 7, dpi = 300, bg = "white")

print("Plots saved as 'accuracy_comparison.png' and 'purity_comparison.png'")

# Display plots
print(p_accuracy)
#print(p_purity)