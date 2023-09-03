# Clear the environment
rm(list = ls())

# Load required libraries
library(googlesheets4)
library(ggplot2)
library(ggExtra)
library(ggrepel)
library(RColorBrewer)
library(tidyverse)
library(pheatmap)

input_url <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/2nd_year/General/R/final/compare_betas_by_l1.csv"

# Fetch data from Google Sheets
sheet_input <- read.csv(input_url)
my_annotation <-
sheet_input <- sheet_input[, -1]
heatmap_plot <- pheatmap(sheet_input,
                         cluster_rows = FALSE,
                         angle_col = 0,
                         show_rownames = FALSE,
                         color = my_palette) #,
                         #breaks = my_breaks,
                         #annotation_row = my_annotation,
                         #annotation_colors = my_annotation_color
                         #)                         

heatmap_plot

# Remove rows with null values
processed_sheets <- processed_sheets %>%
  filter(!is.na(GH), !is.na(log2FoldChange), !is.na(`Gene ID`))

# Filter data based on thresholds
processed_sheets_filtered <- processed_sheets %>%
  group_by(`Gene ID`) %>%
  filter(any(abs(log2FoldChange) > log_fold_change_th & padj < padj_th)) %>%
  ungroup() %>%
  arrange(GH) %>%
  as.data.frame()

# Create annotation data frame
my_annotation <- data.frame(NewGeneID = processed_sheets_filtered$`Gene ID`, GH = processed_sheets_filtered$GH)
my_annotation <- unique(my_annotation)
rownames(my_annotation) <- my_annotation$NewGeneID
my_annotation$NewGeneID <- NULL

# Prepare data for heatmap
sub_plot <- subset(processed_sheets_filtered, select = c("Gene ID", "log2FoldChange", "HMO"))
sub_plot_heatmap <- sub_plot %>%
  pivot_wider(names_from = "Gene ID", values_from = log2FoldChange, names_prefix = "") %>%
  column_to_rownames(var = "HMO")
sub_plot_heatmap_t <- t(sub_plot_heatmap)

# Create color palette
my_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# Set the breaks for the color palette
max_value <- abs(max(abs(sub_plot_heatmap_t)))
my_breaks <- c(seq(-max_value, 0, length.out = 51), seq(10^(-10), max_value, length.out = 50))

# Define annotation colors
tmp_annotation_color <- brewer.pal(length(gh_families), "Set1")
names(tmp_annotation_color) <- unique(my_annotation$GH)
my_annotation_color <- list(GH = tmp_annotation_color)

# Create the heatmap
heatmap_plot <- pheatmap(sub_plot_heatmap_t,
                         cluster_rows = FALSE,
                         angle_col = 0,
                         show_rownames = FALSE,
                         color = my_palette,
                         breaks = my_breaks,
                         annotation_row = my_annotation,
                         annotation_colors = my_annotation_color)

# Save the heatmap to a file
ggsave(filename = output_file, plot = heatmap_plot, width = 15, height = 8)
cat("Plot saved to:\n", output_file, "\n")
