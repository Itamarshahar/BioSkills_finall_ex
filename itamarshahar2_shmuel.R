# Clear the environment
rm(list = ls())

# Load required libraries
library(googlesheets4)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(ggExtra)
library(tidyverse)
library(pheatmap)
library(collapsibleTree)
#install.packages("collapsibleTree")
#input_url <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/2nd_year/General/R/final/compare_betas_by_l1.csv"
input_url <- "/Users/itamar_shahar/Downloads/compare_betas_by_l1 (6).csv"

# Fetch data from Google Sheets
sheet_input <- read.csv(input_url)
head(sheet_input)
################################################################################ 
## generating a meta data 
################################################################################ 

split_row_names <- strsplit(colnames(sheet_input), "\\.")

familiy_name <- character(length(split_row_names))

for (k in 1:28) {
  # Extract the "bodypart" value for the current row and store it in the vector
  familiy_name[k] <- split_row_names[[k]][1]
}

gene_name <- character(length(split_row_names))

for (k in 1:28) {
  # Extract the last two elements of the split row name and paste them together
  gene_name[k] <- paste(tail(split_row_names[[k]], 2), collapse = ".")
}



meta_data <- data.frame(
  row.names = colnames(sheet_input),
  family_name = familiy_name,  
  gene_name = gene_name
)

meta_data_annotaion <- data.frame(
  row.names = gene_name[-1],
  Area = familiy_name[-1]
)

################################################################################ 
## generating the colors for the heatmap 
################################################################################ 

tmp_annotation_color <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00")

names(tmp_annotation_color) <- unique(familiy_name)

my_annotation_color <- list(familiy_n = tmp_annotation_color)


#hput))
#rownames(sheet_input) <- colnames(sheet_input)[-1]

rownames(sheet_input) <- gene_name[-1]
sheet_input <- sheet_input[, -1]

colnames(sheet_input) <- gene_name[-1]

#A <- 6
callback <- function(hc, mat) {
  
  print(hc$order)# <-  seq(1, 27)
  #print(hc$order[3])
}

color_gradient <- colorRampPalette(c("white","blue"))

#pheatmap(test, clustering_callback = callback)
title <- "The Methylation Similarity Between Different Cells at Human Body"#  levels acrros Are epithelial cells across the digestive system have different methylation patterns? Why? 

heatmap_plot <- pheatmap(
                        sheet_input,
                         r_rows = FALSE,
                         cluster_cols = FALSE,
                         #cluster_row = FALSE,
                         angle_col = 90,
                         #show_rownames = FALSE,
                         #show_colnames = FALSE,
                        #color = colorRampPalette(c("blue", , "red"))(50)),
                        color = color_gradient(100),

                        #custom_colors <- c("#fef0d9", "#fdd49e", "#fdbb84", "#fc8d59", "#ef6548", "#d7301f", "#990000")
                        #custom_colors <- c("#fef0d9", "#fdd49e", "#fdbb84", "#fc8d59", "#ef6548", "#d7301f", "#990000"),
                         annotation_row = meta_data_annotaion,
                         #annotation_row= meta_data_n,
                         #breaks = my_breaks,
                         cellwidth = 10,
                         cellheight = 10,
                         fontsize_row = 8,
                         fontsize_col = 8,
                         annotation_colors = my_annotation_color,
                        #clustering_distance_rows = "manhattan",
                        #clustering_method = ""
                         #Rowv = as.dendrogram(hclust(dist(sheet_input))),
                      main =title,
                      fontsize = 8,
                        
                        #clustering_callback = callback,
                         #row_order = row.names(sheet_input)
                         )
output_file <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/2nd_year/General/R/final/heatmap.pdf"
output_file <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/2nd_year/General/R/final"
ggsave(filename = output_file, 
       plot = heatmap_plot, 
       width = 15, 
       height = 8)

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
