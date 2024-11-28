library(dplyr)
library(tidyverse)
library(ggplot2)
#library(GEOquery)
library(pheatmap)
library(readxl)
library(ComplexHeatmap)
#library(InteractiveComplexHeatmap)
library(circlize)

library(dendextend)


#filename <-'/home/damian/Desktop/R/datatest/va_c.csv'
filename2<-'/home/damian/Desktop/R/datatest/cvc_one.csv'


data <- read.table(filename2, sep = "\t", quote = "", stringsAsFactors = FALSE, header = TRUE)
duplicates <- data[duplicated(data[[1]]), ] # Assuming the first column is used as row names
print(duplicates)






# Read the data into a data.frame
my_data1 <- read.table(filename2, sep="\t", quote="", stringsAsFactors=FALSE,header=TRUE,row.names=1)
#my_data <- read.table(filename, sep="\t", quote="", stringsAsFactors=FALSE,header=TRUE)


dim(my_data) # (rows columns)

nrow(my_data) # rows: locations (bins) in genome
ncol(my_data) # columns: cells

my_matrix1 <- as.matrix(my_data1)
#my_matrix <- as.matrix(my_data[c(1:26)  ,c(1:73)]) # [all rows, columns 2-25]

scaled_matrix <- t(scale(t(my_matrix1))) 

# Check the range of values in the data matrix
range(my_matrix, na.rm = TRUE)

# Define the new color function based on the actual range of data
# Assuming the range from the above check is [0, 10], you can modify it like this:
col_fun = colorRamp2(c(0, 5, 10), c("white", "blue", "red"))

# Save gene column for annotating the heatmap later
gene_info <- data.frame(Gene = my_data$Gene)
gene_info



# Default parameters
Heatmap(scaled_matrix )


# Flip rows and columns around
scaled_matrix  <- t(scaled_matrix )  # "transpose"
Heatmap(scaled_matrix )


# Keep genome bins in order, not clustered
Heatmap(scaled_matrix , cluster_columns=FALSE)


fontsize <- 0.6

# Put cell labels on the left side
Heatmap(scaled_matrix , 
        name = "With Cluster",
        cluster_columns=TRUE,
        row_names_side = "left", 
        show_row_dend = TRUE,
        row_names_gp=gpar(cex=fontsize))
# Make the dendrogram wider
Heatmap(scaled_matrix , 
        name = "Sample1",
        cluster_columns=FALSE,
        row_names_side = "left", 
        show_row_dend = TRUE,
        row_names_gp=gpar(cex=fontsize),
        row_dend_width = unit(3, "cm"))

Heatmap(scaled_matrix, 
        cluster_columns = FALSE,
        row_names_side = "left", 
        show_row_dend = TRUE,
        row_names_gp = gpar(cex = fontsize),  # Row font size
        column_names_gp = gpar(cex = 0.5))


Heatmap(scaled_matrix, 
        name = "Cluster",
        cluster_columns = TRUE,
        row_names_side = "left", 
        show_row_dend = TRUE,
        row_names_gp = gpar(cex = fontsize),  # Row font size
        column_names_gp = gpar(cex = 0.5))

# Different distance calculation methods
# "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"
# euclidean is the default

# Different clustering methods
# "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)

# Watch the dendrogram and heatmap change when we change the methods
Heatmap(scaled_matrix , 
        cluster_columns = TRUE,                  # Disable column clustering
        row_names_side = "left",                  # Display row names on the left
        row_names_gp = gpar(cex = fontsize),      # Adjust font size for row names
        clustering_distance_rows = "maximum",    # Use 'maximum' distance metric for clustering rows
        clustering_method_rows = "ward.D"        # Use 'ward.D' method for hierarchical clustering of rows
)

# Need to build dendrogram first so we can use it for the color_brances() function
# 1. calculate distances (method="maximum")
# 2. cluster (method="ward.D")
# Create a dendrogram with colored branches

row_dend <- as.dendrogram(hclust(dist(my_matrix), method = "ward.D"))
row_dend <- color_branches(row_dend, k = 3)  # Color branches into 3 groups

# Generate the heatmap
Heatmap(
  my_matrix,
  cluster_columns = FALSE,                  # Disable column clustering
  cluster_rows = row_dend,                  # Use the modified dendrogram for rows
  row_names_side = "left",                  # Display row names on the left
  row_names_gp = gpar(cex = fontsize),      # Adjust font size for row names
  row_dend_width = unit(3, "cm")            # Set the width of the row dendrogram
)


# We can split the heatmap into clusters

Heatmap(
  my_matrix, 
  cluster_columns = FALSE,                # Disable column clustering
  row_names_side = "left",                # Display row names on the left
  show_row_dend = TRUE,                   # Show row dendrogram
  row_names_gp = gpar(cex = fontsize),    # Adjust font size for row names
  row_dend_width = unit(3, "cm"),         # Set the width of the row dendrogram
  clustering_distance_rows = "maximum",  # Use 'maximum' distance metric for row clustering
  clustering_method_rows = "ward.D",     # Use 'ward.D' for hierarchical clustering
  km = 2                                 # Number of clusters for k-means clustering
)


gene_info
gene_info.colors <- c(rep(c("black","white"),13),"red")
gene_info.colors

names(gene_info.colors) <- paste("chr",c(seq(1,26),"X"),sep="")
gene_info.colors

Heatmap(my_matrix, 
        cluster_columns=FALSE,
        row_names_side = "left", 
        show_row_dend = TRUE, 
        row_names_gp=gpar(cex=fontsize),
        row_dend_width = unit(3, "cm"),
        clustering_distance_rows ="maximum",
        clustering_method_rows = "ward.D",
        km=2, # number of clusters you want
        bottom_annotation = HeatmapAnnotation(gene_info,col = list(chrom=gene_info.colors),show_legend=FALSE)
)
