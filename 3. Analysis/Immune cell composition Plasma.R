# load libraries
library(dplyr)
library(ggpubr)
library(plyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(Matrix)
library(ggeasy)
library(tidyr)
library(grid)
library(gridExtra)
library(glmGamPoi)
library(cowplot)
library(DelayedArray)
library(UCell)
library(Nebulosa)
library(SCpubr)
library(viridis)
library(clusterProfiler)
library(org.Hs.eg.db)
library(readxl)
library(diffcyt)
library(ConsensusClusterPlus)
library(FlowSOM)
library(CATALYST)
library(flowCore)
library(ggplot2)
library(SingleCellExperiment)
library(dittoSeq)
library(scater)
library(patchwork)
library(cowplot)
library(viridis)
library(speckle)
library(limma)
library(ggplot2)
library(edgeR)
library(pheatmap)
library(gt)
library(data.table)
library(magrittr)
library(GOplot)
library(enrichR)
library(decoupleR)
library(OmnipathR)
library(dorothea)
library(fgsea)
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyverse)
library(data.table)
library(presto)
library(ggrepel)
library(gplots)
library(RColorBrewer)


# Get data ----
Plasma_subset_anno <- readRDS(file="D:/DATA/Gastric cancer/GS_subsets/Plasma_subset.rds")
DimPlot(Plasma_subset_anno, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Initial")

view(Plasma_subset_anno@meta.data)

# Subset 
immune_annotation <- subset(x = Plasma_subset_anno, Type == 'Tumor' & Patient_ID_SPC != FALSE)


##################################### Dataframe for Outcome
md.immune_annotation<- immune_annotation@meta.data %>% as.data.table

length(which(md.immune_annotation$Patient_ID==68 & md.immune_annotation$Type=='Tumor'& md.immune_annotation$seurat_clusters==1))

# Vector of sample numbers
sample_numbers <- unique(immune_annotation$Patient_ID_SPC)

# Vector of cluster names
cluster_vec <- unique(md.immune_annotation$seurat_clusters)
type_vec <- c("Tumor", "Normal")


# Create an empty data frame to store the results
output_df <- data.frame(sample_numbers = numeric(),
                        cluster_vec = character(),
                        type_vec = character(),
                        count = numeric(),
                        total_cell = numeric(),
                        subset_fraction = numeric(),
                        stringsAsFactors = FALSE)

# Loop over the variables and add rows to the data frame
for (type in type_vec) {
  for (cluster in cluster_vec) {
    for (number in sample_numbers) {
      # Calculate the length of the subset
      subset_length <- length(which(md.immune_annotation$Patient_ID_SPC == number & md.immune_annotation$Type == type & md.immune_annotation$seurat_clusters == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.immune_annotation$Patient_ID_SPC == number & md.immune_annotation$Type == type))
      subset_fraction <- subset_length / total_count
      # Add a row to the data frame
      new_row <- data.frame(sample_numbers = number,
                            cluster_vec = cluster,
                            type_vec = type,
                            count = subset_length,
                            total_cell = total_count,
                            subset_fraction = subset_fraction,
                            stringsAsFactors = FALSE)
      output_df <- rbind(output_df, new_row)
    }
  }
}


output_df


output_df_tumor <- subset(output_df, type_vec %in% "Tumor")

# Get unique values of sample_numbers and cluster_vec
unique_sample_numbers <- unique(output_df_tumor$sample_numbers)
unique_cluster_vec <- unique(output_df_tumor$cluster_vec)

# Create an empty dataframe with rows and columns based on unique values
new_df <- data.frame(
  sample_numbers = unique_sample_numbers,
  stringsAsFactors = FALSE
)

# Add columns for unique cluster_vec values
for (cluster in unique_cluster_vec) {
  new_df[, cluster] <- 0  # Initialize all values to 0
}

# Fill in the values from output_df_tumor
for (i in 1:nrow(output_df_tumor)) {
  row_index <- match(output_df_tumor$sample_numbers[i], new_df$sample_numbers)
  col_name <- output_df_tumor$cluster_vec[i]
  new_df[row_index, col_name] <- output_df_tumor$subset_fraction[i]
}



# un-scaling
new_df_unscaled <- new_df


# Set row names to be the sample_numbers
row.names(new_df_unscaled) <- new_df_unscaled$sample_numbers
new_df_unscaled$sample_numbers <- NULL

# Create a custom pastel color palette
pastel_palette <- brewer.pal(n = 9, name = "YlGnBu")

# Create the heatmap with a custom color palette
heatmap.2(as.matrix(new_df_unscaled), dendrogram = "row", scale = "none", key = TRUE, col = pastel_palette)



# Min-max scaling
new_df_scaled <- new_df
for (col in unique_cluster_vec) {
  col_data <- new_df_scaled[, col]
  new_df_scaled[, col] <- (col_data - min(col_data)) / (max(col_data) - min(col_data))
}

# Set row names to be the sample_numbers
row.names(new_df_scaled) <- new_df_scaled$sample_numbers
new_df_scaled$sample_numbers <- NULL

# Create a custom pastel color palette
pastel_palette <- brewer.pal(n = 9, name = "YlGnBu")

# Create the heatmap with a custom color palette
heatmap.2(as.matrix(new_df_scaled), dendrogram = "row", trace = "none", scale = "none", key = TRUE, density.info="none", col = pastel_palette, margins=c(7,7), 
          main = "Plasma T cells")