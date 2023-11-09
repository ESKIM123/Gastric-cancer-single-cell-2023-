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
cancer.subset <- readRDS(file="D:/DATA/Gastric cancer/immune_cell_clinical_data_modify.rds")
DimPlot(cancer.subset, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Initial")

view(cancer.subset@meta.data)

# Subset 
immune_annotation <- subset(x = cancer.subset, Type == 'Tumor' & Patient_ID_SPC != FALSE &
                              Patient_ID_SPC != 8 & Patient_ID_SPC != 18 & Patient_ID_SPC != 52)


##################################### Dataframe for Outcome
md.immune_annotation<- immune_annotation@meta.data %>% as.data.table

length(which(md.immune_annotation$Patient_ID==68 & md.immune_annotation$Type=='Tumor'& md.immune_annotation$seurat_clusters==1))

# Vector of sample numbers
sample_numbers <- unique(immune_annotation$Patient_ID_SPC)

# Vector of cluster names
cluster_vec <- unique(md.immune_annotation$clusters)
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
      subset_length <- length(which(md.immune_annotation$Patient_ID_SPC == number & md.immune_annotation$Type == type & md.immune_annotation$clusters == cluster))
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

##################################### Clustering

new2_df <- new_df[, !(names(new_df) %in% c("Proliferating B cells", "Proliferating T cells", "Mast cells", "MAIT cells"))]

unique_cluster_vec2 <- c("CD8 T cells",  "CD4 T cells", "B cells", "Myeloid cells", "NK cells", "Regulatory CD4 T cells","Plasma cells", "Gamma-delta T cells")
colnames(new_df)


# un-scaling
new_df_unscaled <- new2_df


# Set row names to be the sample_numbers
row.names(new_df_unscaled) <- new_df_unscaled$sample_numbers

new_df_unscaled$sample_numbers <- NULL

# Create a custom pastel color palette
pastel_palette <- brewer.pal(n = 9, name = "YlGnBu")

# Create the heatmap with a custom color palette
heatmap.2(as.matrix(new_df_unscaled), dendrogram = "row", scale = "none", key = TRUE, col = pastel_palette, margins=c(11,11))



# Min-max scaling
new_df_scaled <- new2_df

for (col in unique_cluster_vec2) {
  col_data <- new_df_scaled[, col]
  new_df_scaled[, col] <- (col_data - min(col_data)) / (max(col_data) - min(col_data))
}


# Set row names to be the sample_numbers
row.names(new_df_scaled) <- new_df_scaled$sample_numbers
new_df_scaled$sample_numbers <- NULL

# Create a custom pastel color palette
pastel_palette <- brewer.pal(n = 9, name = "YlGnBu")

# Create the heatmap with a custom color palette
heatmap.2(as.matrix(new_df_scaled), dendrogram = "row", scale = "none", key = TRUE, col = pastel_palette, margins=c(15,3),
          density.info="none",  trace = "none",
          cexRow = 1,
          srtCol = 45, 
          sepCol = 0.5)



# z-score scaling
for (col in unique_cluster_vec2) {
  col_data <- new_df_scaled[, col]
  new_df_scaled[, col] <- (col_data - mean(col_data)) / sd(col_data)
}

# Create a custom color palette function
custom_color_palette <- colorRampPalette(c("blue", "white", "red"))

# Create the heatmap with the custom color palette
heatmap.2(as.matrix(new_df_scaled), dendrogram = "row", scale = "none", key = TRUE, col = custom_color_palette(100),
          margins = c(15, 3), density.info = "none", trace = "none",
          cexRow = 1, srtCol = 45, sepCol = 0.5)


############## do_AlluvialPlot ####

view(immune_annotation@meta.data)

immune_annotation[["Immune_type"]]  <- immune_annotation[["Patient_ID_SPC"]]


immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 75] <- "Treg_myeloid"
immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 20] <- "Treg_myeloid"
immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 15] <- "Treg_myeloid"
immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 50] <- "Treg_myeloid"
immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 73] <- "Treg_myeloid"
immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 29] <- "Treg_myeloid"
immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 40] <- "Treg_myeloid"
immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 22] <- "Treg_myeloid"
immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 60] <- "Treg_myeloid"
immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 62] <- "Treg_myeloid"
immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 7] <- "Treg_myeloid"


immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 28] <- "B_Plasma"
immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 5] <- "B_Plasma"
immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 27] <- "B_Plasma"
immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 36] <- "B_Plasma"
immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 64] <- "B_Plasma"
immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 4] <- "B_Plasma"
immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 54] <- "B_Plasma"
immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 68] <- "B_Plasma"
immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 65] <- "B_Plasma"
immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 66] <- "B_Plasma"


immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 23] <- "Myeloid_cell_type"


immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 59] <- "CD8_T_cell"
immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 56] <- "CD8_T_cell"
immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 51] <- "CD8_T_cell"
immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 41] <- "CD8_T_cell"
immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 44] <- "CD8_T_cell"
immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 71] <- "CD8_T_cell"
immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 35] <- "CD8_T_cell"
immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 16] <- "CD8_T_cell"
immune_annotation@meta.data["Immune_type"][immune_annotation@meta.data["Immune_type"] == 53] <- "CD8_T_cell"




do_AlluvialPlot(sample = immune_annotation,
                first_group = "Patient_ID_SPC",
                last_group = "Immune_type",
                middle_groups= "TCGA",
                plot.caption="NA",
                font.size = 5)

do_AlluvialPlot(sample = immune_annotation,
                first_group = "TCGA",
                last_group = "Immune_type")

