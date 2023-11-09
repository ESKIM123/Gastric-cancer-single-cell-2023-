# load libraries
library(dplyr)
library(plyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(Matrix)
library(ggeasy)
library(tidyr)
library(gridExtra)
library(harmony)

# Get data
filtered_integration <- readRDS(file="~/Single_cell_sample_merged/integrated_seurat_filtered.rds")
DimPlot(filtered_integration, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("Unsupervised clustering 1:30 / 2 (base)")

# Change dimensionality and visualize UMAP 

######################################################################################
for (i == 30) {
  filtered_integration <- RunUMAP(filtered_integration, reduction = "pca", dims = 1:i)
  filtered_integration <- FindNeighbors(filtered_integration, reduction = "pca", dims = 1:i)
  filtered_integration <- FindClusters(filtered_integration, resolution = 5)
  
  # plot the results and save as PNG file
  plot_title <- paste("Unsupervised clustering 1:", i, sep = "")
  png_1_i <- DimPlot(filtered_integration, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) +
    ggtitle(plot_title) +
    NoLegend() +
    FeaturePlot(filtered_integration, features = 'HSPB1', min.cutoff = '0', max.cutoff = '3', raster=FALSE)
  
  ggsave(paste0("png_1_", i, ".png"), png_1_i, width = 16, height = 7, units = "in")
}
######################################################################################

















filtered_integration_new <- subset(filtered_integration, idents = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
                                                                    20,21,22,24,25,26,29,32,33,34,35,36))

nrow(filtered_integration@meta.data)
nrow(filtered_integration_new@meta.data)

filtered_integration_new <- RunUMAP(filtered_integration_new, reduction = "pca", dims = 1:30)
filtered_integration_new <- FindNeighbors(filtered_integration_new, reduction = "pca", dims = 1:30)
filtered_integration_new <- FindClusters(filtered_integration_new, resolution = 1.5)
DimPlot(filtered_integration_new, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("Unsupervised clustering (filtered new)")


# Change dimensionality and visualize UMAP 
filtered_integration_new <- RunUMAP(filtered_integration_new, reduction = "pca", dims = 1:40) 
filtered_integration_new <- FindNeighbors(filtered_integration_new, reduction = "pca", dims = 1:40) 
filtered_integration_new <- FindClusters(filtered_integration_new, resolution = 5)
DimPlot(filtered_integration_new, reduction = 'umap', raster=FALSE, label = TRUE) 

filtered_integration_new2 <- subset(filtered_integration_new, idents = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
                                                                    20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,
                                                                    40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,
                                                                    60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,
                                                                    80,81,82,83))
nrow(filtered_integration@meta.data)
nrow(filtered_integration_new2@meta.data)

DimPlot(filtered_integration_new2, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Unsupervised clustering (filtered new2)")

saveRDS(filtered_integration_new2, file="~/Single_cell_sample_merged/integrated_seurat_filtered_new.rds")

