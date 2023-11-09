# load libraries
library(dplyr)
library(plyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(DoubletFinder)
library(Matrix)
library(ggeasy)
library(tidyr)
library(glmGamPoi)
library(gridExtra)

# Get data
integrated_seurat <- readRDS(file="~/Single_cell_sample_merged/integrated_seurat.rds")
DimPlot(integrated_seurat, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("Unsupervised clustering 1:30 (base)")

# Change dimensionality and visualize UMAP 
integrated_seurat <- RunUMAP(integrated_seurat, reduction = "pca", dims = 1:30)
integrated_seurat <- FindNeighbors(integrated_seurat, reduction = "pca", dims = 1:30)
integrated_seurat <- FindClusters(integrated_seurat, resolution = 1)
DimPlot(integrated_seurat, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("Unsupervised clustering 1:30")

cluster_features <- c('PTPRC', 'CD3E', 'CD8A', 'CD4', 'FOXP3', 'IL2RA', 'KLRD1', 'NKG7', 'KLRF1', 'MS4A1', 'CD19', 'MKI67', 'MZB1', 'JCHAIN', 
                             'CD14', 'CFP', 'APOBEC3A', 'LAMP3', 'TREM2', 'KIT', 'IL1RL1', 'DCN', 'FN1', 'COL6A2', 'RGS5', 'PLVAP', 'CDH1', 'EPCAM')
DotPlot(integrated_seurat, features = cluster_features)+ RotatedAxis()

################################################################################## Remove doublet clusters
integrated_seurat_filtered_1 <- subset(integrated_seurat, idents = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
                                                                   20,21,22,23,24,25,26,28,29,30,32,39))

integrated_seurat_filtered_1 <- FindNeighbors(integrated_seurat_filtered_1, reduction = "pca", dims = 1:30)
integrated_seurat_filtered_1 <- FindClusters(integrated_seurat_filtered_1, resolution = 1)
DimPlot(integrated_seurat_filtered_1, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("Unsupervised clustering (filtered_1)")

integrated_seurat_filtered_2 <- subset(integrated_seurat_filtered_1, idents = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
                                                                                20,21,22,23,24,25,26,27,29,30,31,32))

integrated_seurat_filtered_2 <- RunUMAP(integrated_seurat_filtered_2, reduction = "pca", dims = 1:30)
integrated_seurat_filtered_2 <- FindNeighbors(integrated_seurat_filtered_2, reduction = "pca", dims = 1:30)
integrated_seurat_filtered_2 <- FindClusters(integrated_seurat_filtered_2, resolution = 2)
DimPlot(integrated_seurat_filtered_2, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("Unsupervised clustering (filtered_2)")

nrow(integrated_seurat@meta.data)
nrow(integrated_seurat_filtered_1@meta.data)
nrow(integrated_seurat_filtered_2@meta.data)

saveRDS(integrated_seurat_filtered_2, file="~/Single_cell_sample_merged/integrated_seurat_filtered.rds")



