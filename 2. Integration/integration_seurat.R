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
merged_seurat<- readRDS(file="~/Single_cell_sample_merged/merged_seurat_singlet.rds")

# Perform integration to correct for batch effects ------
obj.list <- SplitObject(merged_seurat, split.by = 'Gem')
obj.list <- lapply(X = obj.list, FUN = SCTransform)

# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features)

obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# Find integration anchors 
#obj.anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT", anchor.features = features)
obj.anchors <- FindIntegrationAnchors(object.list = obj.list, reference=c(1),
                                  normalization.method="SCT",
                                  anchor.features=features,
                                  verbose=T,
                                  reduction="rpca")

# integrate data
obj.combined.sct <- IntegrateData(anchorset = obj.anchors, normalization.method = "SCT")

# run PCA and UMAP and visualize integrated data
obj.combined.sct <- RunPCA(obj.combined.sct, verbose = FALSE)

obj.combined.sct <- RunUMAP(obj.combined.sct, reduction = "pca", dims = 1:30)
obj.combined.sct <- FindNeighbors(obj.combined.sct, reduction = "pca", dims = 1:30)
obj.combined.sct <- FindClusters(obj.combined.sct, resolution = 0.3)


DimPlot(obj.combined.sct, reduction = 'umap', group.by = 'seurat_clusters') + ggtitle("Unsupervised clustering")
DimPlot(obj.combined.sct, reduction = 'umap', group.by = 'Gem')
DimPlot(obj.combined.sct, reduction = 'umap', group.by = 'Type')
DimPlot(obj.combined.sct, reduction = 'umap', group.by = 'Patient_ID')

saveRDS(obj.combined.sct, file="~/Single_cell_sample_merged/integrated_seurat.rds")





