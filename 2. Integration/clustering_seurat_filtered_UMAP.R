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

# Get data
filtered_integration <- readRDS(file="~/Single_cell_sample_merged/integrated_seurat_filtered.rds")

# Change dimensionality and visualize UMAP 
DimPlot(filtered_integration, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("Unsupervised clustering 1:30 / 2 (base)")

filtered_integration <- RunUMAP(filtered_integration, reduction = "pca", dims = 1:30)
filtered_integration <- FindNeighbors(filtered_integration, reduction = "pca", dims = 1:30)
filtered_integration <- FindClusters(filtered_integration, resolution = 1.5)
DimPlot(filtered_integration, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("Unsupervised clustering")

################################################################################
FeaturePlot(filtered_integration, features = 'PTPRC', min.cutoff = 'q1', max.cutoff = 'q10', raster=FALSE)
FeaturePlot(filtered_integration, features = 'CD3E', min.cutoff = 'q1', max.cutoff = 'q10', raster=FALSE)
FeaturePlot(filtered_integration, features = 'CD8A', min.cutoff = 'q1', max.cutoff = 'q10', raster=FALSE)
FeaturePlot(filtered_integration, features = 'CD4', min.cutoff = 'q1', max.cutoff = 'q10', raster=FALSE)
FeaturePlot(filtered_integration, features = 'FOXP3', min.cutoff = 'q1', max.cutoff = 'q10', raster=FALSE)
FeaturePlot(filtered_integration, features = 'KLRF1', min.cutoff = 'q1', max.cutoff = 'q10', raster=FALSE)

filtered_integration_annotation <- filtered_integration
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `29` = 'T/NK cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `10`  = 'T/NK cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `0`  = 'T/NK cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `5`  = 'T/NK cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `11`  = 'T/NK cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `19`  = 'T/NK cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `2`  = 'T/NK cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `32`  = 'T/NK cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `6`  = 'T/NK cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `14`  = 'T/NK cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `34`  = 'T/NK cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `13`  = 'T/NK cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `17`  = 'T/NK cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `4`  = 'T/NK cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `3`  = 'T/NK cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `20`  = 'T/NK cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `12`  = 'T/NK cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `33`  = 'T/NK cells')

################################################################################
FeaturePlot(filtered_integration, features = 'CD19', min.cutoff = 'q1', max.cutoff = 'q10', raster=FALSE)
FeaturePlot(filtered_integration, features = 'MS4A1', min.cutoff = 'q1', max.cutoff = 'q10', raster=FALSE)

filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `35`  = 'B cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `1`  = 'B cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `25`  = 'B cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `36`  = 'B cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `37`  = 'B cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `9`  = 'B cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `8`  = 'B cells')

################################################################################
FeaturePlot(filtered_integration, features = 'JCHAIN', min.cutoff = 'q1', max.cutoff = 'q10', raster=FALSE)
FeaturePlot(filtered_integration, features = 'MZB1', min.cutoff = 'q1', max.cutoff = 'q10', raster=FALSE)

filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `24`  = 'Plasma cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `7`  = 'Plasma cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `22`  = 'Plasma cells')

################################################################################
FeaturePlot(filtered_integration, features = 'KIT', min.cutoff = 'q50', max.cutoff = 'q90', raster=FALSE)
FeaturePlot(filtered_integration, features = 'IL1RL1', min.cutoff = 'q50', max.cutoff = 'q90', raster=FALSE)

filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `21`  = 'Mast cells')

################################################################################
FeaturePlot(filtered_integration, features = 'CDH1', min.cutoff = 'q50', max.cutoff = 'q90', raster=FALSE)
FeaturePlot(filtered_integration, features = 'EPCAM', min.cutoff = 'q50', max.cutoff = 'q90', raster=FALSE)
FeaturePlot(filtered_integration, features = 'KRT18', min.cutoff = 'q50', max.cutoff = 'q90', raster=FALSE)

filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `31`  = 'Epithelial cells')


################################################################################
FeaturePlot(filtered_integration, features = 'PECAM1', min.cutoff = 'q50', max.cutoff = 'q90', raster=FALSE)
FeaturePlot(filtered_integration, features = 'VWF', min.cutoff = 'q50', max.cutoff = 'q90', raster=FALSE)
FeaturePlot(filtered_integration, features = 'PLVAP', min.cutoff = 'q85', max.cutoff = 'q90', raster=FALSE)

filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `27`  = 'Endothelial cells')

################################################################################
FeaturePlot(filtered_integration, features = 'DCN', min.cutoff = 'q10', max.cutoff = 'q90', raster=FALSE)
FeaturePlot(filtered_integration, features = 'FN1', min.cutoff = 'q10', max.cutoff = 'q90', raster=FALSE)

FeaturePlot(filtered_integration, features = 'ACTA2', min.cutoff = 'q50', max.cutoff = 'q99', raster=FALSE)
FeaturePlot(filtered_integration, features = 'PDGFRB', min.cutoff = 'q50', max.cutoff = 'q99', raster=FALSE)
FeaturePlot(filtered_integration, features = 'RGS5', min.cutoff = 'q85', max.cutoff = 'q99', raster=FALSE)

filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `28`  = 'Fibroblasts')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `23`  = 'Fibroblasts')

filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `30`  = 'Pericytes')

################################################################################
FeaturePlot(filtered_integration, features = 'CD14', min.cutoff = 'q50', max.cutoff = 'q99', raster=FALSE)
FeaturePlot(filtered_integration, features = 'APOBEC3A', min.cutoff = 'q50', max.cutoff = 'q99', raster=FALSE)
FeaturePlot(filtered_integration, features = 'MRC1', min.cutoff = 'q1', max.cutoff = 'q99', raster=FALSE)
FeaturePlot(filtered_integration, features = 'CD209', min.cutoff = 'q1', max.cutoff = 'q99', raster=FALSE)
FeaturePlot(filtered_integration, features = 'FCGR1A', min.cutoff = 'q1', max.cutoff = 'q99', raster=FALSE)

filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `15`  = 'Myeloid cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `26`  = 'Myeloid cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `18`  = 'Myeloid cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `16`  = 'Myeloid cells')


Idents(filtered_integration_annotation)
view(filtered_integration_annotation@meta.data)

################################################################################# Visualization
DimPlot(filtered_integration_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Unsupervised clustering")
DimPlot(filtered_integration_annotation, reduction = 'umap', raster=FALSE, label = FALSE) + NoLegend() + NoAxes() 

DimPlot(filtered_integration_annotation, reduction = 'umap', group.by = 'Type', raster=FALSE, label = FALSE) + NoLegend() + NoAxes() + ggtitle(NULL)
DimPlot(filtered_integration_annotation, reduction = 'umap', group.by = 'Patient_ID', raster=FALSE, label = FALSE) + NoLegend() + NoAxes() + ggtitle(NULL)


filtered_integration_annotation_sum <- filtered_integration_annotation
Idents(filtered_integration_annotation_sum)

filtered_integration_annotation_sum <- RenameIdents(filtered_integration_annotation_sum, `Mast cells`  = 'Immune cells')
filtered_integration_annotation_sum <- RenameIdents(filtered_integration_annotation_sum, `Myeloid cells`  = 'Immune cells')
filtered_integration_annotation_sum <- RenameIdents(filtered_integration_annotation_sum, `B cells`  = 'Immune cells')
filtered_integration_annotation_sum <- RenameIdents(filtered_integration_annotation_sum, `Plasma cells`  = 'Immune cells')
filtered_integration_annotation_sum <- RenameIdents(filtered_integration_annotation_sum, `T/NK cells`  = 'Immune cells')

filtered_integration_annotation_sum <- RenameIdents(filtered_integration_annotation_sum, `Fibroblasts`  = 'Non-immune cells')
filtered_integration_annotation_sum <- RenameIdents(filtered_integration_annotation_sum, `Pericytes`  = 'Non-immune cells')
filtered_integration_annotation_sum <- RenameIdents(filtered_integration_annotation_sum, `Endothelial cells`  = 'Non-immune cells')
filtered_integration_annotation_sum <- RenameIdents(filtered_integration_annotation_sum, `Epithelial cells`  = 'Non-immune cells')

DimPlot(filtered_integration_annotation_sum, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Unsupervised clustering")
DimPlot(filtered_integration_annotation_sum, reduction = 'umap', raster=FALSE, label = FALSE, cols = c('Immune cells' = 'lightsalmon', 'Non-immune cells' = 'grey'))  + NoLegend() + NoAxes() 


################################################################################ Dotplot

# Define an order of cluster identities
Dotplot_levels <- c('Epithelial cells','Endothelial cells','Pericytes','Fibroblasts','Mast cells','Myeloid cells','Plasma cells','B cells','T/NK cells')

# Relevel object@ident
filtered_integration_annotation@active.ident <- factor(x = filtered_integration_annotation@active.ident, levels = Dotplot_levels)


cluster_features <- c('PTPRC', 'CD3E', 'CD8A', 'CD4', 'FOXP3', 'IL2RA', 'KLRD1', 'NKG7', 'KLRF1', 'MS4A1', 'CD19', 'MZB1', 'JCHAIN', 
                      'CD14', 'CFP', 'APOBEC3A', 'LAMP3', 'TREM2', 'KIT', 'IL1RL1', 'DCN', 'FN1','ACTA2', 'PDGFRB', 'RGS5', 'PLVAP','VWF','PECAM1','KRT18', 'CDH1', 'EPCAM')
DotPlot(filtered_integration_annotation, features = cluster_features)+ RotatedAxis() + xlab(NULL) + ylab(NULL) 







