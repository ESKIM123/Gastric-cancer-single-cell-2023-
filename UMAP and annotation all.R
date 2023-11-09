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


# Get data ----
filtered_integration <- readRDS(file="D:/DATA/Gastric cancer/integrated_seurat_filtered_SPC_ID_singlet.rds")
DimPlot(filtered_integration, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("Unsupervised clustering 1:30 / 2 (base)")


######################## Change dimensionality and visualize UMAP ####################
for (i in 10:20) {
  filtered_integration <- FindNeighbors(filtered_integration, reduction = "pca", dims = 1:i)
  filtered_integration <- FindClusters(filtered_integration, resolution = 5)
  filtered_integration <- RunUMAP(filtered_integration, reduction = "pca", dims = 1:i) 
  
  # plot the results and save as PNG file
  plot_title <- paste("Unsupervised clustering 1:", i, sep = "")
  png_1_i <- DimPlot(filtered_integration, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) +
    ggtitle(plot_title) +
    NoLegend() +
    FeaturePlot(filtered_integration, features = 'HSPB1', min.cutoff = '0', max.cutoff = '3', raster=FALSE)
  
  ggsave(paste0("png_1_", i, ".png"), png_1_i, width = 16, height = 7, units = "in")
}
##################################### Manual #########################################

filtered_integration <- FindNeighbors(filtered_integration, reduction = "pca", dims = 1:10)
filtered_integration <- FindClusters(filtered_integration, resolution = 3)
filtered_integration <- RunUMAP(filtered_integration, reduction = "pca", dims = 1:10)

DimPlot(filtered_integration, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:10 3)")+ NoLegend() 


######################## Change resolution and visualize UMAP ####################
j_value <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2)

for (j in j_value) {
  filtered_integration <- FindNeighbors(filtered_integration, reduction = "pca", dims = 1:10)
  filtered_integration <- FindClusters(filtered_integration, resolution = j)
  filtered_integration <- RunUMAP(filtered_integration, reduction = "pca", dims = 1:10) 
  
  # plot the results and save as PNG file
  plot_title <- paste("Unsupervised clustering 1:10_", j, sep = "")
  png_1_j <- DimPlot(filtered_integration, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) +
    ggtitle(plot_title) +
    NoLegend() +
    FeaturePlot(filtered_integration, features = 'HSPB1', min.cutoff = '0', max.cutoff = '3', raster=FALSE)
  
  ggsave(paste0("png_1_10_", j, ".png"), png_1_j, width = 16, height = 7, units = "in")
}
##################################### Manual(final) #########################################

filtered_integration <- FindNeighbors(filtered_integration, reduction = "pca", dims = 1:10)
filtered_integration <- FindClusters(filtered_integration, resolution = 0.2)
filtered_integration <- RunUMAP(filtered_integration, reduction = "pca", dims = 1:10)

DimPlot(filtered_integration, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:10 0.2)")+ NoLegend() |
FeaturePlot(filtered_integration, features = 'HSPB1', min.cutoff = '0', max.cutoff = '3', raster=FALSE)

##################################### Drawing UMAP #########################################
colnames(filtered_integration@meta.data)
nrow(filtered_integration@meta.data)

DimPlot(filtered_integration, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:10 0.2)")+ NoLegend() 

DimPlot(filtered_integration, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = FALSE)+ NoAxes() + ggtitle(NULL)+ NoLegend() 
DimPlot(filtered_integration, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE)+ NoAxes() + ggtitle(NULL) 
DimPlot(filtered_integration, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE)

DimPlot(filtered_integration, reduction = 'umap', group.by = 'Patient_ID_SPC', raster=FALSE, label = FALSE)+ NoAxes() + ggtitle(NULL)+ NoLegend() 
DimPlot(filtered_integration, reduction = 'umap', group.by = 'Patient_ID_SPC', raster=FALSE, label = F)+ NoAxes() + ggtitle(NULL) 

DimPlot(filtered_integration, reduction = 'umap', group.by = 'Type', raster=FALSE, label = FALSE)+ NoAxes() + ggtitle(NULL)+ NoLegend() 
DimPlot(filtered_integration, reduction = 'umap', group.by = 'Type', raster=FALSE, label = F)+ NoAxes() + ggtitle(NULL) 

# Exported with 2000/2000 size png

##################################### Save data: clustering #########################################

saveRDS(filtered_integration, file="D:/DATA/Gastric cancer/integrated_seurat_filltered_SPC_UMAP.rds")


##################################### Annotation #########################################

# result <- FindMarkers(filtered_integration, ident.1 = "6", ident.2 = c("5","1"), min.pct = 0.5, logfc.threshold = log(2))


cluster_features <- c('PTPRC', 'CD3E', 'CD8A', 'CD4', 'FOXP3', 'IL2RA', 'KLRD1', 'NKG7', 'KLRF1', 'MS4A1', 'CD19', 'MZB1', 'JCHAIN', 
                      'CD14', 'CFP', 'APOBEC3A', 'LAMP3', 'TREM2', 'KIT', 'IL1RL1', 'DCN', 'FN1','ACTA2', 'PDGFRB', 'RGS5', 'PLVAP',
                      'VWF','PECAM1','KRT18', 'CDH1', 'EPCAM', 'HSPB1', 'HSPA6', 'HSPA1A')
DotPlot(filtered_integration, features = cluster_features)+ RotatedAxis() + xlab(NULL) + ylab(NULL) 

filtered_integration_annotation <- filtered_integration
Idents(filtered_integration_annotation)

filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `0`  = '0: CD8 T cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `1`  = '1: CD4 T cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `5`  = '5: Treg cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `6`  = '6: Dying/Stressed T cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `4`  = '4: Plasma cells')

filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `2`  = '2: B cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `3`  = '3: Myeloid cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `8`  = '8: Mast cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `7`  = '7: Fibroblasts')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `10`  = '10: Endothelial cells')
filtered_integration_annotation <- RenameIdents(filtered_integration_annotation, `9`  = '9: Epithelial cells')
Idents(filtered_integration_annotation)

Dotplot_levels <- c('Epithelial cells','Endothelial cells','Pericytes','Fibroblasts','Mast cells','Myeloid cells','Plasma cells','B cells','T/NK cells')
DotPlot(filtered_integration_annotation, features = cluster_features)+ RotatedAxis() + xlab(NULL) + ylab(NULL) 


##################################### Save data: annotation #########################################

saveRDS(filtered_integration_annotation, file="D:/DATA/Gastric cancer/Annotation_all.rds")

##################################### Save data: immune subset #########################################
immune_cell <- subset(x = filtered_integration, idents = c("0", "1", "5", "2", "3", "4", "8"))
DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE)+ NoLegend() 

saveRDS(immune_cell, file="D:/DATA/Gastric cancer/immune_cell.rds")


