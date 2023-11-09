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
immune_cell <- readRDS(file="D:/DATA/Gastric cancer/immune_cell_clinical_data_add.rds")
DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("Initial")

immune_cell <- RunPCA(immune_cell, verbose = FALSE)
immune_cell <- FindNeighbors(immune_cell, reduction = "pca", dims = 1:44)
immune_cell <- FindClusters(immune_cell, resolution = 5)
immune_cell <- RunUMAP(immune_cell, reduction = "pca", dims = 1:44)


DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() |
  FeaturePlot(immune_cell, features = 'TRDC', min.cutoff = '0', max.cutoff = '3', raster=FALSE)




RunPCA(reduction.name = 'apca')

DefaultAssay(immune_cell_subset)<-"integrated"








##################################### MNN Subclustering #########################################
immune_cell_subset<-immune_cell

DefaultAssay(T_NK_subset) <- 'RNA'
T_NK_subset <- NormalizeData(T_NK_subset) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

DefaultAssay(immune_cell_subset) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
VariableFeatures(immune_cell_subset) <- rownames(immune_cell_subset[["ADT"]])
immune_cell_subset <- NormalizeData(immune_cell_subset, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')
DefaultAssay(immune_cell_subset)<-"integrated"

# Identify multimodal neighbors. These will be stored in the neighbors slot, 
# and can be accessed using T_NK_subset[['weighted.nn']]
# The WNN graph can be accessed at T_NK_subset[["wknn"]], 
# and the SNN graph used for clustering at T_NK_subset[["wsnn"]]
# Cell-specific modality weights can be accessed at T_NK_subset$RNA.weight
immune_cell_subset <- FindMultiModalNeighbors(
  immune_cell_subset, reduction.list = list("pca", "apca"), 
  dims.list = list(1:50, 1:43), modality.weight.name = "RNA.weight"
)

immune_cell_subset <- RunUMAP(immune_cell_subset, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
immune_cell_subset <- FindClusters(immune_cell_subset, graph.name = "wsnn", algorithm = 3, resolution = 5, verbose = FALSE)

DimPlot(immune_cell_subset, reduction = 'wnn.umap', raster=F, label = TRUE, repel = TRUE) + NoLegend()


FeaturePlot(immune_cell_subset, features = 'CD8A', min.cutoff = '0', max.cutoff = '5', raster=FALSE, reduction = 'wnn.umap')|
  FeaturePlot(immune_cell_subset, features = 'CD4', min.cutoff = '0', max.cutoff = '5', raster=FALSE, reduction = 'wnn.umap')

FeaturePlot(immune_cell_subset, features = 'adt_CD8', min.cutoff = '0', max.cutoff = '5', raster=FALSE, reduction = 'wnn.umap')
FeaturePlot(immune_cell_subset, features = 'adt_CD4', min.cutoff = '0', max.cutoff = '5', raster=FALSE, reduction = 'wnn.umap')

FeaturePlot(immune_cell_subset, features = 'adt_CD56', min.cutoff = '0', max.cutoff = '3', raster=FALSE, reduction = 'wnn.umap')
FeaturePlot(immune_cell_subset, features = 'adt_TCR-V-alpha', min.cutoff = '0', max.cutoff = '3', raster=FALSE, reduction = 'wnn.umap')

  FeaturePlot(immune_cell_subset, features = 'adt_TCR-gamma-delta', min.cutoff = '0', max.cutoff = '5', raster=FALSE, reduction = 'wnn.umap')

dittoBoxPlot(immune_cell_subset, assay = "integrated", "CD4", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
             boxplot.width = 0.7, max=1.5) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))



VlnPlot(immune_cell_subset, features = "integrated.weight", group.by = 'seurat_clusters', sort = TRUE, pt.size = 0, raster=F) + NoLegend()


NK_T_anno_i <- subset(x = immune_cell_subset, idents = 2)
FeatureScatter(NK_T_anno_i, feature1 = "adt_CD4", feature2 = "adt_CD8")
FeatureScatter(NK_T_anno_i, feature1 = "adt_TCR-gamma-delta", feature2 = "adt_CD8")
FeatureScatter(NK_T_anno_i, feature1 = "adt_TCR-gamma-delta", feature2 = "adt_CD56")
FeatureScatter(NK_T_anno_i, feature1 = "adt_CD8", feature2 = "adt_CD56")
FeatureScatter(NK_T_anno_i, feature1 = "CD4", feature2 = "CD8A")

































immune_cell_subset <- RunUMAP(immune_cell_subset, reduction = 'pca', dims = 1:50, assay = 'RNA', 
                       reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
immune_cell_subset <- RunUMAP(immune_cell_subset, reduction = 'apca', dims = 1:43, assay = 'ADT', 
                       reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')


p3 <- DimPlot(immune_cell_subset, reduction = 'rna.umap', group.by = 'seurat_clusters', label = TRUE, 
              repel = TRUE, label.size = 2.5) + NoLegend()
p4 <- DimPlot(immune_cell_subset, reduction = 'adt.umap', group.by = 'seurat_clusters', label = TRUE, 
              repel = TRUE, label.size = 2.5) + NoLegend()
p3 + p4


p5 <- FeaturePlot(T_NK_subset, features = c("adt_CD4","adt_CD8","adt_CD56", "adt_TCR-V-alpha", "TCR-gamma-delta"),
                  reduction = 'wnn.umap', max.cutoff = 2, 
                  cols = c("lightgrey","darkgreen"), ncol = 3)
p6 <- FeaturePlot(T_NK_subset, features = c("rna_TRDC","rna_MPO","rna_AVP"), 
                  reduction = 'wnn.umap', max.cutoff = 3, ncol = 3)
p5 / p6

VlnPlot(T_NK_subset, features = "RNA.weight", group.by = 'celltype.l2', sort = TRUE, pt.size = 0.1) +
  NoLegend()


FeaturePlot(T_NK_subset, features = "TCR-gamma-delta", reduction = 'wnn.umap', max.cutoff = 10)

################################# MNN Subclustering 完##############################