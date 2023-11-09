# load libraries
library(dplyr)
library(ggplot2)
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
library(SeuratDisk)


# Get data ----
immune_annotation <- readRDS(file="D:/DATA/Gastric cancer/immune_cell_annotation_clinical_data_add.rds")
DimPlot(immune_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Initial")

##################################### Subclustering #########################################
T_NK_subset <- subset(x = immune_annotation, ident = c("CD8 T cells", "CD4 T cells", "Proliferating T cells", "Regulatory CD4 T cells", "Gamma-delta T cells",
                                                          "NK cells", "MAIT cells"))

DimPlot(T_NK_subset, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Initial")

##################################### Subclustering #########################################
DefaultAssay(T_NK_subset) <- 'RNA'
T_NK_subset <- NormalizeData(T_NK_subset) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

DefaultAssay(T_NK_subset) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
VariableFeatures(T_NK_subset) <- rownames(T_NK_subset[["ADT"]])
T_NK_subset <- NormalizeData(T_NK_subset, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')


# Identify multimodal neighbors. These will be stored in the neighbors slot, 
# and can be accessed using T_NK_subset[['weighted.nn']]
# The WNN graph can be accessed at T_NK_subset[["wknn"]], 
# and the SNN graph used for clustering at T_NK_subset[["wsnn"]]
# Cell-specific modality weights can be accessed at T_NK_subset$RNA.weight
T_NK_subset <- FindMultiModalNeighbors(
  T_NK_subset, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
)

T_NK_subset <- RunUMAP(T_NK_subset, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
T_NK_subset <- FindClusters(T_NK_subset, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)

p1 <- DimPlot(T_NK_subset, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p2 <- DimPlot(T_NK_subset, reduction = 'wnn.umap', group.by = 'seurat_clusters', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p1 + p2

T_NK_subset <- RunUMAP(T_NK_subset, reduction = 'pca', dims = 1:30, assay = 'RNA', 
              reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
T_NK_subset <- RunUMAP(T_NK_subset, reduction = 'apca', dims = 1:18, assay = 'ADT', 
              reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')


p3 <- DimPlot(T_NK_subset, reduction = 'rna.umap', group.by = 'seurat_clusters', label = TRUE, 
              repel = TRUE, label.size = 2.5) + NoLegend()
p4 <- DimPlot(T_NK_subset, reduction = 'adt.umap', group.by = 'seurat_clusters', label = TRUE, 
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






