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

######################## Change dimensionality and visualize UMAP ####################
for (i in 10:40) {
  immune_cell <- FindNeighbors(immune_cell, reduction = "pca", dims = 1:i)
  immune_cell <- FindClusters(immune_cell, resolution = 0.2)
  immune_cell <- RunUMAP(immune_cell, reduction = "pca", dims = 1:i) 
  
  # plot the results and save as PNG file
  plot_title <- paste("Unsupervised clustering 1:", i, sep = "")
  png_1_i <- DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) +
    ggtitle(plot_title) +
    NoLegend() +
    FeaturePlot(immune_cell, features = 'CD8A', min.cutoff = '0', max.cutoff = '3', raster=FALSE)

  ggsave(paste0("immune_png_1_", i, ".png"), png_1_i, width = 16, height = 7, units = "in")
}
##################################### Manual #########################################

immune_cell <- FindNeighbors(immune_cell, reduction = "pca", dims = 1:44)
immune_cell <- FindClusters(immune_cell, resolution = 5)
immune_cell <- RunUMAP(immune_cell, reduction = "pca", dims = 1:44)

DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() |
  FeaturePlot(immune_cell, features = 'TRAV1-2', min.cutoff = '0', max.cutoff = '5', raster=FALSE)

DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() |
  FeaturePlot(immune_cell, features = 'TRDV1', min.cutoff = '0', max.cutoff = '5', raster=FALSE)

DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() |
  FeaturePlot(immune_cell, features = 'TRDV2', min.cutoff = '0', max.cutoff = '5', raster=FALSE)

DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() |
  FeaturePlot(immune_cell, features = 'TRGC1', min.cutoff = '0', max.cutoff = '5', raster=FALSE)

DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() |
  FeaturePlot(immune_cell, features = 'TRGC2', min.cutoff = '0', max.cutoff = '5', raster=FALSE)

DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() |
  FeaturePlot(immune_cell, features = 'TRDC', min.cutoff = '0', max.cutoff = '5', raster=FALSE)


DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() |
  FeaturePlot(immune_cell, features = 'NKG7', min.cutoff = '0', max.cutoff = '5', raster=FALSE)

DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() |
  FeaturePlot(immune_cell, features = 'KLRD1', min.cutoff = '0', max.cutoff = '5', raster=FALSE)

DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() |
  FeaturePlot(immune_cell, features = 'KLRF1', min.cutoff = '0', max.cutoff = '5', raster=FALSE)


DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() |
  FeaturePlot(immune_cell, features = 'CD8A', min.cutoff = '0', max.cutoff = '5', raster=FALSE)

DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() |
  FeaturePlot(immune_cell, features = 'FOXP3', min.cutoff = '0', max.cutoff = '5', raster=FALSE)

DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() |
  FeaturePlot(immune_cell, features = 'CD14', min.cutoff = '0', max.cutoff = '5', raster=FALSE)


DefaultAssay(immune_cell) <- "ADT"

DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() |
  FeaturePlot(immune_cell, features = 'CD8', min.cutoff = '0', max.cutoff = '5', raster=FALSE)

DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() |
  FeaturePlot(immune_cell, features = 'CD4', min.cutoff = '0', max.cutoff = '5', raster=FALSE)

DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() |
  FeaturePlot(immune_cell, features = 'CD56', min.cutoff = '0', max.cutoff = '0.5', raster=FALSE)

DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() |
  FeaturePlot(immune_cell, features = 'TCR-V-alpha', min.cutoff = '0', max.cutoff = '3', raster=FALSE)

DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() |
  FeaturePlot(immune_cell, features = 'TCR-gamma-delta', min.cutoff = '1', max.cutoff = '2', raster=FALSE)


DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() |
  FeaturePlot(immune_cell, features = 'TCR-gamma-delta', raster=FALSE)

DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() |
  FeaturePlot(immune_cell, features = 'TCR-gamma-delta', raster=FALSE)


DefaultAssay(immune_cell) <- "integrated"

# result <- FindMarkers(immune_cell, ident.1 = "8", ident.2 = 1, min.pct = 0.5, logfc.threshold = log(2))


######################## Change resolution and visualize UMAP ####################
j_value <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9)

for (j in j_value) {
  filtered_integration <- FindNeighbors(filtered_integration, reduction = "pca", dims = 1:44)
  filtered_integration <- FindClusters(filtered_integration, resolution = j)
  filtered_integration <- RunUMAP(filtered_integration, reduction = "pca", dims = 1:44) 
  
  # plot the results and save as PNG file
  plot_title <- paste("Unsupervised clustering 1:10_", j, sep = "")
  png_1_j <- DimPlot(filtered_integration, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) +
    ggtitle(plot_title) +
    NoLegend() +
    FeaturePlot(filtered_integration, features = 'TRAV1-2', min.cutoff = '0', max.cutoff = '5', raster=FALSE)
  
  ggsave(paste0("png_1_44_", j, ".png"), png_1_j, width = 16, height = 7, units = "in")
}
##################################### Manual(final) #########################################

immune_cell <- FindNeighbors(immune_cell, reduction = "pca", dims = 1:44)
immune_cell <- FindClusters(immune_cell, resolution = 5)
immune_cell <- RunUMAP(immune_cell, reduction = "pca", dims = 1:44)

DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() |
  FeaturePlot(immune_cell, features = 'CD8A', min.cutoff = '0', max.cutoff = '5', raster=FALSE)

##################################### Drawing UMAP #########################################
colnames(immune_cell@meta.data)
nrow(immune_cell@meta.data)

DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() 

DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = FALSE)+ NoAxes() + ggtitle(NULL)+ NoLegend() 
DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE)+ NoAxes() + ggtitle(NULL) 
DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE)

# Exported with 2000/2000 size png

##################################### Save data: clustering #########################################

saveRDS(immune_cell, file="D:/DATA/Gastric cancer/immune_cell_UMAP.rds")

immune_cell_annotation <- immune_cell
immune_cell_annotation <- readRDS(file="D:/DATA/Gastric cancer/immune_cell_UMAP.rds")

##################################### Annotation : marker - non T/NK cells #########################################

immune_cell_annotation <- RenameIdents(immune_cell_annotation, `57`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `43`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `27`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `72`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `12`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `0`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `81`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `13`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `52`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `6`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `41`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `18`  = 'B cells')
DimPlot(immune_cell_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() 

immune_cell_annotation <- RenameIdents(immune_cell_annotation, `77`  = 'Proliferating B cells')
DimPlot(immune_cell_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() 


immune_cell_annotation <- RenameIdents(immune_cell_annotation, `49`  = 'Mast cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `28`  = 'Mast cells')
DimPlot(immune_cell_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() 


immune_cell_annotation <- RenameIdents(immune_cell_annotation, `85`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `58`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `63`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `39`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `80`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `42`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `47`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `56`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `36`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `78`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `64`  = 'Plasma cells')
DimPlot(immune_cell_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() 


immune_cell_annotation <- RenameIdents(immune_cell_annotation, `45`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `10`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `8`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `55`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `50`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `23`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `83`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `75`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `35`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `82`  = 'Myeloid cells')
DimPlot(immune_cell_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() 

##################################### Annotation : marker - T/NK cells #########################################


NK_T_cluster <- c("61", "30", "54", "53", "66", "25", "22", "17", "46", "79", "70", "15", "11", "5", "20", "59", "16",
                    "37", "86", "68", "34", "76", "3", "44", "38", "71", "84", "9", "40", "4", "1", "14", "29", "26", "74",
                    "32", "31", "33", "24", "67", "60", "65", "7", "48", "69", "51", "73", "21", "19", "2", "62")
  

NK_T_anno <- subset(x = immune_cell_annotation, idents = NK_T_cluster)


##################################### MNN Subclustering #########################################
T_NK_subset<-NK_T_anno
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

################################# MNN Subclustering 完##############################

ElbowPlot(NK_T_anno_re)

NK_T_anno_re <- NK_T_anno
NK_T_anno_re <- FindNeighbors(NK_T_anno_re, reduction = "pca", dims = 1:14)
NK_T_anno_re <- FindClusters(NK_T_anno_re, resolution = 5)
NK_T_anno_re <- RunUMAP(NK_T_anno_re, reduction = "pca", dims = 1:14)

DimPlot(NK_T_anno_re, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() |
  FeaturePlot(NK_T_anno_re, features = 'CD8A', min.cutoff = '0', max.cutoff = '5', raster=FALSE)

FeaturePlot(NK_T_anno_re, features = 'CD8A', min.cutoff = '0', max.cutoff = '5', raster=FALSE)|
  FeaturePlot(NK_T_anno_re, features = 'CD4', min.cutoff = '0', max.cutoff = '5', raster=FALSE)

FeaturePlot(NK_T_anno_re, features = 'adt_CD8', min.cutoff = '1', max.cutoff = '4', raster=FALSE)|
  FeaturePlot(NK_T_anno_re, features = 'adt_CD4', min.cutoff = '1', max.cutoff = '3', raster=FALSE)




VlnPlot(CD8_anno, features = "SVIL", pt.size = 0, split.by = "seurat_clusters")

















DimPlot(NK_T_anno, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() |
  DimPlot(NK_T_anno, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend()
FeaturePlot(NK_T_anno, features = 'CD8A', min.cutoff = '0', max.cutoff = '5', raster=FALSE)|
  FeaturePlot(NK_T_anno, features = 'CD4', min.cutoff = '0', max.cutoff = '5', raster=FALSE)


DimPlot(NK_T_anno, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() 
FeaturePlot(NK_T_anno, features = 'adt_CD8', min.cutoff = '1', max.cutoff = '4', raster=FALSE)|
  FeaturePlot(NK_T_anno, features = 'adt_CD4', min.cutoff = '1', max.cutoff = '3', raster=FALSE)


dittoBoxPlot(NK_T_anno_re, assay = "ADT", "TCR-gamma-delta", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
             boxplot.width = 0.7, min=1.2,max=1.6) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))

dittoBoxPlot(NK_T_anno_re, assay = "ADT", "CD56", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
             boxplot.width = 0.7, max=2) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))

dittoBoxPlot(NK_T_anno, assay = "ADT", "CD4", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
             boxplot.width = 0.7, max=1.5) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))


dittoBoxPlot(NK_T_anno, assay = "integrated", "TCR-gamma-delta", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
             boxplot.width = 0.7, min=1.2,max=1.6) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))

dittoBoxPlot(NK_T_anno, assay = "integrated", "CD56", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
             boxplot.width = 0.7, max=2) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))





NK_T_anno_i <- subset(x = NK_T_anno_re, idents = 30)
FeatureScatter(NK_T_anno_i, feature1 = "adt_CD4", feature2 = "adt_CD8")
FeatureScatter(NK_T_anno_i, feature1 = "adt_TCR-gamma-delta", feature2 = "adt_CD8")
FeatureScatter(NK_T_anno_i, feature1 = "adt_TCR-gamma-delta", feature2 = "adt_CD56")
FeatureScatter(NK_T_anno_i, feature1 = "adt_CD8", feature2 = "adt_CD56")
FeatureScatter(NK_T_anno_i, feature1 = "CD4", feature2 = "CD8A")

############################################## Annotation ############################



NK_T_cluster <- c(61, 30, 54, 53, 66, 25, 22, 17, 46, 79, 70, 15, 11, 5, 20, 59, 16, 37, 86, 68, 34, 76, 3, 44, 38, 71, 84, 9, 40, 4, 1, 14, 29, 26, 74, 32, 31, 33, 24, 67, 60, 65, 7, 48, 69, 51, 73, 21, 19, 2, 62)
NK_T_cluster <- c("61", "30", "54", "53", "66", "25", "22", "17", "46", "79", "70", "15", "11", "5", "20", "59", "16", "37", "86", "68", "34", "76", "3", "44", "38", "71", "84", "9", "40", "4", "1", "14", "29", "26", "74", "32", "31", "33", "24", "67", "60", "65", "7", "48", "69", "51", "73", "21", "19", "2", "62")

DimPlot(immune_anno_ADT, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() 
view(immune_cell_annotation@meta.data)
view(immune_anno_ADT@meta.data)


immune_anno_ADT <- subset(x = immune_cell_annotation, seurat_clusters == NK_T_cluster)
FeatureScatter(immune_anno_ADT, feature1 = "adt_CD56", feature2 = "adt_CD8")+ xlim(0, 3) + ylim(0, 4)+geom_point(color = "lightblue", size=1) + geom_density_2d()
FeatureScatter(immune_anno_ADT, feature1 = "adt_TCR-alpha-beta", feature2 = "adt_TCR-V-alpha")+ xlim(0, 2) + ylim(0, 2.5)+geom_point(color = "lightblue", size=1) + geom_density_2d()
FeatureScatter(immune_anno_ADT, feature1 = "adt_TCR-alpha-beta", feature2 = "adt_TCR-gamma-delta")+ xlim(0, 3) + ylim(0.5, 2)
FeatureScatter(immune_cell_annotation, feature1 = "adt_CD56", feature2 = "adt_CD8")+ xlim(0, 3) + ylim(0, 4)+geom_point(color = "lightblue", size=1) + geom_density_2d()


immune_cell_annotation <- RenameIdents(immune_cell_annotation, `61`  = 'MAIT cells')
DimPlot(immune_cell_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() 

immune_cell_annotation <- RenameIdents(immune_cell_annotation, `30`  = 'NK cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `54`  = 'NK cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `53`  = 'NK cells')
DimPlot(immune_cell_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() 

immune_cell_annotation <- RenameIdents(immune_cell_annotation, `66`  = 'Gamma-delta T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `25`  = 'Gamma-delta T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `22`  = 'Gamma-delta T cells')
DimPlot(immune_cell_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() 


immune_cell_annotation <- RenameIdents(immune_cell_annotation, `17`  = 'Regulatory CD4 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `46`  = 'Regulatory CD4 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `79`  = 'Regulatory CD4 T cells')
DimPlot(immune_cell_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() 


immune_cell_annotation <- RenameIdents(immune_cell_annotation, `70`  = 'Proliferating T cells')
DimPlot(immune_cell_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() 


immune_cell_annotation <- RenameIdents(immune_cell_annotation, `15`  = 'CD4 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `11`  = 'CD4 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `5`  = 'CD4 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `20`  = 'CD4 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `59`  = 'CD4 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `16`  = 'CD4 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `37`  = 'CD4 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `86`  = 'CD4 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `68`  = 'CD4 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `34`  = 'CD4 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `76`  = 'CD4 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `3`  = 'CD4 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `44`  = 'CD4 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `38`  = 'CD4 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `71`  = 'CD4 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `84`  = 'CD4 T cells')
DimPlot(immune_cell_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() 

immune_cell_annotation <- RenameIdents(immune_cell_annotation, `9`  = 'CD8 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `40`  = 'CD8 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `4`  = 'CD8 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `1`  = 'CD8 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `14`  = 'CD8 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `29`  = 'CD8 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `26`  = 'CD8 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `74`  = 'CD8 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `32`  = 'CD8 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `31`  = 'CD8 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `33`  = 'CD8 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `24`  = 'CD8 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `67`  = 'CD8 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `60`  = 'CD8 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `65`  = 'CD8 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `7`  = 'CD8 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `48`  = 'CD8 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `69`  = 'CD8 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `51`  = 'CD8 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `73`  = 'CD8 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `21`  = 'CD8 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `19`  = 'CD8 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `2`  = 'CD8 T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `62`  = 'CD8 T cells')
DimPlot(immune_cell_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() 


FeatureScatter(CD8_anno_0, feature1 = "adt_CD4", feature2 = "adt_CD8")

##################################### Annotation : marker Figure #########################################

# result <- FindMarkers(immune_cell, ident.1 = "4", ident.2 = "0", min.pct = 0.5, logfc.threshold = log(2))

DimPlot(immune_cell_annotation, reduction = 'umap', raster=FALSE, label = T) + NoLegend() 
DimPlot(immune_cell_annotation, reduction = 'umap', raster=FALSE, label = F) + NoLegend() + NoAxes()

# Markers
immune_markers <- c("CD3E", "CD3G", "CD8A", "CD4", "FOXP3", "IL2RA", "KLRD1", "KLRF1", "GNLY", "NKG7", "TRAV1-2", "TRDV1", "TRDV2", "TRGC1", "TRDC", "MKI67",
                    "MS4A1", "CD19", "MZB1", "JCHAIN", "KIT", "IL1RL1", "ITGB2", "CD14", "FCN1", "CD86", "CSF1R", "CD1C")


# Dot plots
DotPlot(immune_cell_annotation, features = immune_markers)+ RotatedAxis() + xlab(NULL) + ylab(NULL) 

# Create an empty list to store the plots
plot_list <- list()

for (marker in immune_markers) {
  g1 <- FeaturePlot(immune_cell_annotation, features = marker, min.cutoff = '1', max.cutoff = '4', raster=FALSE)+ NoAxes()
  
  # Add the plot to the list
  plot_list[[which(immune_markers == marker)]] <- g1
}

# Arrange the plots in a grid
g1_arrange <- grid.arrange(grobs = plot_list, nrow = 4)

ggsave(paste0("arrange nrow 4.png"), g1_arrange, width = 20, height = 10, units = "in", limitsize = FALSE)



##################################### Annotation : marker(ADT) #########################################
DefaultAssay(immune_cell_annotation) <-"ADT"

# Markers
immune_markers <- c("CD8", "CD4", "CD56", "CD19", "CD11c", "CD14")

# Dot plots
DotPlot(immune_cell_annotation, features = immune_markers)+ RotatedAxis() + xlab(NULL) + ylab(NULL) 

FeaturePlot(immune_cell_annotation, features = "CD8", min.cutoff = '1', max.cutoff = '3', raster=FALSE)+ NoAxes()

ADT_CD8 <- FeaturePlot(immune_cell_annotation, features = "CD8", min.cutoff = '1', max.cutoff = '3', raster=FALSE)+ NoAxes()
ADT_CD4 <-FeaturePlot(immune_cell_annotation, features = "CD4", min.cutoff = '1', max.cutoff = '2', raster=FALSE)+ NoAxes()
ADT_CD56 <-FeaturePlot(immune_cell_annotation, features = "CD56", min.cutoff = '0.5', max.cutoff = '1', raster=FALSE)+ NoAxes()
ADT_CD19 <-FeaturePlot(immune_cell_annotation, features = "CD19", min.cutoff = '1', max.cutoff = '3', raster=FALSE)+ NoAxes()
ADT_CD11c <-FeaturePlot(immune_cell_annotation, features = "CD11c", min.cutoff = '1', max.cutoff = '3', raster=FALSE)+ NoAxes()
ADT_CD14 <-FeaturePlot(immune_cell_annotation, features = "CD14", min.cutoff = '0.5', max.cutoff = '1.5', raster=FALSE)+ NoAxes()

plot_list2 <- list(ADT_CD8, ADT_CD4, ADT_CD56, ADT_CD19, ADT_CD11c, ADT_CD14)
g2_arrange <- grid.arrange(grobs = plot_list2, nrow = 1)

DefaultAssay(immune_cell_annotation) <-"integrated"


##################################### Cell counts #########################################

view(immune_cell_annotation@meta.data)

immune_cell_annotation@meta.data$clusters <- immune_cell_annotation@active.ident

SCpubr::do_BarPlot(sample = immune_cell_annotation, 
                         group.by = "clusters", 
                         legend.position = "none", 
                         plot.title = "",  xlab = "",  ylab = "Number of cells per cluster", font.size = 11)


##################################### Save data: annotation #########################################

saveRDS(immune_cell_annotation, file="D:/DATA/Gastric cancer/immune_cell_annotation.rds")

