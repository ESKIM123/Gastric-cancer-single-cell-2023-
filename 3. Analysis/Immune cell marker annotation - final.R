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
view(immune_cell@meta.data)
DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("Initial")

immune_cell <- RunPCA(immune_cell, verbose = FALSE)
immune_cell <- FindNeighbors(immune_cell, reduction = "pca", dims = 1:44)
immune_cell <- FindClusters(immune_cell, resolution = 20)
immune_cell <- RunUMAP(immune_cell, reduction = "pca", dims = 1:44)

DimPlot(immune_cell, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:44 20)")+ NoLegend() 


FeaturePlot(immune_cell, features = 'TRDC', min.cutoff = '0', max.cutoff = '5', raster=FALSE)
FeaturePlot(immune_cell, features = 'adt_CD8', min.cutoff = '1', max.cutoff = '4', raster=FALSE)
FeaturePlot(immune_cell, features = 'adt_CD4', min.cutoff = '1', max.cutoff = '3', raster=FALSE)


NK_T_anno_i <- subset(x = NK_T_anno_re, idents = 30)
FeatureScatter(NK_T_anno_i, feature1 = "adt_CD4", feature2 = "adt_CD8")
FeatureScatter(NK_T_anno_i, feature1 = "adt_TCR-gamma-delta", feature2 = "adt_CD8")
FeatureScatter(NK_T_anno_i, feature1 = "adt_TCR-gamma-delta", feature2 = "adt_CD56")
FeatureScatter(NK_T_anno_i, feature1 = "adt_CD8", feature2 = "adt_CD56")
FeatureScatter(NK_T_anno_i, feature1 = "CD4", feature2 = "CD8A")


VlnPlot(immune_cell, features = "CD8", pt.size = 0, split.by = "seurat_clusters")


##################################### Save data: clustering #########################################

saveRDS(immune_cell, file="D:/DATA/Gastric cancer/immune_cell_UMAP.rds")

immune_cell_annotation <- immune_cell
immune_cell_annotation <- readRDS(file="D:/DATA/Gastric cancer/immune_cell_UMAP.rds")


##################################### Annotation : marker  ######################################### OLD
immune_cell_annotation <- FindNeighbors(immune_cell_annotation, reduction = "pca", dims = 1:44)
immune_cell_annotation <- FindClusters(immune_cell_annotation, resolution = 5)
immune_cell_annotation <- RunUMAP(immune_cell_annotation, reduction = "pca", dims = 1:44)

DimPlot(immune_cell_annotation, reduction = 'umap', group.by = 'seurat_clusters', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() 

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







##################################### Annotation : marker - non T/NK cells ######################################### NEW (unfinished)

immune_cell_annotation <- RenameIdents(immune_cell_annotation, `196`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `116`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `152`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `133`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `19`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `38`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `99`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `184`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `197`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `204`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `106`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `237`  = 'Myeloid cells')

immune_cell_annotation <- RenameIdents(immune_cell_annotation, `172`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `101`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `222`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `170`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `46`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `17`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `107`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `51`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `135`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `102`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `189`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `80`  = 'Myeloid cells')

immune_cell_annotation <- RenameIdents(immune_cell_annotation, `92`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `226`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `223`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `163`  = 'Myeloid cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `225`  = 'Myeloid cells')
DimPlot(immune_cell_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("(1:44 20)")+ NoLegend() 


immune_cell_annotation <- RenameIdents(immune_cell_annotation, `138`  = 'Mast cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `180`  = 'Mast cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `55`  = 'Mast cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `60`  = 'Mast cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `232`  = 'Mast cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `112`  = 'Mast cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `241`  = 'Mast cells')
DimPlot(immune_cell_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("(1:44 20)")+ NoLegend() 

immune_cell_annotation <- RenameIdents(immune_cell_annotation, `187`  = 'Proliferating B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `228`  = 'Proliferating B cells')
DimPlot(immune_cell_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("(1:44 20)")+ NoLegend() 



immune_cell_annotation <- RenameIdents(immune_cell_annotation, `108`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `105`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `57`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `198`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `113`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `123`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `167`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `139`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `42`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `3`  = 'B cells')

immune_cell_annotation <- RenameIdents(immune_cell_annotation, `169`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `238`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `47`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `25`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `15`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `4`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `183`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `210`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `21`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `41`  = 'B cells')

immune_cell_annotation <- RenameIdents(immune_cell_annotation, `114`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `73`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `160`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `188`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `174`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `150`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `216`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `86`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `13`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `111`  = 'B cells')

immune_cell_annotation <- RenameIdents(immune_cell_annotation, `219`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `142`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `44`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `29`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `119`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `164`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `98`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `81`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `186`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `144`  = 'B cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `227`  = 'B cells')
DimPlot(immune_cell_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("(1:44 20)")+ NoLegend() 


immune_cell_annotation <- RenameIdents(immune_cell_annotation, `231`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `145`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `158`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `229`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `233`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `211`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `53`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `36`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `217`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `129`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `82`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `193`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `208`  = 'Plasma cells')

immune_cell_annotation <- RenameIdents(immune_cell_annotation, `136`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `56`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `77`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `182`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `240`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `155`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `221`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `94`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `200`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `122`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `161`  = 'Plasma cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `120`  = 'Plasma cells')
DimPlot(immune_cell_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("(1:44 20)")+ NoLegend() 


##################################### Annotation : marker - T/NK cells #########################################

immune_cell_annotation <- RenameIdents(immune_cell_annotation, ``  = 'Proliferating T cells')
immune_cell_annotation <- RenameIdents(immune_cell_annotation, `70`  = 'Proliferating T cells')
DimPlot(immune_cell_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() 


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

##################################### Annotation : marker Figure #########################################

# result <- FindMarkers(immune_cell, ident.1 = "4", ident.2 = "0", min.pct = 0.5, logfc.threshold = log(2))

DimPlot(immune_cell_annotation, reduction = 'umap', raster=FALSE, label = T) + NoLegend() 
DimPlot(immune_cell_annotation, reduction = 'umap', raster=FALSE, label = F) + NoLegend() + NoAxes()

# Markers
immune_markers <- c("CD3E", "CD3G", "CD8A", "CD8B", "CD4", "FOXP3", "IL2RA", "KLRD1", "KLRF1", "GNLY", "NKG7", "TRAV1-2", "TRDV1", "TRDV2", "TRGC1", "TRDC", "MKI67",
                    "MS4A1", "CD19", "MZB1", "JCHAIN", "KIT", "IL1RL1", "ITGB2", "CD14", "FCN1", "CD86", "CSF1R", "CD1C")


# Dot plots
DotPlot(immune_cell_annotation, features = immune_markers)+ RotatedAxis() + xlab(NULL) + ylab(NULL) + theme(axis.text.x = element_text(face = "italic"))
DotPlot(immune_cell_annotation, features = immune_markers)+ RotatedAxis() + xlab(NULL) + ylab(NULL) +NoLegend() + theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))


FeaturePlot(immune_cell_annotation, features = "CD3E", min.cutoff = '1', max.cutoff = '4', raster=FALSE)+ NoAxes()+ theme(plot.title = element_text(face = "italic"))

# Create an empty list to store the plots
plot_list <- list()

for (marker in immune_markers) {
  g1 <- FeaturePlot(immune_cell_annotation, features = marker, min.cutoff = '1', max.cutoff = '4', raster=FALSE)+ NoAxes()+ theme(plot.title = element_text(face = "italic"))
  
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

##################################### Save data: annotation #########################################
saveRDS(immune_cell_annotation, file="D:/DATA/Gastric cancer/immune_cell_annotation_final_OLD.rds")
immune_cell_annotation <- readRDS(file="D:/DATA/Gastric cancer/immune_cell_annotation_final_OLD.rds")
##################################### Cell counts #########################################

view(immune_cell_annotation@meta.data)
immune_cell_annotation@meta.data$clusters <- immune_cell_annotation@active.ident

DimPlot(immune_cell_annotation, reduction = 'umap', group.by = 'clusters', raster=FALSE, label = TRUE) + ggtitle("(1:44 5)")+ NoLegend() 


immune_cell_annotation

dittoBarPlot(immune_cell_annotation, "clusters", group.by = "Patient_ID_SPC", scale = "percent")
dittoBarPlot(immune_cell_annotation, "clusters", group.by = "clusters", scale = "count")

  
dittoBarPlot(immune_cell_annotation, "clusters", group.by = "clusters", scale = "count",  retain.factor.levels = F, main = "")+ 
  scale_y_continuous(expand = c(0, 0)) + NoLegend()+ xlab(NULL)+ ylab(NULL)


do_BarPlot(sample = immune_cell_annotation, 
        group.by = "clusters", 
        legend.position = "none", 
        plot.title = "",  xlab = "",  ylab = "Number of cells per cluster", font.size = 11)



