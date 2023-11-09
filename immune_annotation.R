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


# Get data ----
immune_clusters <- readRDS(file="D:/DATA/Gastric cancer/immune_clusters_SPC_ID.rds")

# DefaultAssay(immune_clusters) <- "integrated"
immune_clusters <- RunUMAP(immune_clusters, reduction = "pca", dims = 1:50) 
immune_clusters <- FindNeighbors(immune_clusters, reduction = "pca", dims = 1:50) 
immune_clusters <- FindClusters(immune_clusters, resolution = 5)

DimPlot(immune_clusters, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes()
DimPlot(immune_clusters, reduction = 'umap', raster=FALSE, label = FALSE) + NoLegend() +NoAxes()


immune_clusters_annotation <- immune_clusters


head(FindMarkers(immune_clusters, ident.1 = "71", ident.2 = "4", min.pct = 0.5, logfc.threshold = log(2)), n=20)

DimPlot(immune_clusters, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
FeaturePlot(immune_clusters, features = "TRAV1-2", min.cutoff =1, max.cutoff =3,raster=FALSE)

DimPlot(immune_clusters, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_clusters, features = "KLRB1", min.cutoff =1, max.cutoff =3,raster=FALSE)

DimPlot(immune_clusters, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_clusters, features = "TRDV1", min.cutoff =1, max.cutoff =3, raster=FALSE)

DimPlot(immune_clusters, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_clusters, features = "TRDV2", min.cutoff =1, max.cutoff =3, raster=FALSE)

DimPlot(immune_clusters, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_clusters, features = "TRGC1", min.cutoff =1, max.cutoff =3, raster=FALSE)

DimPlot(immune_clusters, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_clusters, features = "TRGC2", min.cutoff =1, max.cutoff =3, raster=FALSE)

DimPlot(immune_clusters, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_clusters, features = "TRDC", min.cutoff =1, max.cutoff =3, raster=FALSE)



################################################################################# Annotation


# B cells
#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "CD19", legend.position="none", min.cutoff =1, max.cutoff =3, pt.size = 2)
#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "MS4A1", legend.position="none", min.cutoff =1, max.cutoff =4, pt.size = 2)

immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '66' = "B cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '54' = "B cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '38' = "B cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `17` = "B cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `75` = "B cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `82` = "B cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `25` = "B cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `0` = "B cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `12` = "B cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `6` = "B cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `1` = "B cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `77` = "B cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `73` = "Proliferating B cells")

#DimPlot(immune_clusters_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes()


# Plasma cell
#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "SDC1", legend.position="none", min.cutoff =1, max.cutoff =3, pt.size = 2)
#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "JCHAIN", legend.position="none", min.cutoff =1, max.cutoff =4, pt.size = 2)

immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '70' = "Plasma cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '15' = "Plasma cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '36' = "Plasma cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `48` = "Plasma cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `83` = "Plasma cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `65` = "Plasma cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `68` = "Plasma cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `62` = "Plasma cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `79` = "Plasma cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `57` = "Plasma cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `78` = "Plasma cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `24` = "Plasma cells")

#DimPlot(immune_clusters_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes()


# Mast cell
#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "KIT", legend.position="none", min.cutoff =1, max.cutoff =3, pt.size = 2)
#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "IL1RL1", legend.position="none", min.cutoff =1, max.cutoff =4, pt.size = 2)

immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '32' = "Mast cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '50' = "Mast cells")

#DimPlot(immune_clusters_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes()


# Myeloid cell
#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "ITGB2", legend.position="none", min.cutoff =0, max.cutoff =4, pt.size = 2)
#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "CD14", legend.position="none", min.cutoff =1, max.cutoff =4, pt.size = 2)
#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "FCN1", legend.position="none", min.cutoff =1, max.cutoff =4, pt.size = 2)
#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "CSF1R", legend.position="none", min.cutoff =1, max.cutoff =4, pt.size = 2)
#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "CD86", legend.position="none", min.cutoff =1, max.cutoff =10, pt.size = 2)
#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "CD1C", legend.position="none", min.cutoff =1, max.cutoff =15, pt.size = 2)

immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '9' = "Monocytes")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '26' = "Monocytes")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '37' = "Monocytes")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '22' = "Monocytes")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '55' = "Monocytes")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '81' = "Monocytes")

immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '45' = "Dendritic cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '43' = "Dendritic cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '80' = "Dendritic cells")

#DimPlot(immune_clusters_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes()


# T/NK cell

#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "CD3G", legend.position="none", min.cutoff =0, max.cutoff =3, pt.size = 2)
#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "CD3E", legend.position="none", min.cutoff =0, max.cutoff =3, pt.size = 2)

#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "GNLY", legend.position="none", min.cutoff =0, max.cutoff =15, pt.size = 2)
#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "NCAM1", legend.position="none", min.cutoff =6, max.cutoff =8, pt.size = 2)

immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '31' = "NK cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '42' = "NK cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '23' = "gamma delta T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '29' = "gamma delta T cells")

#DimPlot(immune_clusters_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes()


#SCpubr::do_FeaturePlot(sample = immune_clusters_annotation, 
#                       idents.highlight = levels(immune_clusters_annotation)[!(levels(immune_clusters_annotation) 
#                                                                               %in% c('Dendritic cells','Monocytes','Plasma cells','B cells','Mast cells','Proliferating B cells'))], 
#                       features = c("MKI67"),min.cutoff =1, max.cutoff =10, pt.size = 2, legend.position="none")

#SCpubr::do_FeaturePlot(sample = immune_clusters_annotation, 
#                       idents.highlight = levels(immune_clusters_annotation)[!(levels(immune_clusters_annotation) 
#                                                                               %in% c('Dendritic cells','Monocytes','Plasma cells','B cells','Mast cells','Proliferating B cells'))], 
#                       features = c("HSPB1"),min.cutoff =1, max.cutoff =4, pt.size = 2, legend.position="none")

immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '64' = "Proliferating T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '59' = "Stressed / dying T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '41' = "Stressed / dying T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '71' = "Stressed / dying T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '76' = "Stressed / dying T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '40' = "Stressed / dying T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '30' = "Stressed / dying T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '69' = "Stressed / dying T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '52' = "Stressed / dying T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '53' = "Stressed / dying T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '58' = "Stressed / dying T cells")

#DimPlot(immune_clusters_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes()


# T cells specific

#SCpubr::do_FeaturePlot(sample = immune_clusters_annotation, 
#                       idents.highlight = levels(immune_clusters_annotation)[!(levels(immune_clusters_annotation) 
#                                                                               %in% c('Dendritic cells','Monocytes','Plasma cells','B cells','Mast cells','Proliferating B cells'))], 
#                       features = c("CD4"),min.cutoff =1, max.cutoff =3, pt.size = 2, legend.position="none")

#SCpubr::do_FeaturePlot(sample = immune_clusters_annotation, 
#                       idents.highlight = levels(immune_clusters_annotation)[!(levels(immune_clusters_annotation) 
#                                                                               %in% c('Dendritic cells','Monocytes','Plasma cells','B cells','Mast cells','Proliferating B cells'))], 
#                       features = c("CD8A"),min.cutoff =1, max.cutoff =4, pt.size = 2, legend.position="none")

#SCpubr::do_FeaturePlot(sample = immune_clusters_annotation, 
#                       idents.highlight = levels(immune_clusters_annotation)[!(levels(immune_clusters_annotation) 
#                                                                               %in% c('Dendritic cells','Monocytes','Plasma cells','B cells','Mast cells','Proliferating B cells'))], 
#                       features = c("FOXP3"),min.cutoff =1, max.cutoff =4, pt.size = 2, legend.position="none")


immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '4' = "Regulatory CD4 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '33' = "Regulatory CD4 T cells")
#DimPlot(immune_clusters_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes()

immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '49' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '47' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '72' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '46' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '3' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '7' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '8' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '61' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '16' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '21' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '19' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '14' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '28' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '51' = "MAIT cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '67' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '10' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '63' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '13' = "CD8 T cells")

#DimPlot(immune_clusters_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes()


immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '2' = "CD4 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '56' = "CD4 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '39' = "CD4 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '60' = "CD4 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '27' = "CD4 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '20' = "CD4 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '35' = "CD4 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '5' = "CD4 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '18' = "CD4 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '74' = "CD4 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '44' = "CD4 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '84' = "CD4 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '34' = "CD4 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '11' = "CD4 T cells")

DimPlot(immune_clusters_annotation, reduction = 'umap', raster=FALSE, label = TRUE, label.size = 6) + NoLegend()+NoAxes() +RotatedAxis()

immune_clusters_annotation[["immune_cluster"]]  <- immune_clusters_annotation@active.ident

saveRDS(immune_clusters_annotation, file="D:/DATA/Gastric cancer/immune_clusters_annotation.rds")


################################################################################# Remove Stressed cells ----

# Get data 
immune_clusters_annotation <- readRDS(file="D:/DATA/Gastric cancer/immune_clusters_annotation.rds")

# Remove Stressed cells & no patient information 

immune_no_stress <- subset(x = immune_clusters_annotation, subset = Patient_ID_SPC != 'FALSE')

immune_no_stress <- subset(x = immune_clusters_annotation, idents = c("CD4 T cells", "CD8 T cells", "Regulatory CD4 T cells", "B cells", 
                                                          "Plasma cells", "Mast cells", "Monocytes", "Dendritic cells",
                                                          "MAIT cells", "Proliferating T cells", "gamma delta T cells", "NK cells", "Proliferating B cells"))

DimPlot(immune_no_stress, reduction = 'umap', raster=FALSE, label = TRUE, label.size = 4) + NoLegend()+NoAxes()

# Re-cluster by SCT

#immune_no_stress_SCT <- SCTransform(immune_no_stress, method = "glmGamPoi", verbose = FALSE)
#immune_no_stress_SCT <- RunPCA(immune_no_stress_SCT, verbose = FALSE)
#immune_no_stress_SCT <- RunUMAP(immune_no_stress_SCT, dims = 1:30, verbose = FALSE)
#immune_no_stress_SCT <- FindNeighbors(immune_no_stress_SCT, dims = 1:30, verbose = FALSE)
#immune_no_stress_SCT <- FindClusters(immune_no_stress_SCT, verbose = FALSE)
#DimPlot(immune_no_stress_SCT, reduction = 'umap', raster=FALSE, label = TRUE, label.size = 4) + NoLegend()+NoAxes()

#saveRDS(immune_no_stress_SCT, file="D:/DATA/Gastric cancer/immune_no_stress_SCT.rds")

# Recluster by w/o SCT

immune_no_stress_PCA <- RunPCA(immune_no_stress, features = VariableFeatures(object = immune_no_stress))
immune_no_stress_PCA <- FindNeighbors(immune_no_stress_PCA, reduction = "pca", dims = 1:50) 
immune_no_stress_PCA <- FindClusters(immune_no_stress_PCA, resolution = 5)
immune_no_stress_PCA <- RunUMAP(immune_no_stress_PCA, reduction = "pca", dims = 1:50) 
DimPlot(immune_no_stress_PCA, reduction = 'umap', raster=FALSE, label = TRUE, label.size = 4) + NoLegend()+NoAxes()

saveRDS(immune_no_stress_PCA, file="D:/DATA/Gastric cancer/immune_no_stress_PCA.rds")


################################################################################# Annotation ----

immune_no_stress_PCA <- subset(x = immune_no_stress_PCA, idents != c(81, 56, 53, 64))
DimPlot(immune_no_stress_PCA, reduction = 'umap', raster=FALSE, label = TRUE, label.size = 4) + NoLegend()+NoAxes()

view(immune_no_stress_PCA@meta.data)

immune_no_stress_PCA_1 <- subset(x = immune_no_stress_PCA, subset = seurat_clusters == c("81", "56", "53", "64"))


immune_no_stress_PCA_1 <- subset(x = immune_no_stress_PCA, idents = 
                                   c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", 
                                     "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", 
                                     "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", 
                                     "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", 
                                     "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", 
                                     "50", "51", "52", "54", "55", "57", "58", "59", 
                                     "60", "61", "62", "63", "65", "66", "67", "68", "69", 
                                     "70", "71", "72", "73", "74", "75", "76", "77", "78", "79", 
                                     "80" ))


DimPlot(immune_no_stress_PCA_1, reduction = 'umap', raster=FALSE, label = TRUE, label.size = 4) + NoLegend()+NoAxes()


immune_no_stress_PCA_1 <- FindNeighbors(immune_no_stress_PCA_1, reduction = "pca", dims = 1:50) 
immune_no_stress_PCA_1 <- FindClusters(immune_no_stress_PCA_1, resolution = 5)
immune_no_stress_PCA_1 <- RunUMAP(immune_no_stress_PCA_1, reduction = "pca", dims = 1:50) 
DimPlot(immune_no_stress_PCA_1, reduction = 'umap', raster=FALSE, label = TRUE, label.size = 4) + NoLegend()+NoAxes()


##################################################################################
DefaultAssay(immune_T_cell) <- "integrated"

DefaultAssay(immune_no_stress_PCA_1) <- "ADT"

DimPlot(immune_no_stress_PCA_1, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_no_stress_PCA_1, features = "CD8", min.cutoff =1, max.cutoff =3,raster=FALSE)

DimPlot(immune_no_stress_PCA_1, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_no_stress_PCA_1, features = "TCR-gamma-delta", min.cutoff =1, max.cutoff =3,raster=FALSE)

DimPlot(immune_no_stress_PCA_1, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_no_stress_PCA_1, features = "TCR-V-alpha", min.cutoff =1, max.cutoff =3,raster=FALSE)




immune_T_cell <- subset(x = immune_no_stress_PCA_1, idents = c("11", "6", "1", "0", "5", "9", "8"))
immune_T_cell <- FindNeighbors(immune_T_cell, reduction = "pca", dims = 1:50) 
immune_T_cell <- FindClusters(immune_T_cell, resolution = 5)
immune_T_cell <- RunUMAP(immune_T_cell, reduction = "pca", dims = 1:50) 
DimPlot(immune_T_cell, reduction = 'umap', raster=FALSE, label = TRUE, label.size = 4) + NoLegend()+NoAxes()

DefaultAssay(immune_T_cell) <- "ADT"

DimPlot(immune_T_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_T_cell, features = "CD8", min.cutoff =1, max.cutoff =3,raster=FALSE)

DimPlot(immune_T_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_T_cell, features = "CD4", min.cutoff =1, max.cutoff =3,raster=FALSE)

DimPlot(immune_T_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_T_cell, features = "CD56", min.cutoff =1, max.cutoff =3,raster=FALSE)



result <- FindMarkers(immune_no_stress_PCA_1, ident.1 = "32", ident.2 = "8", min.pct = 0.5, logfc.threshold = log(2))



DimPlot(immune_T_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_T_cell, features = "TRAV1-2", min.cutoff =1, max.cutoff =3,raster=FALSE)

immune_T_cell <- RenameIdents(immune_T_cell, '55' = "MAIT cells")


DimPlot(immune_T_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_T_cell, features = "TRDC", min.cutoff =1, max.cutoff =3, raster=FALSE)

DimPlot(immune_T_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_T_cell, features = "TRDV1", min.cutoff =1, max.cutoff =3, raster=FALSE)

DimPlot(immune_T_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_T_cell, features = "TRDV2", min.cutoff =1, max.cutoff =3, raster=FALSE)

DimPlot(immune_T_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_T_cell, features = "TRGC1", min.cutoff =1, max.cutoff =3, raster=FALSE)

DimPlot(immune_T_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_T_cell, features = "TRGC2", min.cutoff =1, max.cutoff =3, raster=FALSE)



immune_T_cell <- RenameIdents(immune_T_cell, '11' = "gamma delta T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '66' = "gamma delta T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '52' = "gamma delta T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '12' = "gamma delta T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '10' = "gamma delta T cells")


DimPlot(immune_T_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_T_cell, features = "NCAM1", min.cutoff =1, max.cutoff =3, raster=FALSE)

DimPlot(immune_T_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_T_cell, features = "NKG7", raster=FALSE)

DimPlot(immune_T_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_T_cell, features = "KLRD1", raster=FALSE)

DimPlot(immune_T_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_T_cell, features = "GNLY", raster=FALSE)


immune_T_cell <- RenameIdents(immune_T_cell, '58' = "NK cells")
immune_T_cell <- RenameIdents(immune_T_cell, '21' = "NK cells")
immune_T_cell <- RenameIdents(immune_T_cell, '20' = "NK cells")

DimPlot(immune_T_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_T_cell, features = "MKI67", raster=FALSE)

immune_T_cell <- RenameIdents(immune_T_cell, '53' = "Proliferating T cells")

DimPlot(immune_T_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_T_cell, features = "FOXP3", raster=FALSE)

immune_T_cell <- RenameIdents(immune_T_cell, '35' = "Regulatory T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '17' = "Regulatory T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '25' = "Regulatory T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '40' = "Regulatory T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '63' = "Regulatory T cells")

DimPlot(immune_T_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_T_cell, features = "CD4", min.cutoff =1, max.cutoff =3, raster=FALSE)

immune_T_cell <- RenameIdents(immune_T_cell, '54' = "CD4 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '50' = "CD4 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '23' = "CD4 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '37' = "CD4 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '60' = "CD4 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '49' = "CD4 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '31' = "CD4 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '38' = "CD4 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '4' = "CD4 T cells")

immune_T_cell <- RenameIdents(immune_T_cell, '3' = "CD4 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '0' = "CD4 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '15' = "CD4 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '22' = "CD4 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '13' = "CD4 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '47' = "CD4 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '42' = "CD4 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '1' = "CD4 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '44' = "CD4 T cells")

DimPlot(immune_T_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_T_cell, features = "CD8A", min.cutoff =1, max.cutoff =3, raster=FALSE)

immune_T_cell <- RenameIdents(immune_T_cell, '19' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '26' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '30' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '18' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '14' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '46' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '56' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '59' = "CD8 T cells")

immune_T_cell <- RenameIdents(immune_T_cell, '5' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '65' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '41' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '7' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '32' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '51' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '61' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '28' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '9' = "CD8 T cells")

immune_T_cell <- RenameIdents(immune_T_cell, '16' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '43' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '24' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '33' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '36' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '6' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '2' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '45' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '48' = "CD8 T cells")

immune_T_cell <- RenameIdents(immune_T_cell, '27' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '39' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '62' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '67' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '57' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '29' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '64' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '8' = "CD8 T cells")
immune_T_cell <- RenameIdents(immune_T_cell, '34' = "CD8 T cells")


immune_T_cell[["immune_cluster"]]  <- immune_T_cell@active.ident


################################################################################# Remove Cells with no patient data

immune_T_cell[["sample"]]  <- immune_T_cell[["Patient_ID_SPC"]]
immune_T_cell[["group"]]  <- immune_T_cell[["Type"]]

md.immune_T_cell <- immune_T_cell@meta.data %>% as.data.table

length(which(md.immune_T_cell$sample==59 & md.immune_T_cell$group=='Tumor'& md.immune_T_cell$immune_cluster=='B cells'))

# Vector of sample numbers
sample_numbers <- c(4, 5, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 51, 52, 53, 54, 56, 59, 60, 62, 64, 65, 66, 68, 71, 73, 75)

# Vector of cluster names
cluster_vec <- c("CD4 T cells", "CD8 T cells", "Regulatory T cells", "NK cells", "gamma delta T cells", "MAIT cells", "Proliferating T cells")

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
    for (sample_num in sample_numbers) {
      # Calculate the length of the subset
      subset_length <- length(which(md.immune_T_cell$sample == sample_num & md.immune_T_cell$group == type & md.immune_T_cell$immune_cluster == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.immune_T_cell$sample == sample_num & md.immune_T_cell$group == type))
      subset_fraction <- subset_length / total_count
      # Add a row to the data frame
      new_row <- data.frame(sample_numbers = sample_num,
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

# Save the resulting data frame as a CSV file in the "D:/DATA/Gastric cancer/" directory
write.csv(output_df, file = "D:/DATA/Gastric cancer/cell_fraction_output.csv", row.names = FALSE)

################################################################################## Graph

# Create an empty list to store the plots
plot_list <- list()

for (cluster in cluster_vec) {
  # Subset the data for each type_vec
  subset_df <- output_df %>% subset(cluster_vec == cluster)
  normal_df <- subset_df %>% subset(type_vec == "Normal")
  tumor_df <- subset_df %>% subset(type_vec == "Tumor")
  # Perform the unpaired t-test
  ttest <- t.test(normal_df$subset_fraction, tumor_df$subset_fraction, var.equal=TRUE)
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
    geom_point(size = 4, shape =21, color = "black")+
    scale_fill_manual(values = c("grey", "brown")) +
    #geom_line(aes(group = sample_numbers), color = "grey") +
    labs(title = paste(cluster),
         x = "Type", y = "Subset Fraction") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          axis.ticks.length = unit(0.1, "cm"),
          #        axis.text.y = element_blank(),
          #        axis.text.x = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          legend.position = "none")+ 
    
    stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
    stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", size = 0.7, width=0.5, color="#03839C", alpha = 0.5)+
    annotate("text", x = 1.5, y = 0.9, label = paste("p =", format(ttest$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[cluster]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 2)

## ----


























################################################################################# subset-Tumor-Patient

immune_T_cell <- subset(x = immune_T_cell, subset = Patient_ID_SPC != 'FALSE')


cancer.subset <- subset(x = immune_T_cell, subset = Type == 'Tumor')

cancer.subset[["Molecular_type"]]  <- cancer.subset[["Patient_ID_SPC"]]
cancer.subset[["Stage"]]  <- cancer.subset[["Patient_ID_SPC"]]
cancer.subset[["Pathology_WHO"]]  <- cancer.subset[["Patient_ID_SPC"]]
cancer.subset[["Pathology_Lauren"]]  <- cancer.subset[["Patient_ID_SPC"]]

cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 71] <- "MSI"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 51] <- "MSI"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 35] <- "MSI"

cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 16] <- "CIN"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 50] <- "CIN"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 73] <- "CIN"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 53] <- "CIN"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 75] <- "CIN"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 23] <- "CIN"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 22] <- "CIN"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 27] <- "CIN"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 62] <- "CIN"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 40] <- "CIN"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 60] <- "CIN"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 15] <- "CIN"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 65] <- "CIN"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 5] <- "CIN"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 20] <- "CIN"

cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 7] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 4] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 29] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 28] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 8] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 18] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 68] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 36] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 44] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 56] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 52] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 66] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 64] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 54] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 59] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 41] <- "GS"

cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 44] <- "Stage II"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 68] <- "Stage II"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 71] <- "Stage II"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 73] <- "Stage II"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 18] <- "Stage II"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 56] <- "Stage II"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 15] <- "Stage II"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 40] <- "Stage II"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 52] <- "Stage II"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 60] <- "Stage II"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 75] <- "Stage II"

cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 28] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 35] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 4] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 7] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 29] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 53] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 54] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 8] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 16] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 23] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 27] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 36] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 20] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 22] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 41] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 50] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 59] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 65] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 66] <- "Stage III"

cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 5] <- "Stage IV"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 51] <- "Stage IV"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 62] <- "Stage IV"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 64] <- "Stage IV"


cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 23] <- "SRC"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 56] <- "SRC"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 41] <- "SRC"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 59] <- "SRC"

cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 65] <- "Mucinous AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 66] <- "Mucinous AD"

cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 44] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 68] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 71] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 73] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 18] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 15] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 40] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 52] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 60] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 75] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 28] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 35] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 4] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 7] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 29] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 53] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 54] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 8] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 16] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 27] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 36] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 20] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 22] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 50] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 5] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 51] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 62] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 64] <- "Tubular AD"



cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 23] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 56] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 41] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 59] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 65] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 66] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 44] <- "Intestinal"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 68] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 71] <- "Intestinal"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 73] <- "Intestinal"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 18] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 15] <- "Intestinal"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 40] <- "Intestinal"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 52] <- "Indeterminate"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 60] <- "Intestinal"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 75] <- "Intestinal"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 28] <- "Intestinal"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 35] <- "Intestinal"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 4] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 7] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 29] <- "Intestinal"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 53] <- "Intestinal"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 54] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 8] <- "Indeterminate"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 16] <- "Intestinal"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 27] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 36] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 20] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 22] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 50] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 5] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 51] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 62] <- "Intestinal"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 64] <- "Diffuse"

################################################################################## do_BarPlot


SCpubr::do_BarPlot(cancer.subset,
                   group.by = 'immune_cluster',
                   split.by = "Patient_ID_SPC",
                   position = "fill",
                   legend.position = "none",
                   flip = FALSE)

SCpubr::do_BarPlot(cancer.subset,
                   group.by = 'immune_cluster',
                   split.by = "Molecular_type",
                   position = "fill",
                   legend.position = "none",
                   flip = FALSE)

SCpubr::do_BarPlot(cancer.subset,
                   group.by = 'immune_cluster',
                   split.by = "Stage",
                   position = "fill",
                   legend.position = "none",
                   flip = FALSE)

SCpubr::do_BarPlot(cancer.subset,
                   group.by = 'immune_cluster',
                   split.by = "Outcome",
                   position = "fill",
                   legend.position = "none",
                   flip = FALSE)

SCpubr::do_BarPlot(cancer.subset,
                   group.by = 'immune_cluster',
                   split.by = "Pathology_WHO",
                   position = "fill",
                   legend.position = "none",
                   flip = FALSE)

SCpubr::do_BarPlot(cancer.subset,
                   group.by = 'immune_cluster',
                   split.by = "Pathology_Lauren",
                   position = "fill",
                   legend.position = "none",
                   flip = FALSE)

################################################################################## Dataframe for Molecular_type
md.cancer.subset <- cancer.subset@meta.data %>% as.data.table

length(which(md.cancer.subset$sample==4 & md.cancer.subset$Molecular_type=='GS'& md.cancer.subset$immune_cluster=='Regulatory CD4 T cells'))

# Vector of sample numbers
sample_numbers <- c(4, 5, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 51, 52, 53, 54, 56, 59, 60, 62, 64, 65, 66, 68, 71, 73, 75)

# Vector of cluster names
cluster_vec <- c("CD4 T cells", "CD8 T cells", "Regulatory T cells", "NK cells", "gamma delta T cells", "MAIT cells", "Proliferating T cells")

type_vec <- c("CIN", "GS", "MSI")

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
    for (sample_num in sample_numbers) {
      # Calculate the length of the subset
      subset_length <- length(which(md.cancer.subset$sample == sample_num & md.cancer.subset$Molecular_type == type & md.cancer.subset$immune_cluster == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.cancer.subset$sample == sample_num & md.cancer.subset$Molecular_type == type))
      subset_fraction <- subset_length / total_count
      # Add a row to the data frame
      new_row <- data.frame(sample_numbers = sample_num,
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

# Save the resulting data frame as a CSV file in the "D:/DATA/Gastric cancer/" directory
write.csv(output_df, file = "D:/DATA/Gastric cancer/cell_fraction_output(Molecular_type).csv", row.names = FALSE)



################################################################################## Graph

# Create an empty list to store the plots
plot_list <- list()

for (cluster in cluster_vec) {
  # Subset the data for each type_vec
  subset_df <- output_df %>% subset(cluster_vec == cluster)
  GS_df <- subset_df %>% subset(type_vec == "GS")
  CIN_df <- subset_df %>% subset(type_vec == "CIN")
  MSI_df <- subset_df %>% subset(type_vec == "MSI")  
  # Perform the unpaired t-test
  ttest_GS_CIN <- t.test(GS_df$subset_fraction, CIN_df$subset_fraction, var.equal=TRUE)
  ttest_GS_MSI <- t.test(GS_df$subset_fraction, MSI_df$subset_fraction, var.equal=TRUE)
  ttest_CIN_MSI <- t.test(CIN_df$subset_fraction, MSI_df$subset_fraction, var.equal=TRUE)
  
  
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
    geom_point(size = 4, shape =21, color = "black")+
    scale_fill_manual(values = c("blue", "grey", "brown")) +
    #geom_line(aes(group = sample_numbers), color = "grey") +
    labs(title = paste(cluster),
         x = "Type", y = "Subset Fraction") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          axis.ticks.length = unit(0.1, "cm"),
          #        axis.text.y = element_blank(),
          #        axis.text.x = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          legend.position = "none")+ 
    
    stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
    stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", size = 0.7, width=0.5, color="#03839C", alpha = 0.5)+
    annotate("text", x = 1.5, y = 0.6, label = paste("p =", format(ttest_GS_CIN$p.value, digits = 3)))+
    annotate("text", x = 2.5, y = 0.7, label = paste("p =", format(ttest_GS_MSI$p.value, digits = 3)))+
    annotate("text", x = 2, y = 0.9, label = paste("p =", format(ttest_CIN_MSI$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[cluster]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 2)


################################################################################## Dataframe for Stage
md.cancer.subset <- cancer.subset@meta.data %>% as.data.table

length(which(md.cancer.subset$sample==73 & md.cancer.subset$Stage=='Stage II'& md.cancer.subset$immune_cluster=='Regulatory CD4 T cells'))

# Vector of sample numbers
sample_numbers <- c(4, 5, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 51, 52, 53, 54, 56, 59, 60, 62, 64, 65, 66, 68, 71, 73, 75)

# Vector of cluster names
cluster_vec <- c("CD4 T cells", "CD8 T cells", "Regulatory T cells", "NK cells", "gamma delta T cells", "MAIT cells", "Proliferating T cells")

type_vec <- c("Stage II", "Stage III", "Stage IV")

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
    for (sample_num in sample_numbers) {
      # Calculate the length of the subset
      subset_length <- length(which(md.cancer.subset$sample == sample_num & md.cancer.subset$Stage == type & md.cancer.subset$immune_cluster == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.cancer.subset$sample == sample_num & md.cancer.subset$Stage == type))
      subset_fraction <- subset_length / total_count
      # Add a row to the data frame
      new_row <- data.frame(sample_numbers = sample_num,
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

# Save the resulting data frame as a CSV file in the "D:/DATA/Gastric cancer/" directory
write.csv(output_df, file = "D:/DATA/Gastric cancer/cell_fraction_output(Stage).csv", row.names = FALSE)



################################################################################## Graph

# Create an empty list to store the plots
plot_list <- list()

for (cluster in cluster_vec) {
  # Subset the data for each type_vec
  subset_df <- output_df %>% subset(cluster_vec == cluster)
  II_df <- subset_df %>% subset(type_vec == "Stage II")
  III_df <- subset_df %>% subset(type_vec == "Stage III")
  IV_df <- subset_df %>% subset(type_vec == "Stage IV")  
  # Perform the unpaired t-test
  ttest_II_III <- t.test(II_df$subset_fraction, III_df$subset_fraction, var.equal=TRUE)
  ttest_II_IV <- t.test(II_df$subset_fraction, IV_df$subset_fraction, var.equal=TRUE)
  ttest_III_IV <- t.test(III_df$subset_fraction, IV_df$subset_fraction, var.equal=TRUE)
  
  
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
    geom_point(size = 4, shape =21, color = "black")+
    scale_fill_manual(values = c("blue", "grey", "brown")) +
    #geom_line(aes(group = sample_numbers), color = "grey") +
    labs(title = paste(cluster),
         x = "Type", y = "Subset Fraction") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          axis.ticks.length = unit(0.1, "cm"),
          #        axis.text.y = element_blank(),
          #        axis.text.x = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          legend.position = "none")+ 
    
    stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
    stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", size = 0.7, width=0.5, color="#03839C", alpha = 0.5)+
    annotate("text", x = 1.5, y = 0.6, label = paste("p =", format(ttest_II_III$p.value, digits = 3)))+
    annotate("text", x = 2.5, y = 0.7, label = paste("p =", format(ttest_III_IV$p.value, digits = 3)))+
    annotate("text", x = 2, y = 0.9, label = paste("p =", format(ttest_II_IV$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[cluster]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 2)

################################################################################## Dataframe for Pathology_Lauren
md.cancer.subset <- cancer.subset@meta.data %>% as.data.table

length(which(md.cancer.subset$sample==68 & md.cancer.subset$Pathology_Lauren=='Diffuse'& md.cancer.subset$immune_cluster=='Regulatory CD4 T cells'))

# Vector of sample numbers
sample_numbers <- c(4, 5, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 51, 52, 53, 54, 56, 59, 60, 62, 64, 65, 66, 68, 71, 73, 75)

# Vector of cluster names
cluster_vec <- c("CD4 T cells", "CD8 T cells", "Regulatory T cells", "NK cells", "gamma delta T cells", "MAIT cells", "Proliferating T cells")

type_vec <- c("Diffuse", "Intestinal", "Indeterminate")

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
    for (sample_num in sample_numbers) {
      # Calculate the length of the subset
      subset_length <- length(which(md.cancer.subset$sample == sample_num & md.cancer.subset$Pathology_Lauren == type & md.cancer.subset$immune_cluster == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.cancer.subset$sample == sample_num & md.cancer.subset$Pathology_Lauren == type))
      subset_fraction <- subset_length / total_count
      # Add a row to the data frame
      new_row <- data.frame(sample_numbers = sample_num,
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

# Save the resulting data frame as a CSV file in the "D:/DATA/Gastric cancer/" directory
write.csv(output_df, file = "D:/DATA/Gastric cancer/cell_fraction_output(Pathology_Lauren).csv", row.names = FALSE)



################################################################################## Graph

# Create an empty list to store the plots
plot_list <- list()

for (cluster in cluster_vec) {
  # Subset the data for each type_vec
  subset_df <- output_df %>% subset(cluster_vec == cluster)
  Diff_df <- subset_df %>% subset(type_vec == "Diffuse")
  Inter_df <- subset_df %>% subset(type_vec == "Indeterminate")  
  Intesti_df <- subset_df %>% subset(type_vec == "Intestinal")
  
  # Perform the unpaired t-test
  ttest_Diff_Inter <- t.test(Diff_df$subset_fraction, Inter_df$subset_fraction, var.equal=TRUE)
  ttest_Inter_Intesti <- t.test(Inter_df$subset_fraction, Intesti_df$subset_fraction, var.equal=TRUE)
  ttest_Diff_Intesti <- t.test(Diff_df$subset_fraction, Intesti_df$subset_fraction, var.equal=TRUE)
  
  
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
    geom_point(size = 4, shape =21, color = "black")+
    scale_fill_manual(values = c("blue", "grey", "brown")) +
    #geom_line(aes(group = sample_numbers), color = "grey") +
    labs(title = paste(cluster),
         x = "Type", y = "Subset Fraction") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          axis.ticks.length = unit(0.1, "cm"),
          #        axis.text.y = element_blank(),
          #        axis.text.x = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          legend.position = "none")+ 
    
    stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
    stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", size = 0.7, width=0.5, color="#03839C", alpha = 0.5)+
    annotate("text", x = 1.5, y = 0.6, label = paste("p =", format(ttest_Diff_Inter$p.value, digits = 3)))+
    annotate("text", x = 2.5, y = 0.7, label = paste("p =", format(ttest_Inter_Intesti$p.value, digits = 3)))+
    annotate("text", x = 2, y = 0.9, label = paste("p =", format(ttest_Diff_Intesti$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[cluster]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 2)



################################################################################## Dataframe for Outcome
md.cancer.subset <- cancer.subset@meta.data %>% as.data.table

length(which(md.cancer.subset$sample==68 & md.cancer.subset$Outcome=='Diffuse'& md.cancer.subset$immune_cluster=='Regulatory CD4 T cells'))

# Vector of sample numbers
sample_numbers <- c(4, 5, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 51, 52, 53, 54, 56, 59, 60, 62, 64, 65, 66, 68, 71, 73, 75)

# Vector of cluster names
cluster_vec <- c("CD4 T cells", "CD8 T cells", "Regulatory T cells", "NK cells", "gamma delta T cells", "MAIT cells", "Proliferating T cells")

type_vec <- c("Recur", "Non_recur", "Meta")

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
    for (sample_num in sample_numbers) {
      # Calculate the length of the subset
      subset_length <- length(which(md.cancer.subset$sample == sample_num & md.cancer.subset$Outcome == type & md.cancer.subset$immune_cluster == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.cancer.subset$sample == sample_num & md.cancer.subset$Outcome == type))
      subset_fraction <- subset_length / total_count
      # Add a row to the data frame
      new_row <- data.frame(sample_numbers = sample_num,
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

# Save the resulting data frame as a CSV file in the "D:/DATA/Gastric cancer/" directory
write.csv(output_df, file = "D:/DATA/Gastric cancer/cell_fraction_output(Outcome).csv", row.names = FALSE)



################################################################################## Graph

# Create an empty list to store the plots
plot_list <- list()

for (cluster in cluster_vec) {
  # Subset the data for each type_vec
  subset_df <- output_df %>% subset(cluster_vec == cluster)
  R_df <- subset_df %>% subset(type_vec == "Recur")
  N_df <- subset_df %>% subset(type_vec == "Non_recur")  
  M_df <- subset_df %>% subset(type_vec == "Meta")
  
  # Perform the unpaired t-test
  ttest_R_N <- t.test(R_df$subset_fraction, N_df$subset_fraction, var.equal=TRUE)
  ttest_N_M <- t.test(N_df$subset_fraction, M_df$subset_fraction, var.equal=TRUE)
  ttest_R_M <- t.test(R_df$subset_fraction, M_df$subset_fraction, var.equal=TRUE)
  
  
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
    geom_point(size = 4, shape =21, color = "black")+
    scale_fill_manual(values = c("blue", "grey", "brown")) +
    #geom_line(aes(group = sample_numbers), color = "grey") +
    labs(title = paste(cluster),
         x = "Type", y = "Subset Fraction") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          axis.ticks.length = unit(0.1, "cm"),
          #        axis.text.y = element_blank(),
          #        axis.text.x = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          legend.position = "none")+ 
    
    stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
    stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", size = 0.7, width=0.5, color="#03839C", alpha = 0.5)+
    annotate("text", x = 1.5, y = 0.6, label = paste("p =", format(ttest_N_M$p.value, digits = 3)))+
    annotate("text", x = 2.5, y = 0.7, label = paste("p =", format(ttest_R_N$p.value, digits = 3)))+
    annotate("text", x = 2, y = 0.9, label = paste("p =", format(ttest_R_M$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[cluster]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 2)



################################################################################################### CD4 ----

CD4_cell <- subset(x = immune_T_cell, idents = "CD4 T cells")


CD4_cell <- RunPCA(CD4_cell, features = VariableFeatures(object = immune_no_stress))
CD4_cell <- FindNeighbors(CD4_cell, reduction = "pca", dims = 1:50) 
CD4_cell <- FindClusters(CD4_cell, resolution = 0.5)
CD4_cell <- RunUMAP(CD4_cell, reduction = "pca", dims = 1:50) 
DimPlot(CD4_cell, reduction = 'umap', raster=FALSE, label = TRUE, label.size = 4) + NoLegend()+NoAxes()















CD4_cell[["immune_cluster"]]  <- CD4_cell@active.ident



################################################################################# Remove Cells with no patient data

CD4_cell[["sample"]]  <- CD4_cell[["Patient_ID_SPC"]]
CD4_cell[["group"]]  <- CD4_cell[["Type"]]

md.CD4_cell <- CD4_cell@meta.data %>% as.data.table

length(which(md.CD4_cell$sample==59 & md.CD4_cell$group=='Tumor'& md.CD4_cell$immune_cluster=='B cells'))

# Vector of sample numbers
sample_numbers <- c(4, 5, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 51, 52, 53, 54, 56, 59, 60, 62, 64, 65, 66, 68, 71, 73, 75)

# Vector of cluster names
cluster_vec <- c(0, 1, 2, 3, 4, 5, 6, 7, 8)

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
    for (sample_num in sample_numbers) {
      # Calculate the length of the subset
      subset_length <- length(which(md.CD4_cell$sample == sample_num & md.CD4_cell$group == type & md.CD4_cell$immune_cluster == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.CD4_cell$sample == sample_num & md.CD4_cell$group == type))
      subset_fraction <- subset_length / total_count
      # Add a row to the data frame
      new_row <- data.frame(sample_numbers = sample_num,
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

# Save the resulting data frame as a CSV file in the "D:/DATA/Gastric cancer/" directory
write.csv(output_df, file = "D:/DATA/Gastric cancer/cell_fraction_output.csv", row.names = FALSE)

################################################################################## Graph

# Create an empty list to store the plots
plot_list <- list()

for (cluster in cluster_vec) {
  # Subset the data for each type_vec
  subset_df <- output_df %>% subset(cluster_vec == cluster)
  normal_df <- subset_df %>% subset(type_vec == "Normal")
  tumor_df <- subset_df %>% subset(type_vec == "Tumor")
  # Perform the unpaired t-test
  ttest <- t.test(normal_df$subset_fraction, tumor_df$subset_fraction, var.equal=TRUE)
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
    geom_point(size = 4, shape =21, color = "black")+
    scale_fill_manual(values = c("grey", "brown")) +
    #geom_line(aes(group = sample_numbers), color = "grey") +
    labs(title = paste(cluster),
         x = "Type", y = "Subset Fraction") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          axis.ticks.length = unit(0.1, "cm"),
          #        axis.text.y = element_blank(),
          #        axis.text.x = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          legend.position = "none")+ 
    
    stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
    stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", size = 0.7, width=0.5, color="#03839C", alpha = 0.5)+
    annotate("text", x = 1.5, y = 0.9, label = paste("p =", format(ttest$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[cluster]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 2)

## ----

DimPlot(CD4_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(CD4_cell, features = "CCR7", min.cutoff =1, max.cutoff =3,raster=FALSE)

DimPlot(CD4_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(CD4_cell, features = "TCF7", min.cutoff =1, max.cutoff =3,raster=FALSE)

DimPlot(CD4_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(CD4_cell, features = "SELL", min.cutoff =1, max.cutoff =3,raster=FALSE)

DimPlot(CD4_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(CD4_cell, features = "IL7R", min.cutoff =1, max.cutoff =3,raster=FALSE)




DimPlot(CD4_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(CD4_cell, features = "CXCR5", min.cutoff =1, max.cutoff =3,raster=FALSE)

DimPlot(CD4_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(CD4_cell, features = "BCL6", min.cutoff =1, max.cutoff =3,raster=FALSE)

DimPlot(CD4_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(CD4_cell, features = "CD200", min.cutoff =1, max.cutoff =3,raster=FALSE)

DimPlot(CD4_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(CD4_cell, features = "BTLA", min.cutoff =1, max.cutoff =3,raster=FALSE)

DimPlot(CD4_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(CD4_cell, features = "CXCL13", min.cutoff =1, max.cutoff =3,raster=FALSE)


DimPlot(CD4_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(CD4_cell, features = "IFNG", min.cutoff =1, max.cutoff =3,raster=FALSE)

DimPlot(CD4_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(CD4_cell, features = "CXCR3", min.cutoff =1, max.cutoff =3,raster=FALSE)

DimPlot(CD4_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(CD4_cell, features = "PDCD1", min.cutoff =1, max.cutoff =3,raster=FALSE)




DimPlot(CD4_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(CD4_cell, features = "IL17A", min.cutoff =1, max.cutoff =3,raster=FALSE)

DimPlot(CD4_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(CD4_cell, features = "IL22", min.cutoff =1, max.cutoff =3,raster=FALSE)

DimPlot(CD4_cell, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(CD4_cell, features = "IL23R", min.cutoff =1, max.cutoff =3,raster=FALSE)








CD8_cell <- subset(x = immune_T_cell, idents = "CD4 T cells")


CD8_cell <- RunPCA(CD4_cell, features = VariableFeatures(object = immune_no_stress))
CD8_cell <- FindNeighbors(CD8_cell, reduction = "pca", dims = 1:50) 
CD8_cell <- FindClusters(CD8_cell, resolution = 0.5)
CD8_cell <- RunUMAP(CD8_cell, reduction = "pca", dims = 1:50) 
DimPlot(CD8_cell, reduction = 'umap', raster=FALSE, label = TRUE, label.size = 4) + NoLegend()+NoAxes()








CD8_cell <- subset(x = immune_T_cell, idents = "CD8 T cells")














































DimPlot(immune_no_stress_PCA_1, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_no_stress_PCA_1, features = "NCAM1", min.cutoff =1, max.cutoff =3, raster=FALSE)

DimPlot(immune_no_stress_PCA_1, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_no_stress_PCA_1, features = "CD4", min.cutoff =1, max.cutoff =3, raster=FALSE)

DimPlot(immune_no_stress_PCA_1, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_no_stress_PCA_1, features = "CD8A", min.cutoff =1, max.cutoff =3, raster=FALSE)

FeaturePlot(immune_no_stress_PCA_1, features = "CD4", min.cutoff =1, max.cutoff =3, raster=FALSE) |
  FeaturePlot(immune_no_stress_PCA_1, features = "CD8A", min.cutoff =1, max.cutoff =3, raster=FALSE)

# B cells
#SCpubr::do_FeaturePlot(sample = immune_no_stress_PCA, features = "CD19", legend.position="none", min.cutoff =1, max.cutoff =3, pt.size = 2)
#SCpubr::do_FeaturePlot(sample = immune_no_stress_PCA, features = "MS4A1", legend.position="none", min.cutoff =1, max.cutoff =4, pt.size = 2)

immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '73' = "B cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '66' = "B cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '0' = "B cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, `8` = "B cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, `7` = "B cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, `79` = "B cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, `22` = "B cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, `24` = "B cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, `10` = "B cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, `58` = "B cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, `11` = "B cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, `53` = "B cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, `64` = "B cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, `71` = "Proliferating B cells")

#DimPlot(immune_no_stress_PCA, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes()


# Plasma cell
#SCpubr::do_FeaturePlot(sample = immune_no_stress_PCA, features = "SDC1", legend.position="none", min.cutoff =1, max.cutoff =3, pt.size = 2)
#SCpubr::do_FeaturePlot(sample = immune_no_stress_PCA, features = "JCHAIN", legend.position="none", min.cutoff =1, max.cutoff =4, pt.size = 2)

immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '62' = "Plasma cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '76' = "Plasma cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '21' = "Plasma cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, `77` = "Plasma cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, `60` = "Plasma cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, `55` = "Plasma cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, `49` = "Plasma cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, `68` = "Plasma cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, `57` = "Plasma cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, `58` = "Plasma cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, `40` = "Plasma cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, `67` = "Plasma cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, `12` = "Plasma cells")

#DimPlot(immune_no_stress_PCA, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes()


# Mast cell
#SCpubr::do_FeaturePlot(sample = immune_no_stress_PCA, features = "KIT", legend.position="none", min.cutoff =1, max.cutoff =3, pt.size = 2)
#SCpubr::do_FeaturePlot(sample = immune_no_stress_PCA, features = "IL1RL1", legend.position="none", min.cutoff =1, max.cutoff =4, pt.size = 2)

immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '50' = "Mast cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '31' = "Mast cells")

#DimPlot(immune_no_stress_PCA, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes()


# Myeloid cell
#SCpubr::do_FeaturePlot(sample = immune_no_stress_PCA, features = "ITGB2", legend.position="none", min.cutoff =0, max.cutoff =4, pt.size = 2)
#SCpubr::do_FeaturePlot(sample = immune_no_stress_PCA, features = "CD14", legend.position="none", min.cutoff =1, max.cutoff =4, pt.size = 2)
#SCpubr::do_FeaturePlot(sample = immune_no_stress_PCA, features = "FCN1", legend.position="none", min.cutoff =1, max.cutoff =4, pt.size = 2)
#SCpubr::do_FeaturePlot(sample = immune_no_stress_PCA, features = "CSF1R", legend.position="none", min.cutoff =1, max.cutoff =4, pt.size = 2)
#SCpubr::do_FeaturePlot(sample = immune_no_stress_PCA, features = "CD86", legend.position="none", min.cutoff =1, max.cutoff =10, pt.size = 2)
#SCpubr::do_FeaturePlot(sample = immune_no_stress_PCA, features = "CD1C", legend.position="none", min.cutoff =1, max.cutoff =15, pt.size = 2)

immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '78' = "Myeloid cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '37' = "Myeloid cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '44' = "Myeloid cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '45' = "Myeloid cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '23' = "Myeloid cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '42' = "Myeloid cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '19' = "Myeloid cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '34' = "Myeloid cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '39' = "Myeloid cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '74' = "Myeloid cells")

#DimPlot(immune_no_stress_PCA, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes()


# T/NK cell


result <- FindMarkers(immune_no_stress_PCA, ident.1 = "32", ident.2 = "8", min.pct = 0.5, logfc.threshold = log(2))


result[result$avg_log2FC > 0,] <- result[result$avg_log2FC > 0,][order(-result[result$avg_log2FC > 0,]$avg_log2FC),]

head(result, n=20)

DimPlot(immune_no_stress_PCA, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_no_stress_PCA, features = "TRAV1-2", min.cutoff =1, max.cutoff =3,raster=FALSE)

DimPlot(immune_no_stress_PCA, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_no_stress_PCA, features = "KLRB1", min.cutoff =1, max.cutoff =3,raster=FALSE)

immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '65' = "MAIT cells")


DimPlot(immune_no_stress_PCA, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_no_stress_PCA, features = "TRDV1", min.cutoff =1, max.cutoff =3, raster=FALSE)

DimPlot(immune_no_stress_PCA, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_no_stress_PCA, features = "TRDV2", min.cutoff =1, max.cutoff =3, raster=FALSE)

DimPlot(immune_no_stress_PCA, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_no_stress_PCA, features = "TRGC1", min.cutoff =1, max.cutoff =3, raster=FALSE)

DimPlot(immune_no_stress_PCA, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_no_stress_PCA, features = "TRGC2", min.cutoff =1, max.cutoff =3, raster=FALSE)

DimPlot(immune_no_stress_PCA, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_no_stress_PCA, features = "TRDC", min.cutoff =1, max.cutoff =3, raster=FALSE)


DimPlot(immune_no_stress_PCA, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_no_stress_PCA, features = "GNLY", min.cutoff =1, max.cutoff =3, raster=FALSE)

DimPlot(immune_no_stress_PCA, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_no_stress_PCA, features = "NCAM1", min.cutoff =1, max.cutoff =3, raster=FALSE)


FeaturePlot(immune_no_stress_PCA, features = "PDCD1", raster=FALSE)
FeaturePlot(immune_no_stress_PCA, features = "HAVCR2", raster=FALSE)
FeaturePlot(immune_no_stress_PCA, features = "CTLA4", raster=FALSE)
FeaturePlot(immune_no_stress_PCA, features = "LAG3", raster=FALSE)

FeaturePlot(immune_no_stress_PCA, features = "TNFRSF9", raster=FALSE)
FeaturePlot(immune_no_stress_PCA, features = "TNFRSF4", raster=FALSE)
FeaturePlot(immune_no_stress_PCA, features = "ICOS", raster=FALSE)

FeaturePlot(immune_no_stress_PCA, features = "CD69", raster=FALSE)
FeaturePlot(immune_no_stress_PCA, features = "ENTPD1", raster=FALSE)
FeaturePlot(immune_no_stress_PCA, features = "ITGAE", raster=FALSE)

FeaturePlot(immune_no_stress_PCA, features = "CD28", raster=FALSE, min.cutoff =0, max.cutoff =3)
FeaturePlot(immune_no_stress_PCA, features = "B3GAT1", raster=FALSE, min.cutoff =0, max.cutoff =3)

DimPlot(immune_no_stress_PCA, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
FeaturePlot(immune_no_stress_PCA, features = "GZMB", raster=FALSE, min.cutoff =0, max.cutoff =3)

FeaturePlot(immune_no_stress_PCA, features = "KIR2DL1", raster=FALSE, min.cutoff =0, max.cutoff =3)
FeaturePlot(immune_no_stress_PCA, features = "KIR2DL2", raster=FALSE, min.cutoff =0, max.cutoff =3)
FeaturePlot(immune_no_stress_PCA, features = "CCR7", raster=FALSE, min.cutoff =0, max.cutoff =3)
FeaturePlot(immune_no_stress_PCA, features = "SELL", raster=FALSE, min.cutoff =0, max.cutoff =3)




SCpubr::do_FeaturePlot(sample = immune_no_stress_PCA, features = "CD3G", legend.position="none", min.cutoff =0, max.cutoff =3, pt.size = 2)
SCpubr::do_FeaturePlot(sample = immune_no_stress_PCA, features = "CD3E", legend.position="none", min.cutoff =0, max.cutoff =3, pt.size = 2)

FeaturePlot(sample = immune_no_stress_PCA_1, features = "GNLY", legend.position="none", min.cutoff =0, max.cutoff =15, pt.size = 2)
FeaturePlot(sample = immune_no_stress_PCA_1, features = "NCAM1", legend.position="none", min.cutoff =6, max.cutoff =8, pt.size = 2)



immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '23' = "gamma delta T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '29' = "gamma delta T cells")

immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '31' = "NK cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '42' = "NK cells")



DimPlot(immune_no_stress_PCA, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_no_stress_PCA, features = "HSPB1", min.cutoff =1, max.cutoff =3, raster=FALSE)

DimPlot(immune_no_stress_PCA, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes() |
  FeaturePlot(immune_no_stress_PCA, features = "FOXP3", min.cutoff =1, max.cutoff =3, raster=FALSE)






#DimPlot(immune_no_stress_PCA, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes()


#SCpubr::do_FeaturePlot(sample = immune_no_stress_PCA, 
#                       idents.highlight = levels(immune_no_stress_PCA)[!(levels(immune_no_stress_PCA) 
#                                                                               %in% c('Dendritic cells','Monocytes','Plasma cells','B cells','Mast cells','Proliferating B cells'))], 
#                       features = c("MKI67"),min.cutoff =1, max.cutoff =10, pt.size = 2, legend.position="none")

#SCpubr::do_FeaturePlot(sample = immune_no_stress_PCA, 
#                       idents.highlight = levels(immune_no_stress_PCA)[!(levels(immune_no_stress_PCA) 
#                                                                               %in% c('Dendritic cells','Monocytes','Plasma cells','B cells','Mast cells','Proliferating B cells'))], 
#                       features = c("HSPB1"),min.cutoff =1, max.cutoff =4, pt.size = 2, legend.position="none")

immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '64' = "Proliferating T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '59' = "Stressed / dying T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '41' = "Stressed / dying T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '71' = "Stressed / dying T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '76' = "Stressed / dying T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '40' = "Stressed / dying T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '30' = "Stressed / dying T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '69' = "Stressed / dying T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '52' = "Stressed / dying T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '53' = "Stressed / dying T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '58' = "Stressed / dying T cells")

#DimPlot(immune_no_stress_PCA, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes()


# T cells specific

#SCpubr::do_FeaturePlot(sample = immune_no_stress_PCA, 
#                       idents.highlight = levels(immune_no_stress_PCA)[!(levels(immune_no_stress_PCA) 
#                                                                               %in% c('Dendritic cells','Monocytes','Plasma cells','B cells','Mast cells','Proliferating B cells'))], 
#                       features = c("CD4"),min.cutoff =1, max.cutoff =3, pt.size = 2, legend.position="none")

#SCpubr::do_FeaturePlot(sample = immune_no_stress_PCA, 
#                       idents.highlight = levels(immune_no_stress_PCA)[!(levels(immune_no_stress_PCA) 
#                                                                               %in% c('Dendritic cells','Monocytes','Plasma cells','B cells','Mast cells','Proliferating B cells'))], 
#                       features = c("CD8A"),min.cutoff =1, max.cutoff =4, pt.size = 2, legend.position="none")

#SCpubr::do_FeaturePlot(sample = immune_no_stress_PCA, 
#                       idents.highlight = levels(immune_no_stress_PCA)[!(levels(immune_no_stress_PCA) 
#                                                                               %in% c('Dendritic cells','Monocytes','Plasma cells','B cells','Mast cells','Proliferating B cells'))], 
#                       features = c("FOXP3"),min.cutoff =1, max.cutoff =4, pt.size = 2, legend.position="none")


immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '4' = "Regulatory CD4 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '33' = "Regulatory CD4 T cells")
#DimPlot(immune_no_stress_PCA, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes()

immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '49' = "CD8 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '47' = "CD8 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '72' = "CD8 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '46' = "CD8 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '3' = "CD8 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '7' = "CD8 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '8' = "CD8 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '61' = "CD8 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '16' = "CD8 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '21' = "CD8 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '19' = "CD8 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '14' = "CD8 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '28' = "CD8 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '67' = "CD8 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '10' = "CD8 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '63' = "CD8 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '13' = "CD8 T cells")

#DimPlot(immune_no_stress_PCA, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes()


immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '2' = "CD4 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '56' = "CD4 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '39' = "CD4 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '60' = "CD4 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '27' = "CD4 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '20' = "CD4 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '35' = "CD4 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '5' = "CD4 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '18' = "CD4 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '74' = "CD4 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '44' = "CD4 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '84' = "CD4 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '34' = "CD4 T cells")
immune_no_stress_PCA <- RenameIdents(immune_no_stress_PCA, '11' = "CD4 T cells")

DimPlot(immune_no_stress_PCA, reduction = 'umap', raster=FALSE, label = TRUE, label.size = 6) + NoLegend()+NoAxes() +RotatedAxis()

immune_no_stress_PCA[["immune_cluster"]]  <- immune_no_stress_PCA@active.ident

saveRDS(immune_no_stress_PCA, file="D:/DATA/Gastric cancer/immune_no_stress_PCA.rds")






































































































################################################################################# Dot plot

genes <- c("CD3G","CD3E","CD4","CD8A","CCR7","FOXP3","GNLY","NCAM1","NKG7","HSPB1","HSPH1", "GADD45B","MKI67",
           "CD86","CD1C","ITGB2","CD14","FCN1","CSF1R","KIT","IL1RL1","SDC1","JCHAIN","CD19","MS4A1")

SCpubr::do_ExpressionHeatmap(sample = immune_clusters_annotation,
                             features = genes,
                             flip = F,
                             cluster_cols = FALSE,
                             cluster_rows = TRUE,
                             enforce_symmetry = TRUE,
                             use_viridis = FALSE,
                             max.cutoff = 6,
                             min.cutoff = 0
)


SCpubr::do_DotPlot(sample = immune_clusters_annotation, 
                   features = genes)

DefaultAssay(immune_clusters_annotation) <-"ADT"
Proteins <- c('CD45','CD8','CD4','CD25','CD56','CD11c','CD14','CD16','CD19')

SCpubr::do_DotPlot(sample = immune_clusters_annotation, 
                   features = Proteins)

VlnPlot(immune_clusters_annotation, "CD4", raster=FALSE, pt.size = 0) + NoLegend()
VlnPlot(immune_clusters_annotation, "CD8", raster=FALSE, pt.size = 0) + NoLegend()
VlnPlot(immune_clusters_annotation, "CD25", raster=FALSE, pt.size = 0, y.max = 2) + NoLegend()
VlnPlot(immune_clusters_annotation, "CD56", raster=FALSE, pt.size = 0, y.max = 2) + NoLegend()

VlnPlot(immune_clusters_annotation, "CD14", raster=FALSE, pt.size = 0, y.max = 3) + NoLegend()
VlnPlot(immune_clusters_annotation, "CD19", raster=FALSE, pt.size = 0, y.max = 5) + NoLegend()

SCpubr::do_ViolinPlot(sample = immune_clusters_annotation,features = "CD4", plot_boxplot = FALSE)
SCpubr::do_ViolinPlot(sample = immune_clusters_annotation,features = "CD8", plot_boxplot = FALSE)
SCpubr::do_ViolinPlot(sample = immune_clusters_annotation,features = "CD25", plot_boxplot = FALSE, y_cut = rep(2, length(immune_clusters_annotation)))


DefaultAssay(immune_clusters_annotation) <-"integrated"


################################################################################# Tumor vs Normal 

SCpubr::do_BarPlot(sample = immune_clusters_annotation, 
                   group.by = "immune_cluster", 
                   legend.position = "none", 
                   plot.title = "Number of cells per cluster")


SCpubr::do_BarPlot(immune_clusters_annotation,
                   group.by = 'immune_cluster',
                   split.by = "Type",
                   plot.title = "Frequency of immune clusters",
                   position = "fill",
                   legend.position = "none",
                   flip = FALSE)


SCpubr::do_BarPlot(immune_clusters_annotation,
                   group.by = 'immune_cluster',
                   split.by = "Patient_ID",
                   plot.title = "Frequency of immune clusters",
                   position = "fill",
                   legend.position = "none",
                   flip = FALSE)



################################################################################# Remove Cells with no patient data

test.subset <- subset(x = immune_clusters_annotation, subset = Patient_ID_SPC != 'FALSE')

test.subset[["sample"]]  <- test.subset[["Patient_ID_SPC"]]
test.subset[["group"]]  <- test.subset[["Type"]]

md.test.subset <- test.subset@meta.data %>% as.data.table

length(which(md.test.subset$sample==59 & md.test.subset$group=='Tumor'& md.test.subset$immune_cluster=='B cells'))

# Vector of sample numbers
sample_numbers <- c(4, 5, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 51, 52, 53, 54, 56, 59, 60, 62, 64, 65, 66, 68, 71, 73, 75)

# Vector of cluster names
cluster_vec <- c("CD4 T cells", "CD8 T cells", "Regulatory CD4 T cells", "B cells", "Plasma cells", "Mast cells", "Monocytes", "Dendritic cells")

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
    for (sample_num in sample_numbers) {
      # Calculate the length of the subset
      subset_length <- length(which(md.test.subset$sample == sample_num & md.test.subset$group == type & md.test.subset$immune_cluster == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.test.subset$sample == sample_num & md.test.subset$group == type))
      subset_fraction <- subset_length / total_count
      # Add a row to the data frame
      new_row <- data.frame(sample_numbers = sample_num,
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

# Save the resulting data frame as a CSV file in the "D:/DATA/Gastric cancer/" directory
write.csv(output_df, file = "D:/DATA/Gastric cancer/cell_fraction_output.csv", row.names = FALSE)

################################################################################## Graph

# Create an empty list to store the plots
plot_list <- list()

for (cluster in cluster_vec) {
  # Subset the data for each type_vec
  subset_df <- output_df %>% subset(cluster_vec == cluster)
  normal_df <- subset_df %>% subset(type_vec == "Normal")
  tumor_df <- subset_df %>% subset(type_vec == "Tumor")
  # Perform the unpaired t-test
  ttest <- t.test(normal_df$subset_fraction, tumor_df$subset_fraction, var.equal=TRUE)
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
    geom_point(size = 4, shape =21, color = "black")+
    scale_fill_manual(values = c("grey", "brown")) +
    #geom_line(aes(group = sample_numbers), color = "grey") +
    labs(title = paste(cluster),
         x = "Type", y = "Subset Fraction") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          axis.ticks.length = unit(0.1, "cm"),
          #        axis.text.y = element_blank(),
          #        axis.text.x = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          legend.position = "none")+ 
    
    stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
    stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", size = 0.7, width=0.5, color="#03839C", alpha = 0.5)+
    annotate("text", x = 1.5, y = 0.9, label = paste("p =", format(ttest$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[cluster]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 2)



################################################################################# subset-Tumor-Patient

cancer.subset <- subset(x = test.subset, subset = Type == 'Tumor')

cancer.subset[["Molecular_type"]]  <- cancer.subset[["Patient_ID"]]
cancer.subset[["Stage"]]  <- cancer.subset[["Patient_ID"]]
cancer.subset[["Pathology_WHO"]]  <- cancer.subset[["Patient_ID"]]
cancer.subset[["Pathology_Lauren"]]  <- cancer.subset[["Patient_ID"]]

cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 71] <- "MSI"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 51] <- "MSI"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 35] <- "MSI"

cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 16] <- "CIN"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 50] <- "CIN"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 73] <- "CIN"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 53] <- "CIN"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 75] <- "CIN"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 23] <- "CIN"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 22] <- "CIN"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 27] <- "CIN"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 62] <- "CIN"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 40] <- "CIN"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 60] <- "CIN"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 15] <- "CIN"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 65] <- "CIN"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 5] <- "CIN"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 20] <- "CIN"

cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 7] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 4] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 29] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 28] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 8] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 18] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 68] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 36] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 44] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 56] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 52] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 66] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 64] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 54] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 59] <- "GS"
cancer.subset@meta.data["Molecular_type"][cancer.subset@meta.data["Molecular_type"] == 41] <- "GS"

cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 44] <- "Stage II"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 68] <- "Stage II"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 71] <- "Stage II"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 73] <- "Stage II"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 18] <- "Stage II"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 56] <- "Stage II"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 15] <- "Stage II"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 40] <- "Stage II"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 52] <- "Stage II"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 60] <- "Stage II"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 75] <- "Stage II"

cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 28] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 35] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 4] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 7] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 29] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 53] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 54] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 8] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 16] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 23] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 27] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 36] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 20] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 22] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 41] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 50] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 59] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 65] <- "Stage III"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 66] <- "Stage III"

cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 5] <- "Stage IV"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 51] <- "Stage IV"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 62] <- "Stage IV"
cancer.subset@meta.data["Stage"][cancer.subset@meta.data["Stage"] == 64] <- "Stage IV"


cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 23] <- "SRC"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 56] <- "SRC"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 41] <- "SRC"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 59] <- "SRC"

cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 65] <- "Mucinous AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 66] <- "Mucinous AD"

cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 44] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 68] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 71] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 73] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 18] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 15] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 40] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 52] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 60] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 75] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 28] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 35] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 4] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 7] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 29] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 53] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 54] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 8] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 16] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 27] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 36] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 20] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 22] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 50] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 5] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 51] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 62] <- "Tubular AD"
cancer.subset@meta.data["Pathology_WHO"][cancer.subset@meta.data["Pathology_WHO"] == 64] <- "Tubular AD"



cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 23] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 56] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 41] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 59] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 65] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 66] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 44] <- "Intestinal"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 68] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 71] <- "Intestinal"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 73] <- "Intestinal"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 18] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 15] <- "Intestinal"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 40] <- "Intestinal"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 52] <- "Indeterminate"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 60] <- "Intestinal"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 75] <- "Intestinal"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 28] <- "Intestinal"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 35] <- "Intestinal"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 4] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 7] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 29] <- "Intestinal"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 53] <- "Intestinal"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 54] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 8] <- "Indeterminate"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 16] <- "Intestinal"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 27] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 36] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 20] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 22] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 50] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 5] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 51] <- "Diffuse"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 62] <- "Intestinal"
cancer.subset@meta.data["Pathology_Lauren"][cancer.subset@meta.data["Pathology_Lauren"] == 64] <- "Diffuse"

################################################################################## do_BarPlot


SCpubr::do_BarPlot(cancer.subset,
                   group.by = 'immune_cluster',
                   split.by = "Patient_ID",
                   position = "fill",
                   legend.position = "none",
                   flip = FALSE)

SCpubr::do_BarPlot(cancer.subset,
                   group.by = 'immune_cluster',
                   split.by = "Molecular_type",
                   position = "fill",
                   legend.position = "none",
                   flip = FALSE)

SCpubr::do_BarPlot(cancer.subset,
                   group.by = 'immune_cluster',
                   split.by = "Stage",
                   position = "fill",
                   legend.position = "none",
                   flip = FALSE)

SCpubr::do_BarPlot(cancer.subset,
                   group.by = 'immune_cluster',
                   split.by = "Outcome",
                   position = "fill",
                   legend.position = "none",
                   flip = FALSE)

SCpubr::do_BarPlot(cancer.subset,
                   group.by = 'immune_cluster',
                   split.by = "Pathology_WHO",
                   position = "fill",
                   legend.position = "none",
                   flip = FALSE)

SCpubr::do_BarPlot(cancer.subset,
                   group.by = 'immune_cluster',
                   split.by = "Pathology_Lauren",
                   position = "fill",
                   legend.position = "none",
                   flip = FALSE)

################################################################################## Dataframe for Molecular_type
md.cancer.subset <- cancer.subset@meta.data %>% as.data.table

length(which(md.cancer.subset$sample==4 & md.cancer.subset$Molecular_type=='GS'& md.cancer.subset$immune_cluster=='Regulatory CD4 T cells'))

# Vector of sample numbers
sample_numbers <- c(4, 5, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 51, 52, 53, 54, 56, 59, 60, 62, 64, 65, 66, 68, 71, 73, 75)

# Vector of cluster names
cluster_vec <- c("CD4 T cells", "CD8 T cells", "Regulatory CD4 T cells", "B cells", "Plasma cells", "Mast cells", "Monocytes", "Dendritic cells")

type_vec <- c("CIN", "GS", "MSI")

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
    for (sample_num in sample_numbers) {
      # Calculate the length of the subset
      subset_length <- length(which(md.cancer.subset$sample == sample_num & md.cancer.subset$Molecular_type == type & md.cancer.subset$immune_cluster == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.cancer.subset$sample == sample_num & md.cancer.subset$Molecular_type == type))
      subset_fraction <- subset_length / total_count
      # Add a row to the data frame
      new_row <- data.frame(sample_numbers = sample_num,
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

# Save the resulting data frame as a CSV file in the "D:/DATA/Gastric cancer/" directory
write.csv(output_df, file = "D:/DATA/Gastric cancer/cell_fraction_output(Molecular_type).csv", row.names = FALSE)



################################################################################## Graph

# Create an empty list to store the plots
plot_list <- list()

for (cluster in cluster_vec) {
  # Subset the data for each type_vec
  subset_df <- output_df %>% subset(cluster_vec == cluster)
  GS_df <- subset_df %>% subset(type_vec == "GS")
  CIN_df <- subset_df %>% subset(type_vec == "CIN")
  MSI_df <- subset_df %>% subset(type_vec == "MSI")  
  # Perform the unpaired t-test
  ttest_GS_CIN <- t.test(GS_df$subset_fraction, CIN_df$subset_fraction, var.equal=TRUE)
  ttest_GS_MSI <- t.test(GS_df$subset_fraction, MSI_df$subset_fraction, var.equal=TRUE)
  ttest_CIN_MSI <- t.test(CIN_df$subset_fraction, MSI_df$subset_fraction, var.equal=TRUE)
  
  
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
    geom_point(size = 4, shape =21, color = "black")+
    scale_fill_manual(values = c("blue", "grey", "brown")) +
    #geom_line(aes(group = sample_numbers), color = "grey") +
    labs(title = paste(cluster),
         x = "Type", y = "Subset Fraction") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          axis.ticks.length = unit(0.1, "cm"),
          #        axis.text.y = element_blank(),
          #        axis.text.x = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          legend.position = "none")+ 
    
    stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
    stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", size = 0.7, width=0.5, color="#03839C", alpha = 0.5)+
    annotate("text", x = 1.5, y = 0.6, label = paste("p =", format(ttest_GS_CIN$p.value, digits = 3)))+
    annotate("text", x = 2.5, y = 0.7, label = paste("p =", format(ttest_GS_MSI$p.value, digits = 3)))+
    annotate("text", x = 2, y = 0.9, label = paste("p =", format(ttest_CIN_MSI$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[cluster]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 2)


################################################################################## Dataframe for Stage
md.cancer.subset <- cancer.subset@meta.data %>% as.data.table

length(which(md.cancer.subset$sample==73 & md.cancer.subset$Stage=='Stage II'& md.cancer.subset$immune_cluster=='Regulatory CD4 T cells'))

# Vector of sample numbers
sample_numbers <- c(4, 5, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 51, 52, 53, 54, 56, 59, 60, 62, 64, 65, 66, 68, 71, 73, 75)

# Vector of cluster names
cluster_vec <- c("CD4 T cells", "CD8 T cells", "Regulatory CD4 T cells", "B cells", "Plasma cells", "Mast cells", "Monocytes", "Dendritic cells")

type_vec <- c("Stage II", "Stage III", "Stage IV")

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
    for (sample_num in sample_numbers) {
      # Calculate the length of the subset
      subset_length <- length(which(md.cancer.subset$sample == sample_num & md.cancer.subset$Stage == type & md.cancer.subset$immune_cluster == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.cancer.subset$sample == sample_num & md.cancer.subset$Stage == type))
      subset_fraction <- subset_length / total_count
      # Add a row to the data frame
      new_row <- data.frame(sample_numbers = sample_num,
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

# Save the resulting data frame as a CSV file in the "D:/DATA/Gastric cancer/" directory
write.csv(output_df, file = "D:/DATA/Gastric cancer/cell_fraction_output(Stage).csv", row.names = FALSE)



################################################################################## Graph

# Create an empty list to store the plots
plot_list <- list()

for (cluster in cluster_vec) {
  # Subset the data for each type_vec
  subset_df <- output_df %>% subset(cluster_vec == cluster)
  II_df <- subset_df %>% subset(type_vec == "Stage II")
  III_df <- subset_df %>% subset(type_vec == "Stage III")
  IV_df <- subset_df %>% subset(type_vec == "Stage IV")  
  # Perform the unpaired t-test
  ttest_II_III <- t.test(II_df$subset_fraction, III_df$subset_fraction, var.equal=TRUE)
  ttest_II_IV <- t.test(II_df$subset_fraction, IV_df$subset_fraction, var.equal=TRUE)
  ttest_III_IV <- t.test(III_df$subset_fraction, IV_df$subset_fraction, var.equal=TRUE)
  
  
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
    geom_point(size = 4, shape =21, color = "black")+
    scale_fill_manual(values = c("blue", "grey", "brown")) +
    #geom_line(aes(group = sample_numbers), color = "grey") +
    labs(title = paste(cluster),
         x = "Type", y = "Subset Fraction") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          axis.ticks.length = unit(0.1, "cm"),
          #        axis.text.y = element_blank(),
          #        axis.text.x = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          legend.position = "none")+ 
    
    stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
    stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", size = 0.7, width=0.5, color="#03839C", alpha = 0.5)+
    annotate("text", x = 1.5, y = 0.6, label = paste("p =", format(ttest_II_III$p.value, digits = 3)))+
    annotate("text", x = 2.5, y = 0.7, label = paste("p =", format(ttest_III_IV$p.value, digits = 3)))+
    annotate("text", x = 2, y = 0.9, label = paste("p =", format(ttest_II_IV$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[cluster]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 2)

################################################################################## Dataframe for Pathology_Lauren
md.cancer.subset <- cancer.subset@meta.data %>% as.data.table

length(which(md.cancer.subset$sample==68 & md.cancer.subset$Pathology_Lauren=='Diffuse'& md.cancer.subset$immune_cluster=='Regulatory CD4 T cells'))

# Vector of sample numbers
sample_numbers <- c(4, 5, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 51, 52, 53, 54, 56, 59, 60, 62, 64, 65, 66, 68, 71, 73, 75)

# Vector of cluster names
cluster_vec <- c("CD4 T cells", "CD8 T cells", "Regulatory CD4 T cells", "B cells", "Plasma cells", "Mast cells", "Monocytes", "Dendritic cells")

type_vec <- c("Diffuse", "Intestinal", "Indeterminate")

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
    for (sample_num in sample_numbers) {
      # Calculate the length of the subset
      subset_length <- length(which(md.cancer.subset$sample == sample_num & md.cancer.subset$Pathology_Lauren == type & md.cancer.subset$immune_cluster == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.cancer.subset$sample == sample_num & md.cancer.subset$Pathology_Lauren == type))
      subset_fraction <- subset_length / total_count
      # Add a row to the data frame
      new_row <- data.frame(sample_numbers = sample_num,
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

# Save the resulting data frame as a CSV file in the "D:/DATA/Gastric cancer/" directory
write.csv(output_df, file = "D:/DATA/Gastric cancer/cell_fraction_output(Pathology_Lauren).csv", row.names = FALSE)



################################################################################## Graph

# Create an empty list to store the plots
plot_list <- list()

for (cluster in cluster_vec) {
  # Subset the data for each type_vec
  subset_df <- output_df %>% subset(cluster_vec == cluster)
  Diff_df <- subset_df %>% subset(type_vec == "Diffuse")
  Inter_df <- subset_df %>% subset(type_vec == "Indeterminate")  
  Intesti_df <- subset_df %>% subset(type_vec == "Intestinal")
  
  # Perform the unpaired t-test
  ttest_Diff_Inter <- t.test(Diff_df$subset_fraction, Inter_df$subset_fraction, var.equal=TRUE)
  ttest_Inter_Intesti <- t.test(Inter_df$subset_fraction, Intesti_df$subset_fraction, var.equal=TRUE)
  ttest_Diff_Intesti <- t.test(Diff_df$subset_fraction, Intesti_df$subset_fraction, var.equal=TRUE)
  
  
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
    geom_point(size = 4, shape =21, color = "black")+
    scale_fill_manual(values = c("blue", "grey", "brown")) +
    #geom_line(aes(group = sample_numbers), color = "grey") +
    labs(title = paste(cluster),
         x = "Type", y = "Subset Fraction") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          axis.ticks.length = unit(0.1, "cm"),
          #        axis.text.y = element_blank(),
          #        axis.text.x = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          legend.position = "none")+ 
    
    stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
    stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", size = 0.7, width=0.5, color="#03839C", alpha = 0.5)+
    annotate("text", x = 1.5, y = 0.6, label = paste("p =", format(ttest_Diff_Inter$p.value, digits = 3)))+
    annotate("text", x = 2.5, y = 0.7, label = paste("p =", format(ttest_Inter_Intesti$p.value, digits = 3)))+
    annotate("text", x = 2, y = 0.9, label = paste("p =", format(ttest_Diff_Intesti$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[cluster]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 2)



################################################################################## Dataframe for Outcome
md.cancer.subset <- cancer.subset@meta.data %>% as.data.table

length(which(md.cancer.subset$sample==68 & md.cancer.subset$Outcome=='Diffuse'& md.cancer.subset$immune_cluster=='Regulatory CD4 T cells'))

# Vector of sample numbers
sample_numbers <- c(4, 5, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 51, 52, 53, 54, 56, 59, 60, 62, 64, 65, 66, 68, 71, 73, 75)

# Vector of cluster names
cluster_vec <- c("CD4 T cells", "CD8 T cells", "Regulatory CD4 T cells", "B cells", "Plasma cells", "Mast cells", "Monocytes", "Dendritic cells")

type_vec <- c("Recur", "Non_recur", "Meta")

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
    for (sample_num in sample_numbers) {
      # Calculate the length of the subset
      subset_length <- length(which(md.cancer.subset$sample == sample_num & md.cancer.subset$Outcome == type & md.cancer.subset$immune_cluster == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.cancer.subset$sample == sample_num & md.cancer.subset$Outcome == type))
      subset_fraction <- subset_length / total_count
      # Add a row to the data frame
      new_row <- data.frame(sample_numbers = sample_num,
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

# Save the resulting data frame as a CSV file in the "D:/DATA/Gastric cancer/" directory
write.csv(output_df, file = "D:/DATA/Gastric cancer/cell_fraction_output(Outcome).csv", row.names = FALSE)



################################################################################## Graph

# Create an empty list to store the plots
plot_list <- list()

for (cluster in cluster_vec) {
  # Subset the data for each type_vec
  subset_df <- output_df %>% subset(cluster_vec == cluster)
  R_df <- subset_df %>% subset(type_vec == "Recur")
  N_df <- subset_df %>% subset(type_vec == "Non_recur")  
  M_df <- subset_df %>% subset(type_vec == "Meta")
  
  # Perform the unpaired t-test
  ttest_R_N <- t.test(R_df$subset_fraction, N_df$subset_fraction, var.equal=TRUE)
  ttest_N_M <- t.test(N_df$subset_fraction, M_df$subset_fraction, var.equal=TRUE)
  ttest_R_M <- t.test(R_df$subset_fraction, M_df$subset_fraction, var.equal=TRUE)
  
  
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
    geom_point(size = 4, shape =21, color = "black")+
    scale_fill_manual(values = c("blue", "grey", "brown")) +
    #geom_line(aes(group = sample_numbers), color = "grey") +
    labs(title = paste(cluster),
         x = "Type", y = "Subset Fraction") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          axis.ticks.length = unit(0.1, "cm"),
          #        axis.text.y = element_blank(),
          #        axis.text.x = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          legend.position = "none")+ 
    
    stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
    stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", size = 0.7, width=0.5, color="#03839C", alpha = 0.5)+
    annotate("text", x = 1.5, y = 0.6, label = paste("p =", format(ttest_N_M$p.value, digits = 3)))+
    annotate("text", x = 2.5, y = 0.7, label = paste("p =", format(ttest_R_N$p.value, digits = 3)))+
    annotate("text", x = 2, y = 0.9, label = paste("p =", format(ttest_R_M$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[cluster]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 2)

################################################################################## Save data ----

saveRDS(cancer.subset, file="D:/DATA/Gastric cancer/cancer_subset_ID.rds")


