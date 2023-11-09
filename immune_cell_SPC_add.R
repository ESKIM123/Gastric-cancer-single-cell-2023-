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
immune_clusters <- readRDS(file="D:/DATA/Gastric cancer/integrated_seurat_filtered.rds")
view(immune_clusters@meta.data)

# DefaultAssay(immune_clusters) <- "integrated"
immune_clusters <- RunUMAP(immune_clusters, reduction = "pca", dims = 1:50) 
immune_clusters <- FindNeighbors(immune_clusters, reduction = "pca", dims = 1:50) 
immune_clusters <- FindClusters(immune_clusters, resolution = 5)

DimPlot(immune_clusters, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes()
DimPlot(immune_clusters, reduction = 'umap', raster=FALSE, label = FALSE) + NoLegend() +NoAxes()

#immune_clusters_annotation <- immune_clusters

############################## Add Souporcell_output data ##############################
souporcell_output <- readRDS(file="D:/DATA/Gastric cancer/integrated_seurat_metadata_added_spc_info.rds")
view(souporcell_output)

# Add a new column in souporcell_output
souporcell_output$Patient_ID_SPC <- vector(length=nrow(souporcell_output))
colnames(souporcell_output)

# Add patient ID data in soupercell output
#


souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0516_1" | souporcell_output$Gem == "GC0516_2" | souporcell_output$Gem == "GC0516_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag1", 
                                           "28", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0516_1" | souporcell_output$Gem == "GC0516_2" | souporcell_output$Gem == "GC0516_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag2", 
                                           "35", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0516_1" | souporcell_output$Gem == "GC0516_2" | souporcell_output$Gem == "GC0516_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag3", 
                                           "44", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0516_1" | souporcell_output$Gem == "GC0516_2" | souporcell_output$Gem == "GC0516_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag4", 
                                           "68", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0516_1" | souporcell_output$Gem == "GC0516_2" | souporcell_output$Gem == "GC0516_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag5", 
                                           "71", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0516_1" | souporcell_output$Gem == "GC0516_2" | souporcell_output$Gem == "GC0516_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag6", 
                                           "73", souporcell_output$Patient_ID_SPC)

#

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0517_1" | souporcell_output$Gem == "GC0517_2" | souporcell_output$Gem == "GC0517_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag1", 
                                           "4", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0517_1" | souporcell_output$Gem == "GC0517_2" | souporcell_output$Gem == "GC0517_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag2", 
                                           "7", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0517_1" | souporcell_output$Gem == "GC0517_2" | souporcell_output$Gem == "GC0517_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag3", 
                                           "18", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0517_1" | souporcell_output$Gem == "GC0517_2" | souporcell_output$Gem == "GC0517_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag4", 
                                           "29", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0517_1" | souporcell_output$Gem == "GC0517_2" | souporcell_output$Gem == "GC0517_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag5", 
                                           "53", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0517_1" | souporcell_output$Gem == "GC0517_2" | souporcell_output$Gem == "GC0517_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag6", 
                                           "54", souporcell_output$Patient_ID_SPC)

#

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0530_1" | souporcell_output$Gem == "GC0530_2" | souporcell_output$Gem == "GC0530_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag1", 
                                           "8", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0530_1" | souporcell_output$Gem == "GC0530_2" | souporcell_output$Gem == "GC0530_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag2", 
                                           "16", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0530_1" | souporcell_output$Gem == "GC0530_2" | souporcell_output$Gem == "GC0530_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag3", 
                                           "23", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0530_1" | souporcell_output$Gem == "GC0530_2" | souporcell_output$Gem == "GC0530_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag4", 
                                           "27", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0530_1" | souporcell_output$Gem == "GC0530_2" | souporcell_output$Gem == "GC0530_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag5", 
                                           "36", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0530_1" | souporcell_output$Gem == "GC0530_2" | souporcell_output$Gem == "GC0530_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag6", 
                                           "56", souporcell_output$Patient_ID_SPC)

#

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0531_1" | souporcell_output$Gem == "GC0531_2" | souporcell_output$Gem == "GC0531_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag1", 
                                           "15", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0531_1" | souporcell_output$Gem == "GC0531_2" | souporcell_output$Gem == "GC0531_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag2", 
                                           "20", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0531_1" | souporcell_output$Gem == "GC0531_2" | souporcell_output$Gem == "GC0531_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag3", 
                                           "22", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0531_1" | souporcell_output$Gem == "GC0531_2" | souporcell_output$Gem == "GC0531_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag4", 
                                           "40", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0531_1" | souporcell_output$Gem == "GC0531_2" | souporcell_output$Gem == "GC0531_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag5", 
                                           "41", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0531_1" | souporcell_output$Gem == "GC0531_2" | souporcell_output$Gem == "GC0531_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag6", 
                                           "50", souporcell_output$Patient_ID_SPC)

#

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0711_1_1" | souporcell_output$Gem == "GC0711_1_2" | souporcell_output$Gem == "GC0711_1_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag1", 
                                           "52", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0711_1_1" | souporcell_output$Gem == "GC0711_1_2" | souporcell_output$Gem == "GC0711_1_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag2", 
                                           "59", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0711_1_1" | souporcell_output$Gem == "GC0711_1_2" | souporcell_output$Gem == "GC0711_1_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag3", 
                                           "60", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0711_1_1" | souporcell_output$Gem == "GC0711_1_2" | souporcell_output$Gem == "GC0711_1_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag4", 
                                           "65", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0711_1_1" | souporcell_output$Gem == "GC0711_1_2" | souporcell_output$Gem == "GC0711_1_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag5", 
                                           "66", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0711_1_1" | souporcell_output$Gem == "GC0711_1_2" | souporcell_output$Gem == "GC0711_1_3") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag6", 
                                           "75", souporcell_output$Patient_ID_SPC)



souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0711_2_1" | souporcell_output$Gem == "GC0711_2_2") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag1", 
                                           "5", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0711_2_1" | souporcell_output$Gem == "GC0711_2_2") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag2", 
                                           "51", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0711_2_1" | souporcell_output$Gem == "GC0711_2_2") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag3", 
                                           "62", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0711_2_1" | souporcell_output$Gem == "GC0711_2_2") & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag4", 
                                           "64", souporcell_output$Patient_ID_SPC)



souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0712_1_1" | souporcell_output$Gem == "GC0712_1_2" ) & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag1", 
                                           "7", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0712_1_1" | souporcell_output$Gem == "GC0712_1_2" ) & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag2", 
                                           "8", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0712_1_1" | souporcell_output$Gem == "GC0712_1_2" ) & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag3", 
                                           "18", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0712_1_1" | souporcell_output$Gem == "GC0712_1_2" ) & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag4", 
                                           "29", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0712_1_1" | souporcell_output$Gem == "GC0712_1_2" ) & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag5", 
                                           "36", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0712_1_1" | souporcell_output$Gem == "GC0712_1_2" ) & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag6", 
                                           "54", souporcell_output$Patient_ID_SPC)

#

souporcell_output$Patient_ID_SPC <- ifelse(souporcell_output$Gem == "GC0712_2" & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag1", 
                                           "56", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse(souporcell_output$Gem == "GC0712_2" & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag2", 
                                           "5", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse(souporcell_output$Gem == "GC0712_2"  & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag3", 
                                           "64", souporcell_output$Patient_ID_SPC)

#

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0718_1_1" | souporcell_output$Gem == "GC0718_1_2" ) & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag1", 
                                           "15", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0718_1_1" | souporcell_output$Gem == "GC0718_1_2" ) & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag2", 
                                           "22", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0718_1_1" | souporcell_output$Gem == "GC0718_1_2" ) & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag3", 
                                           "28", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0718_1_1" | souporcell_output$Gem == "GC0718_1_2" ) & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag4", 
                                           "41", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0718_1_1" | souporcell_output$Gem == "GC0718_1_2" ) & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag5", 
                                           "59", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse((souporcell_output$Gem == "GC0718_1_1" | souporcell_output$Gem == "GC0718_1_2" ) & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag6", 
                                           "65", souporcell_output$Patient_ID_SPC)


#

souporcell_output$Patient_ID_SPC <- ifelse(souporcell_output$Gem == "GC0718_2" & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag1", 
                                           "66", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse(souporcell_output$Gem == "GC0718_2" & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag2", 
                                           "68", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse(souporcell_output$Gem == "GC0718_2"  & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag3", 
                                           "71", souporcell_output$Patient_ID_SPC)

souporcell_output$Patient_ID_SPC <- ifelse(souporcell_output$Gem == "GC0718_2"  & 
                                             souporcell_output$spc_status == "singlet" & 
                                             souporcell_output$spc_assignment == "Hashtag4", 
                                           "73", souporcell_output$Patient_ID_SPC)

###################################### Integrate SPC data #################################
view(souporcell_output)
immune_clusters@meta.data <- cbind(souporcell_output, souporcell_output[, c("sample_info", "spc_status", "spc_assignment", "Patient_ID_SPC")])

view(immune_clusters@meta.data)

immune_clusters@meta.data <- immune_clusters@meta.data[, 1:(ncol(immune_clusters@meta.data)-1)]

# Count SPC Hashtag data
ncol(immune_clusters)
unique(immune_clusters@meta.data$spc_status)

sum(immune_clusters@meta.data$spc_status == "singlet")
sum(immune_clusters@meta.data$spc_status == "doublet")
sum(immune_clusters@meta.data$spc_status == "unassigned")

sum(immune_clusters@meta.data$Patient_ID_SPC == "FALSE")

# Compare with original Hashtag data
sum(is.na(immune_clusters@meta.data$Patient_ID))

###################################### Save data #################################

saveRDS(immune_clusters, file="D:/DATA/Gastric cancer/integrated_seurat_filtered_SPC_ID.rds")

###################################### Remove doublet data #################################

immune_clusters_singlet <- subset(x = immune_clusters, spc_status == "singlet" | spc_status == "unassigned")

nrow(immune_clusters_singlet@meta.data)
ncol(immune_clusters_singlet)

immune_clusters_singlet <- RunPCA(immune_clusters_singlet, verbose = FALSE)
immune_clusters_singlet <- FindNeighbors(immune_clusters_singlet, reduction = "pca", dims = 1:50) 
immune_clusters_singlet <- FindClusters(immune_clusters_singlet, resolution = 0.4)
immune_clusters_singlet <- RunUMAP(immune_clusters_singlet, reduction = "pca", dims = 1:50) 
DimPlot(immune_clusters_singlet, reduction = 'umap', raster=FALSE) + NoLegend()

###################################### Save data #################################
saveRDS(immune_clusters_singlet, file="D:/DATA/Gastric cancer/integrated_seurat_filtered_SPC_ID_singlet.rds")


