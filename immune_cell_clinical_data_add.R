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
immune_annotation <- readRDS(file="D:/DATA/Gastric cancer/immune_cell.rds")
DimPlot(immune_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Initial")

cancer.subset <- immune_annotation

cancer.subset[["Molecular_type"]]  <- cancer.subset[["Patient_ID_SPC"]]
cancer.subset[["Stage"]]  <- cancer.subset[["Patient_ID_SPC"]]
cancer.subset[["Pathology_WHO"]]  <- cancer.subset[["Patient_ID_SPC"]]
cancer.subset[["Pathology_Lauren"]]  <- cancer.subset[["Patient_ID_SPC"]]
cancer.subset[["MSI_type"]]  <- cancer.subset[["Patient_ID_SPC"]]




cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 71] <- "MSI"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 51] <- "MSI"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 35] <- "MSI"

cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 16] <- "MSS"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 50] <- "MSS"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 73] <- "MSS"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 53] <- "MSS"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 75] <- "MSS"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 23] <- "MSS"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 22] <- "MSS"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 27] <- "MSS"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 62] <- "MSS"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 40] <- "MSS"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 60] <- "MSS"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 15] <- "MSS"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 65] <- "MSS"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 5] <- "MSS"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 20] <- "MSS"

cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 7] <- "MSS"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 4] <- "MSS"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 29] <- "MSS"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 28] <- "MSS"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 8] <- "MSS"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 18] <- "MSS"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 68] <- "MSS"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 36] <- "MSS"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 44] <- "MSS"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 56] <- "MSS"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 52] <- "MSS"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 66] <- "MSS"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 64] <- "MSS"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 54] <- "MSS"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 59] <- "MSS"
cancer.subset@meta.data["MSI_type"][cancer.subset@meta.data["MSI_type"] == 41] <- "MSS"







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

view(cancer.subset@meta.data)

##################################### Save data: annotation #########################################

saveRDS(cancer.subset, file="D:/DATA/Gastric cancer/immune_cell_clinical_data_add.rds")
