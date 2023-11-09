
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
immune_annotation <- readRDS(file="D:/DATA/Gastric cancer/immune_cell_annotation_final_OLD.rds")
DimPlot(immune_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Initial")

TCR.subset <- immune_annotation

TCR.subset[["TCR_type"]]  <- TCR.subset[["Patient_ID_SPC"]]


TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 22] <- "C0_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 59] <- "C0_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 8] <- "C0_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 66] <- "C0_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 28] <- "C0_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 18] <- "C0_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 65] <- "C0_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 16] <- "C0_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 44] <- "C0_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 53] <- "C0_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 56] <- "C0_EXP"

TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 73] <- "C1_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 15] <- "C1_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 71] <- "C1_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 54] <- "C1_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 23] <- "C1_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 60] <- "C1_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 20] <- "C1_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 36] <- "C1_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 40] <- "C1_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 29] <- "C1_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 68] <- "C1_EXP"

TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 52] <- "C2_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 35] <- "C2_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 50] <- "C2_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 7] <- "C2_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 27] <- "C2_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 75] <- "C2_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 4] <- "C2_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 41] <- "C2_EXP"

view(TCR.subset@meta.data)



saveRDS(TCR.subset, file="D:/DATA/Gastric cancer/immune_cell_annotation_final_OLD_TCRtype.rds")


























