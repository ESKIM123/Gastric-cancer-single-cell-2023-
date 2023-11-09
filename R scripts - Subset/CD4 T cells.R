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
library(decoupleR)
library(OmnipathR)
library(dorothea)
library(Trex)


# Get data ----
immune_annotation <- readRDS(file="D:/DATA/Gastric cancer/immune_cell_annotation_final_OLD.rds")
DimPlot(immune_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Initial")

CD4_subset <- subset(x = immune_annotation, idents = c("CD4 T cells", "Regulatory CD4 T cells"))

##################################### CD4 T cells clustering ######################################### OLD
CD4_subset <- RunPCA(CD4_subset)
CD4_subset <- FindNeighbors(CD4_subset, reduction = "pca", dims = 1:20)
CD4_subset <- FindClusters(CD4_subset, resolution = 0.25)
CD4_subset <- RunUMAP(CD4_subset, reduction = "pca", dims = 1:20)

DimPlot(CD4_subset, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("CD4 T cells")

CD4_subset_anno <- subset(x = CD4_subset, seurat_clusters != 6)
CD4_subset_anno <- subset(x = CD4_subset_anno, seurat_clusters != 5)
CD4_subset_anno <- RunUMAP(CD4_subset_anno, reduction = "pca", dims = 1:20)
DimPlot(CD4_subset_anno, reduction = 'umap', raster=FALSE, label = F) + ggtitle("")
DimPlot(CD4_subset_anno, reduction = 'umap', raster=FALSE, label = F) + ggtitle("")+ NoLegend() +NoAxes()

CD4_subset_anno <- RenameIdents(CD4_subset_anno, `0`  = 'C0')
CD4_subset_anno <- RenameIdents(CD4_subset_anno, `1`  = 'C1')
CD4_subset_anno <- RenameIdents(CD4_subset_anno, `2`  = 'C2')
CD4_subset_anno <- RenameIdents(CD4_subset_anno, `3`  = 'C3')
CD4_subset_anno <- RenameIdents(CD4_subset_anno, `4`  = 'C4')

view(CD4_subset_anno@meta.data)
CD4_subset_anno

CD4_subset_anno@meta.data$seurat_clusters <- CD4_subset_anno@active.ident

saveRDS(CD4_subset_anno, file="D:/DATA/Gastric cancer/GS_subsets/CD4_subset.rds")
CD4_subset_anno <- readRDS(file="D:/DATA/Gastric cancer/GS_subsets/CD4_subset.rds")
DimPlot(CD4_subset_anno, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("CD4 T cells") 
DimPlot(CD4_subset_anno, reduction = 'umap', raster=FALSE, label = F) + ggtitle("")+NoLegend()+NoAxes()

CD4_subset_anno


##################################### CD4 T cells clustering ######################################### NEW
CD4_subset_filtered_Trex <- RunPCA(CD4_subset_filtered_Trex)
CD4_subset_filtered_Trex <- FindNeighbors(CD4_subset_filtered_Trex, reduction = "pca", dims = 1:50)
CD4_subset_filtered_Trex <- FindClusters(CD4_subset_filtered_Trex, resolution = 0.3)
CD4_subset_filtered_Trex <- RunUMAP(CD4_subset_filtered_Trex, reduction = "pca", dims = 1:50)

DimPlot(CD4_subset_filtered_Trex, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("CD4 T cells") 



############################################################################## Barplot #####

dittoBarPlot(CD4_subset_anno, "seurat_clusters", group.by = "seurat_clusters", scale = "count",  retain.factor.levels = F, main = "")+ 
  scale_y_continuous(expand = c(0, 0)) + NoLegend()+ xlab(NULL)+ ylab(NULL)

cluster0.markers <- FindMarkers(CD4_subset_anno, ident.1 = "C0", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(CD4_subset_anno, features = "CD8A", split.by = "seurat_clusters", pt.size = 0)


signatures <- list(
  Naive = c("CCR7", "IL7R", "TCF7", "SELL"),
  Tfh = c("STAT3", "BCL6", "CXCR5", "CXCR4", "MAF"),
  Th1 = c("IFNG", "TBX21", "STAT4", "IL12RB2", "TNF", "IFNGR1", "STAT1"),
  Th17 =c("IL23R", "RORC", "TGFBR1", "IL17A", "IL22", "STAT3"),
  Treg = c("FOXP3", "IL2RA", "IL10")
)
CD4_subset_anno <- AddModuleScore_UCell(CD4_subset_anno,features=signatures, name=NULL)
CD4_subset_anno <- SmoothKNN(CD4_subset_anno, signature.names = names(signatures), reduction="pca")


FeaturePlot(CD4_subset_anno, features = "Naive_kNN", min.cutoff = "q30", max.cutoff = "q90")
FeaturePlot(CD4_subset_anno, features = "Tfh_kNN", min.cutoff = "q30", max.cutoff = "q90")
FeaturePlot(CD4_subset_anno, features = "Th1_kNN", min.cutoff = "q30", max.cutoff = "q90")
FeaturePlot(CD4_subset_anno, features = "Th17_kNN", min.cutoff = "q30", max.cutoff = "q90")
FeaturePlot(CD4_subset_anno, features = "Treg_kNN", min.cutoff = "q30", max.cutoff = "q90")

FeaturePlot(CD4_subset_anno, features = "CXCR3", min.cutoff = "q70", max.cutoff = "q90", blend = F, cols = c("lightgrey", "red"))




# Dot plots
genes <- c("IFNG", "TBX21", "CXCR3", "IL12RB2", "TGFB1", "IFNGR1", "STAT4",
           "CCR7", "IL7R", "TCF7", "SELL",
           "IL23R", "RORC", "TGFBR1", "IL17A", "IL22", "STAT3",
           "FOXP3", "IL2RA", "IL10", "TNFRSF9", "CTLA4",
           "BCL6", "CXCR5", "CXCR4", "MAF")

DotPlot(CD4_subset_anno, features = genes)+ RotatedAxis() + xlab(NULL) + ylab(NULL) + theme(axis.text.x = element_text(face = "italic"))
DotPlot(CD4_subset_anno, features = genes)+ RotatedAxis() + xlab(NULL) + ylab(NULL) +NoLegend() + 
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))

do_ExpressionHeatmap(sample = CD4_subset_anno,
                     features = genes,
                     viridis_direction = -1, max.cutoff = 1.5)























































CD4_subset_anno <- CD4_subset

# Markers
CD4_markers <- c("CCR7", "TCF7", "SELL", 
                 "ISG15", "IFIT3", "STAT1",
                 "GNLY", "PRF1",
                 "ENTPD1", "ITGAE",
                 "PDCD1", "TIGIT", "CTLA4",
                 "FOXP3", "EOMES", "TBX21", "GZMK", "GZMB"
)


# Dot plots
DotPlot(CD4_subset_anno, features = CD4_markers, col.min = 0,  col.max = 2)+ RotatedAxis() + xlab(NULL) + ylab(NULL) + theme(axis.text.x = element_text(face = "italic"))
SCpubr::do_ExpressionHeatmap(sample = CD4_subset_anno,
                             features = CD4_markers,
                             flip = F,
                             cluster_cols = F,
                             cluster_rows = F,
                             enforce_symmetry = TRUE,
                             use_viridis = FALSE,
                             max.cutoff = 3,
                             min.cutoff = 0)

signatures <- list(
  Naive_like_cells = c("TCF7", "LEF1", "TXK", "KLF2", "SCML1", "CCR7", "IL7R", "EEF1A1", "SELL", "EEF1B2", "PLAC8", "TPT1", "ACTN1", "GIMAP7", "MAL", "LINC00861", "FHIT"),
  ISG_cells = c("PLSCR1", "IRF7", "STAT1", "SP100", "IFI16", "TNFSF10", "TNFSF13B", "IFIT3", "IFIT1", "RSAD2", "IFI44L", "MX1", "ISG15", "IFI6", "MX2", "IFIT2", "CMPK2"),
  CCR6_Th17_cells = c("BHLHE40", "HOPX", "ZFP36L2", "MYBL1", "RORA", "TNFSF13B", "CCL20", "0LG", "CCL5", "TIMP1", "FLT3LG", "IL7R", "CXCR4", "CCR6", "IFNGR1", 
                      "GPR35", "ANXA1", "KLRB1", "IL4I1", "VIM", "FKBP11", "CAPG", "ABCB1", "TPT1", "PERP", "DPP4"),
  IL26_Th17_cells = c("RORC", "RORA", "CEBPD", "RUNX2", "ID2", "IL26", "IL17F", "IL17A", "TNFSF13B", "CCL20", "0LG", "GZMA", "CKLF", "IL22", "TNFSF14", "CXCR6", 
                      "CCR6", "IL23R", "IL17RE", "IL7R", "KLRB1", "CTSH", "CA10", "CAPG", "PTPN13", "MGAT4A", "ERN1", "IL4I1", "ANKRD28", "TMIGD2"),
  GZMK_Tem_cells = c("EOMES", "LITAF", "NFATC2", "RUNX3", "NR4A2", "GZMK", "CCL4", "GZMA", "CCL5", "GZMH", "TNFSF9", "CCL3L3", "IFNG", "FASLG", "CCL3", "CD74", "CCR5", 
                     "CXCR4", "CXCR3", "NKG7", "CST7", "CRTAM", "CCL4L2", "ENC1", "ITM2C", "SLAMF7", "AOAH", "F2R", "DTHD1"),
  IL21_Tfh_cells = c("NR3C1", "TOX2", "TOX", "TSHZ2", "RBPJ", "CXCL13", "IL21", "TNFSF8", "TNFSF11", "TNFSF4", "LIF", "0LG", "IL6ST", "CXCR5", "IL6R", "GNG4", "CD200",
                     "NMB", "IGFL2", "LHFP", "CPM", "BTLA", "FKBP5", "PDE7B", "ITM2A"),
  IFNG_Tfh_Th1_cells = c("ZBED2", "ZEB2", "ID2", "BHLHE40", "RBPJ", "CXCL13", "CCL4", "IFNG", "GZMB", "CCL3", "GZMA", "CCL5", "CSF2", "IL21", "FAM3C", "CXCR6", "CCR5", 
                         "", "IL2RG", "CD74", "LAG3", "HAVCR2", "PDCD1", "KRT86", "MYO7A", "PTMS", "RDH10", "CCL4L2", "PDE7B", "DUSP4"),
  TEMRA_cells = c("GZMA", "XCL2", "GZMM", "CCL3", "FASLG", "GZMK", "CX3CR1", "CMKLR1", "IL5RA", "CXCR2", "IL18RAP", "FGFBP2", "NKG7", "GNLY", "S1PR5", "KLRD1", "C1orf21", 
                  "PRF1", "PLEK", "CTSW", "PRSS23"),
  CREM_Tm_cells = c("CREM","ZFP36L2","ZNF331","NR4A2","SKIL", "CCL5","GZMA","TGFB1","AREG","CXCL16","GZMK", "CXCR4","IFNGR1","GPR35","4","IL7R",
                    "ZFP36","FTH1","SRGN","FAM177A1","TNFAIP3","YPEL5","LMNA","ANXA1","PTGER4","CDKN1A"),
  CCL5_Tm_cells = c("ZFP36L2", "HOPX", "ID2", "BHLHE40", "NR4A2", "CCL5", "GZMA", "GZMK", "XCL1", "0LG", "CCL4", "XCL2", "GZMH", "IFNG",
                    "GZMM", "CXCR4", "IL7R", "IFNGR1", "4", "IL18RAP", "ANXA1", "ALOX5AP", "GLUL", "ITGA1", "PARP8", "PTGER4", "CD69", "ZFP36", "CD99", "CLU"),
  CAPG_Tm_cells = c("HOPX", "BHLHE40", "ID2", "RBPJ", "PPARG", "GZMA", "CCL5", "CKLF", "TIMP1", "IL32", "TGFB1", "0LG", "HMGB1", "CSF2", "CXCR3", "CD74", 
                    "IL2RG", "CXCR6", "SH3BGRL3", "TMSB4X", "ACTB", "CAPG", "ALOX5AP", "CD52", "PFN1", "S100A11", "S100A4", "COTL1"),
  Th2_cells = c("AREG", "IL4", "IL5", "IL13","CXCR4", "CCR4", "CCR8", "PTGDR2", "HAVCR1", "IL17RB", "IL1RL1","BATF", "GATA3", "IRF4", "STAT6"),
  Tfh_cells = c("IL21","CXCR3", "CXCR5", "ICOS", "PDCD1","BATF", "BCL6", "MAF", "IRF4", "STAT3"), 
  Cytotoxic = c("GZMH", "GZMB", "GZMA", "PRF1"), 
  Exhaustion = c("PDCD1", "HAVCR2", "LAG3", "TIGIT")
)

CD4_subset_anno <- AddModuleScore_UCell(CD4_subset_anno,features=signatures, name=NULL)
CD4_subset_anno <- SmoothKNN(CD4_subset_anno, signature.names = names(signatures), reduction="pca")

FeaturePlot(CD4_subset_anno, reduction = "umap", features = "Cytotoxic_kNN")
VlnPlot(CD4_subset_anno, features = "Cytotoxic_kNN", pt.size = 0, split.by = "seurat_clusters")+ NoLegend()

FeaturePlot(CD4_subset_anno, reduction = "umap", features = "Exhaustion_kNN")
VlnPlot(CD4_subset_anno, features = "Exhaustion_kNN", pt.size = 0, split.by = "seurat_clusters")+ NoLegend()


CD4.markers <- FindAllMarkers(CD4_subset_anno, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CD4.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(CD4_subset_anno, features = top10$gene) + NoLegend()


##################################### Stats #########################################

CD4_subset_cancer <- subset(x = CD4_subset_anno, Type == 'Tumor' & Outcome != 'Meta' & Patient_ID_SPC != F)
CD4_subset_cancer_MSI <- subset(x = CD4_subset_cancer, MSI_type == 'MSI')
CD4_subset_cancer_MSS <- subset(x = CD4_subset_cancer, MSI_type == 'MSS')

SCpubr::do_BarPlot(sample = CD4_subset_cancer, 
                   group.by = "seurat_clusters", 
                   legend.position = "none", 
                   plot.title = "",  xlab = "",  ylab = "Number of cells per cluster", font.size = 11)


# 전체
SCpubr::do_BarPlot(sample = CD4_subset_cancer, 
                   group.by = "seurat_clusters", 
                   split.by = "Patient_ID_SPC",
                   position = "Fill",legend.position = "bottom",
                   plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)

# 부분
SCpubr::do_BarPlot(sample = CD4_subset_cancer_MSI, 
                   group.by = "seurat_clusters", 
                   split.by = "Patient_ID_SPC",
                   position = "Fill",legend.position = "bottom",
                   plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)

SCpubr::do_BarPlot(sample = CD4_subset_cancer_MSS, 
                   group.by = "seurat_clusters", 
                   split.by = "Patient_ID_SPC",
                   position = "Fill",legend.position = "bottom",
                   plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)
