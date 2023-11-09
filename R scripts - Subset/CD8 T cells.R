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



# Get data ----
immune_annotation <- readRDS(file="D:/DATA/Gastric cancer/immune_cell_annotation_final_OLD.rds")
DimPlot(immune_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Initial")

CD8_subset <- subset(x = immune_annotation, ident = "CD8 T cells")

##################################### CD8 T cells clustering ######################################### OLD
CD8_subset <- RunPCA(CD8_subset)
CD8_subset <- FindNeighbors(CD8_subset, reduction = "pca", dims = 1:20)
CD8_subset <- FindClusters(CD8_subset, resolution = 0.25)
CD8_subset <- RunUMAP(CD8_subset, reduction = "pca", dims = 1:20)

DimPlot(CD8_subset, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("CD8 T cells")
VlnPlot(CD8_subset, features = "CD8A", pt.size = 0, split.by = "seurat_clusters")+ NoLegend()

CD8_subset_filtered <- subset(x = CD8_subset, seurat_clusters != 2)
CD8_subset_filtered <- subset(x = CD8_subset_filtered, seurat_clusters != 7)
CD8_subset_filtered <- subset(x = CD8_subset_filtered, seurat_clusters != 5)

CD8_subset_filtered <- RunUMAP(CD8_subset_filtered, reduction = "pca", dims = 1:20)
DimPlot(CD8_subset_filtered, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("CD8 T cells")

CD8_subset_anno <- CD8_subset_filtered 

CD8_subset_anno <- RenameIdents(CD8_subset_anno, `0`  = 'C0')
CD8_subset_anno <- RenameIdents(CD8_subset_anno, `1`  = 'C1')
CD8_subset_anno <- RenameIdents(CD8_subset_anno, `3`  = 'C2')
CD8_subset_anno <- RenameIdents(CD8_subset_anno, `4`  = 'C3')
CD8_subset_anno <- RenameIdents(CD8_subset_anno, `6`  = 'C4')


CD8_subset_anno@meta.data$seurat_clusters <- CD8_subset_anno@active.ident
DimPlot(CD8_subset_anno, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("CD8 T cells")
DimPlot(CD8_subset_anno, reduction = 'umap', raster=FALSE, label = F) + ggtitle("")+NoLegend()+NoAxes()

CD8_subset_anno

saveRDS(CD8_subset_anno, file="D:/DATA/Gastric cancer/GS_subsets/CD8_subset.rds")
CD8_subset_anno <- readRDS(file="D:/DATA/Gastric cancer/GS_subsets/CD8_subset.rds")


##################################### CD8 T cells clustering ######################################### 

do_NebulosaPlot(sample = CD8_subset_anno, features = "CD226", legend.position = "none", pt.size = 0.5, plot.title = "")
do_NebulosaPlot(sample = CD8_subset_anno, features = "EOMES", legend.position = "none", pt.size = 0.5, plot.title = "")
do_NebulosaPlot(sample = CD8_subset_anno, features = "ITGAE", legend.position = "none", pt.size = 0.5, plot.title = "", viridis_direction = -1)
do_NebulosaPlot(sample = CD8_subset_anno, features = "ISG15", legend.position = "none", pt.size = 0.5, plot.title = "")

dittoBarPlot(CD8_subset_anno, "seurat_clusters", group.by = "seurat_clusters", scale = "count",  retain.factor.levels = F, main = "")+ 
  scale_y_continuous(expand = c(0, 0)) + NoLegend()+ xlab(NULL)+ ylab(NULL)


genes <- c("CD226", "GZMB", "GNLY", "CD69", "PRF1",
           "EOMES", "TCF7", "B3GAT1", "CD28", "CX3CR1", "GZMK", "GZMM",
           "ITGAE", "ENTPD1", "PDCD1", "TIGIT", "HAVCR2", "LAG3", "CTLA4", "TNFRSF9", "TOX",
           "ISG15", "MX1", "IRF7", "IFNG", 
           "HLA-DRA", "CD40LG", "FAS", "IL7R", "TBX21")

# Dot plots
DotPlot(CD8_subset_anno, features = genes)+ RotatedAxis() + xlab(NULL) + ylab(NULL) + theme(axis.text.x = element_text(face = "italic"))
DotPlot(CD8_subset_anno, features = genes)+ RotatedAxis() + xlab(NULL) + ylab(NULL) +NoLegend() + 
theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))

do_ExpressionHeatmap(sample = CD8_subset_anno,
                     features = genes,
                     viridis_direction = -1)


##################################### Singiture ######################################### 



signatures <- list(
  Cytotoxic = c("GZMB", "PRF1"), 
  Exhaustion = c("PDCD1", "TIGIT", "HAVCR2", "LAG3", "CTLA4"),
  Polyfunction = c("IFNG", "TNFSF4"),
  IFNG = "IFNG",
  PDCD1 = "PDCD1", 
  TNFSF4 = "TNFSF4", 
  GZMB = "GZMB", 
  PRF1 = "PRF1"
)



CD8_subset_anno <- AddModuleScore_UCell(CD8_subset_anno,features=signatures, name=NULL)
CD8_subset_anno <- SmoothKNN(CD8_subset_anno, signature.names = names(signatures), reduction="pca")

FeaturePlot(CD8_subset_anno, reduction = "umap", features = "Cytotoxic_kNN")
VlnPlot(CD8_subset_anno, features = "Exhaustion_kNN", pt.size = 0, split.by = "seurat_clusters", sort = T)+ NoLegend()
VlnPlot(CD8_subset_anno, features = "Polyfunction_kNN", pt.size = 0, split.by = "seurat_clusters", sort = "increasing")+ NoLegend()
VlnPlot(CD8_subset_anno, features = "Cytotoxic_kNN", pt.size = 0, split.by = "seurat_clusters", sort = T)+ NoLegend()

VlnPlot(CD8_subset_anno, features = "PDCD1_kNN", pt.size = 0, split.by = "seurat_clusters", sort = "increasing")+ NoLegend()
VlnPlot(CD8_subset_anno, features = "IFNG_kNN", pt.size = 0, split.by = "seurat_clusters", sort = "increasing")+ NoLegend()
VlnPlot(CD8_subset_anno, features = "GZMB_kNN", pt.size = 0, split.by = "seurat_clusters", sort = "increasing")+ NoLegend()


VlnPlot(CD8_subset_anno, features = "Polyfunction", pt.size = 0, split.by = "seurat_clusters")+ NoLegend()
VlnPlot(CD8_subset_anno, features = "Cytotoxic", pt.size = 0, split.by = "seurat_clusters")+ NoLegend()


VlnPlot(CD8_subset_anno, features = "IFNG", pt.size = 0, split.by = "seurat_clusters")+ NoLegend()+ylim(-1, 10)
VlnPlot(CD8_subset_anno, features = "PDCD1", pt.size = 0, split.by = "seurat_clusters")+ NoLegend()+ylim(-1, 10)


##################################### Normal vs Tumors #########################################

DimPlot(CD8_subset_anno, reduction = 'umap', group.by = "Type", raster=FALSE, label = F, cols = c("#6495ED", "lightgrey"))+
  NoAxes()+NoLegend()+ ggtitle("")


dittoBarPlot(CD8_subset_anno, "seurat_clusters", group.by = "Type", scale = "count", main = "")+ 
  NoLegend()+ xlab(NULL)+ ylab(NULL)+
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))

dittoBarPlot(CD8_subset_anno, "seurat_clusters", group.by = "Type", scale = "percent", main = "")+ 
  NoLegend()+ xlab(NULL)+ ylab(NULL)+
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))


cancer_subset <- subset(x = CD8_subset_anno, Type == 'Tumor' & Patient_ID_SPC != FALSE)
SCpubr::do_BarPlot(sample = cancer_subset, 
                   group.by = "seurat_clusters", 
                   split.by = "Patient_ID_SPC",
                   position = "Fill",legend.position = "bottom",
                   plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)

Normal_subset <- subset(x = CD8_subset_anno, Type == 'Normal' & Patient_ID_SPC != FALSE)
SCpubr::do_BarPlot(sample = Normal_subset, 
                   group.by = "seurat_clusters", 
                   split.by = "Patient_ID_SPC",
                   position = "Fill",legend.position = "bottom",
                   plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)


CD8_subset_anno_notMETA <- subset(x = CD8_subset_anno, Stage == 'Stage II' | Stage == 'Stage III')
view(CD8_subset_anno_notMETA@meta.data)

SCpubr::do_BarPlot(sample = CD8_subset_anno_notMETA, 
                   group.by = "seurat_clusters", 
                   split.by = "Type",
                   position = "stack",legend.position = "bottom",
                   plot.title = "", xlab = "",  ylab = "Number of cells per cluster", font.size = 11)

SCpubr::do_BarPlot(sample = CD8_subset_anno_notMETA, 
                   group.by = "seurat_clusters", 
                   split.by = "Type",
                   position = "Fill",legend.position = "bottom",
                   plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)


cancer_subset_stage <- subset(x = CD8_subset_anno_notMETA, Type == 'Tumor' & Patient_ID_SPC != FALSE)
SCpubr::do_BarPlot(sample = cancer_subset_stage, 
                   group.by = "seurat_clusters", 
                   split.by = "Patient_ID_SPC",
                   position = "Fill",legend.position = "bottom",
                   plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)


dittoBarPlot(cancer_subset_stage, "seurat_clusters", group.by = "Patient_ID_SPC", scale = "percent", main = "")+ 
  NoLegend()+ xlab(NULL)+ ylab(NULL)+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))


Normal_subset_stage <- subset(x = CD8_subset_anno_notMETA, Type == 'Normal' & Patient_ID_SPC != FALSE)
SCpubr::do_BarPlot(sample = Normal_subset_stage, 
                   group.by = "seurat_clusters", 
                   split.by = "Patient_ID_SPC",
                   position = "Fill",legend.position = "bottom",
                   plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)

dittoBarPlot(Normal_subset_stage, "seurat_clusters", group.by = "Patient_ID_SPC", scale = "percent", main = "")+ 
  NoLegend()+ xlab(NULL)+ ylab(NULL)+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))

dittoBarPlot(Normal_subset_stage, "seurat_clusters", group.by = "Patient_ID_SPC", scale = "percent", main = "")





##################################### Stats PAIRED ★ ######################################### 
##################################### Dataframe for Outcome
md.immune_annotation<- CD8_subset_anno@meta.data %>% as.data.table

length(which(md.immune_annotation$Patient_ID==68 & md.immune_annotation$Type=='Tumor'& md.immune_annotation$seurat_clusters==1))

# Vector of sample numbers
sample_numbers <- c(15, 18, 22, 28, 29, 36, 41, 5, 54, 56, 59, 64, 65, 66, 68, 7, 71, 73, 8)

# Vector of cluster names
cluster_vec <- unique(md.immune_annotation$seurat_clusters)
type_vec <- c("Tumor", "Normal")

unique(cancer_subset@meta.data$Patient_ID_SPC)



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
    for (number in sample_numbers) {
      # Calculate the length of the subset
      subset_length <- length(which(md.immune_annotation$Patient_ID_SPC == number & md.immune_annotation$Type == type & md.immune_annotation$seurat_clusters == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.immune_annotation$Patient_ID_SPC == number & md.immune_annotation$Type == type))
      subset_fraction <- subset_length / total_count
      # Add a row to the data frame
      new_row <- data.frame(sample_numbers = number,
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

################################################################################## Graph

# Create an empty list to store the plots
plot_list <- list()

for (cluster in cluster_vec) {
  # Subset the data for each type_vec
  subset_df <- output_df %>% subset(cluster_vec == cluster)
  normal_df <- subset_df %>% subset(type_vec == "Normal")
  tumor_df <- subset_df %>% subset(type_vec == "Tumor")
  # Perform the unpaired t-test
  ttest <- t.test(normal_df$subset_fraction, tumor_df$subset_fraction)
  wtest <- wilcox.test(normal_df$subset_fraction, tumor_df$subset_fraction)
  
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction)) +
    geom_boxplot(color = "black", outlier.shape = NA, fill = c("#E1D4BB", "#537188"), alpha=0.8)+
    geom_point(shape = 21, color = "black", alpha=0.8, aes(fill=type_vec)) + 
    scale_fill_manual(values = c("grey", "brown")) +
    scale_y_continuous(labels = function(x) paste0(x * 100))+
    geom_line(aes(group = sample_numbers), alpha=0.5)+
    ggtitle("") + 
    xlab("") + ylab("") +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          panel.border = element_rect(color = "black", fill = NA, size = 0.2),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          axis.ticks.length = unit(0.1, "cm"),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          legend.position = "none")
  
  #stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
  #annotate("text", x = 1.5, y = 0.7, label = paste("p =", format(ttest$p.value, digits = 3)))+
  #annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[which(cluster_vec == cluster)]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 1)


############## Stat add
plot_list <- list()
for (cluster in cluster_vec) {
  # Subset the data for each type_vec
  subset_df <- output_df %>% subset(cluster_vec == cluster)
  normal_df <- subset_df %>% subset(type_vec == "Normal")
  tumor_df <- subset_df %>% subset(type_vec == "Tumor")
  # Perform the unpaired t-test
  ttest <- t.test(normal_df$subset_fraction, tumor_df$subset_fraction)
  wtest <- wilcox.test(normal_df$subset_fraction, tumor_df$subset_fraction)
  
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction)) +
    geom_boxplot(color = "black", outlier.shape = NA, fill = c("grey", "brown"), alpha=0.5)+
    geom_point(shape = 21, color = "black", alpha=0.8, aes(fill=type_vec)) + 
    scale_fill_manual(values = c("grey", "brown")) +
    geom_line(aes(group = sample_numbers), alpha=0.5)+
    ggtitle(cluster) + 
    xlab("") + ylab("") +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          axis.ticks.length = unit(0.1, "cm"),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          legend.position = "none")+
    
    #stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
    annotate("text", x = 1.5, y = 0.7, label = paste("p =", format(ttest$p.value, digits = 3)))+
    annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[which(cluster_vec == cluster)]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 1)




























































VlnPlot(CD8_subset_filtered_xTCR, features = "PDCD1", pt.size = 0, split.by = "seurat_clusters")+ NoLegend()




Exhaustion = c("PDCD1", "HAVCR2", "LAG3", "TIGIT")















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


do_NebulosaPlot(sample = CD8_subset_anno, features = "EOMES", legend.position = "none", pt.size = 0.5, plot.title = "")
do_NebulosaPlot(sample = CD8_subset_anno, features = "CD226", legend.position = "none", pt.size = 0.5, plot.title = "")
do_NebulosaPlot(sample = CD8_subset_anno, features = "ITGAE", legend.position = "none", pt.size = 0.5, plot.title = "")
do_NebulosaPlot(sample = CD8_subset_anno, features = "ISG15", legend.position = "none", pt.size = 0.5, plot.title = "")







do_NebulosaPlot(sample = CD8_subset_anno, features = c("EOMES", "TCF7"), joint = TRUE)



plot_density(CD8_subset_anno, c("CD28", "B3GAT1"), joint = TRUE)
plot_density(CD8_subset_anno, c("CD69", "HLA-DRA"), joint = TRUE)
plot_density(CD8_subset_anno, c("GZMA", "GZMB", "GZMK", "GZMM", "GZMH", "GNLY"), joint = TRUE)
plot_density(CD8_subset_anno, c("EOMES", "TCF7"), joint = TRUE)
plot_density(CD8_subset_anno, c("CD40LG", "FAS"), joint = TRUE)
plot_density(CD8_subset_anno, c("ENTPD1", "ITGAE"), joint = TRUE)
plot_density(CD8_subset_anno, c("CD226", "EOMES"), joint = TRUE)
plot_density(CD8_subset_anno, c("CX3CR1", "IL7R"), joint = TRUE)
plot_density(CD8_subset_anno, c("TOX", "CD74"), joint = TRUE)


VlnPlot(CD8_subset_anno, features = "CD69", pt.size = 0, split.by = "seurat_clusters")+ NoLegend()
VlnPlot(CD8_subset_anno, features = "ENTPD1", pt.size = 0, split.by = "seurat_clusters")+ NoLegend()
VlnPlot(CD8_subset_anno, features = "CD226", pt.size = 0, split.by = "seurat_clusters")+ NoLegend()

VlnPlot(CD8_subset_anno, features = "CX3CR1", pt.size = 0, split.by = "seurat_clusters")+ NoLegend()


























CD8_subset_filtered_xTCR <- RunPCA(CD8_subset_filtered_xTCR)
CD8_subset_filtered_xTCR <- FindNeighbors(CD8_subset_filtered_xTCR, reduction = "pca", dims = 1:20)
CD8_subset_filtered_xTCR <- FindClusters(CD8_subset_filtered_xTCR, resolution = 0.25)
CD8_subset_filtered_xTCR <- RunUMAP(CD8_subset_filtered_xTCR, reduction = "pca", dims = 1:20)

DimPlot(CD8_subset_filtered_xTCR, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("CD8 T cells")





VlnPlot(CD8_subset_filtered_xTCR, features = "TCF7", pt.size = 0, split.by = "seurat_clusters")+ NoLegend()
cluster0.markers <- FindMarkers(CD8_subset_filtered_xTCR, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster5.markers <- FindMarkers(CD8_subset_filtered_xTCR, ident.1 = 5, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)


cluster0.markers <- FindMarkers(CD8_subset_filtered_xTCR, ident.1 = 0, ident.2 = 1, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)







CD8_subset_filtered_xTCR

DimPlot(CD8_subset_filtered_xTCR, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("CD8 T cells") |
  FeaturePlot(CD8_subset_filtered_xTCR, reduction = "umap", features = "ISG15")

DimPlot(CD8_subset_filtered_xTCR, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("CD8 T cells") |
  FeaturePlot(CD8_subset_filtered_xTCR, reduction = "umap", features = "CXCR6")

DimPlot(CD8_subset_filtered_xTCR, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("CD8 T cells") |
  FeaturePlot(CD8_subset_filtered_xTCR, reduction = "umap", features = "GZMK")

DimPlot(CD8_subset_filtered_xTCR, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("CD8 T cells") |
  FeaturePlot(CD8_subset_filtered_xTCR, reduction = "umap", features = "TCF7")

DimPlot(CD8_subset_filtered_xTCR, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("CD8 T cells") |
  FeaturePlot(CD8_subset_filtered_xTCR, reduction = "umap", features = "CCR9")

DimPlot(CD8_subset_filtered_xTCR, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("CD8 T cells") |
  FeaturePlot(CD8_subset_filtered_xTCR, reduction = "umap", features = "CD27")

VlnPlot(CD8_subset_filtered_xTCR, features = "CCL5", pt.size = 0, split.by = "seurat_clusters")+ NoLegend()


DimPlot(CD8_subset_filtered_xTCR, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("CD8 T cells") |
  FeaturePlot(CD8_subset_filtered_xTCR, reduction = "umap", features = "KLF2")



DimPlot(CD8_subset_filtered, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("CD8 T cells") |
FeaturePlot(CD8_subset_filtered, reduction = "umap", features = "ISG15")

DimPlot(CD8_subset_filtered, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("CD8 T cells") |
  FeaturePlot(CD8_subset_filtered, reduction = "umap", features = "GNLY")

DimPlot(CD8_subset_filtered, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("CD8 T cells") |
  FeaturePlot(CD8_subset_filtered, reduction = "umap", features = "GZMK")

DimPlot(CD8_subset_filtered, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("CD8 T cells") |
  FeaturePlot(CD8_subset_filtered, reduction = "umap", features = "GZMH")

DimPlot(CD8_subset_filtered, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("CD8 T cells") |
  FeaturePlot(CD8_subset_filtered, reduction = "umap", features = "GZMA")

DimPlot(CD8_subset_filtered, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("CD8 T cells") |
  FeaturePlot(CD8_subset_filtered, reduction = "umap", features = "GZMB")

DimPlot(CD8_subset_filtered, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("CD8 T cells") |
  FeaturePlot(CD8_subset_filtered, reduction = "umap", features = "CCR9")

DimPlot(CD8_subset_filtered, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("CD8 T cells") |
  FeaturePlot(CD8_subset_filtered, reduction = "umap", features = "EOMES")

DimPlot(CD8_subset_filtered, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("CD8 T cells") |
  FeaturePlot(CD8_subset_filtered, reduction = "umap", features = "TBX21")


DimPlot(CD8_subset_filtered, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("CD8 T cells") |
  FeaturePlot(CD8_subset_filtered, reduction = "umap", features = "TRA")

DimPlot(CD8_subset_filtered, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("CD8 T cells") |
  FeaturePlot(CD8_subset_filtered, reduction = "umap", features = "adt_CD4")


CD8_subset_anno <- CD8_subset_filtered

CD8_subset_anno <- RenameIdents(CD8_subset_anno, `0`  = 'C0')
CD8_subset_anno <- RenameIdents(CD8_subset_anno, `1`  = 'C1')
CD8_subset_anno <- RenameIdents(CD8_subset_anno, `3`  = 'C2')
CD8_subset_anno <- RenameIdents(CD8_subset_anno, `4`  = 'C3')
CD8_subset_anno <- RenameIdents(CD8_subset_anno, `5`  = 'C4')
CD8_subset_anno <- RenameIdents(CD8_subset_anno, `6`  = 'C5')
CD8_subset_anno <- RenameIdents(CD8_subset_anno, `7`  = 'C6')

CD8_subset_anno@meta.data$seurat_clusters <- CD8_subset_anno@active.ident

saveRDS(CD8_subset_anno, file="D:/DATA/Gastric cancer/GS_subsets/CD8_subset.rds")
CD8_subset_anno <- readRDS(file="D:/DATA/Gastric cancer/GS_subsets/CD8_subset.rds")

#######################################################################################################


FeaturePlot(CD8_subset, reduction = "umap", features = "ISG15")









# Markers
CD8_markers <- c("CCR7", "TCF7", "SELL", 
                 "ISG15", "IFIT3", "STAT1",
                 "GNLY", "PRF1",
                 "ENTPD1", "ITGAE",
                 "PDCD1", "TIGIT", "CTLA4",
                 "GZMA", "EOMES", "TBX21", "GZMK", "GZMB", "CX3CR1", "GZMH", "CD8A", "NKG7", "CCR9", "CCL5", "CD27", "KLF2", "CXCR5", "adt_CD8", "adt_CD4"
)

FeaturePlot(CD8_subset, reduction = "umap", features = "TRAV1-2")
FeaturePlot(CD8_subset, reduction = "umap", features = "ISG15")

# Dot plots

DotPlot(pos_cells, features = CD8_markers)+ RotatedAxis()

DotPlot(pos_cells, assay = "integrated", features = CD8_markers, col.min = -2,  col.max = 2)+ RotatedAxis() + xlab(NULL) + ylab(NULL) + 
  theme(axis.text.x = element_text(face = "italic"))

CD8_subset


CD8_subset
g$layers

g <- DotPlot(object = CD8_subset, features = CD8_markers, assay="integrated")
g <- g + 
  geom_point(mapping = aes_string(size = 'pct.exp', color = 'avg.exp')) +
  guides(color = guide_colorbar(title = 'Average Expression')) + 
  theme(axis.text.x = element_text(angle=90))

plot(g)










VlnPlot(CD8_subset_anno, features = "CD8A", pt.size = 0, split.by = "seurat_clusters")+ NoLegend()

CD8_subset_anno <- RunPCA(CD8_subset_anno)
CD8_subset_anno <- FindNeighbors(CD8_subset_anno, reduction = "pca", dims = 1:20)
CD8_subset_anno <- FindClusters(CD8_subset_anno, resolution = 0.3)
CD8_subset_anno <- RunUMAP(CD8_subset_anno, reduction = "pca", dims = 1:20)

DimPlot(CD8_subset_anno, reduction = 'umap', raster=FALSE, label = T) + ggtitle("")
VlnPlot(CD8_subset_anno, features = "CD8A", pt.size = 0, split.by = "seurat_clusters")+ NoLegend()

DimPlot(CD8_subset_anno, reduction = 'umap', raster=FALSE, label = F) + ggtitle("")+ NoLegend() +NoAxes()


CD8_subset_anno@meta.data$seurat_clusters <- CD8_subset_anno@active.ident


##############################################################################
CD8_subset_anno
saveRDS(CD8_subset_anno, file="D:/DATA/Gastric cancer/GS_subsets/CD8_subset.rds")
CD8_subset_anno <- readRDS(file="D:/DATA/Gastric cancer/GS_subsets/CD8_subset.rds")












CD8_subset_anno <- RenameIdents(CD8_subset_anno, `0`  = 'C0')
CD8_subset_anno <- RenameIdents(CD8_subset_anno, `1`  = 'C1')
CD8_subset_anno <- RenameIdents(CD8_subset_anno, `3`  = 'C2')
CD8_subset_anno <- RenameIdents(CD8_subset_anno, `4`  = 'C3')
CD8_subset_anno <- RenameIdents(CD8_subset_anno, `5`  = 'C4')
CD8_subset_anno <- RenameIdents(CD8_subset_anno, `6`  = 'C5')
CD8_subset_anno <- RenameIdents(CD8_subset_anno, `7`  = 'C6')

view(CD8_subset_anno@meta.data)
CD8_subset_anno

CD8_subset_anno@meta.data$seurat_clusters <- CD8_subset_anno@active.ident

saveRDS(CD8_subset_anno, file="D:/DATA/Gastric cancer/GS_subsets/CD8_subset.rds")
CD8_subset_anno <- readRDS(file="D:/DATA/Gastric cancer/GS_subsets/CD8_subset.rds")

##################################### CD8 T cells clustering ######################################### NEW
CD8_subset_filtered_Trex <- RunPCA(CD8_subset_filtered_Trex)
CD8_subset_filtered_Trex <- FindNeighbors(CD8_subset_filtered_Trex, reduction = "pca", dims = 1:50)
CD8_subset_filtered_Trex <- FindClusters(CD8_subset_filtered_Trex, resolution = 0.3)
CD8_subset_filtered_Trex <- RunUMAP(CD8_subset_filtered_Trex, reduction = "pca", dims = 1:50)

DimPlot(CD8_subset_filtered_Trex, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("CD8 T cells") 



############################################################################## Annotation #####

CD8_subset_anno <- CD8_subset

# Markers
CD8_markers <- c("CCR7", "TCF7", "SELL", 
                 "ISG15", "IFIT3", "STAT1",
                 "GNLY", "PRF1",
                 "ENTPD1", "ITGAE",
                 "PDCD1", "TIGIT", "CTLA4",
                 "GZMA", "EOMES", "TBX21", "GZMK", "GZMB", "CX3CR1", "GZMH", "CD8A", "NKG7", "CCR9", "CCL5", "CD27", "KLF2", "CXCR5", "adt_CD8", "adt_CD4"
                )

FeaturePlot(CD8_subset_anno, reduction = "umap", features = "adt_CD4")

# Dot plots
DotPlot(CD8_subset, features = CD8_markers, col.min = -2,  col.max = 2)+ RotatedAxis() + xlab(NULL) + ylab(NULL) + 
  theme(axis.text.x = element_text(face = "italic"))+
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")

SCpubr::do_ExpressionHeatmap(sample = CD8_subset_anno,
                             features = CD8_markers,
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

CD8_subset_anno <- AddModuleScore_UCell(CD8_subset_anno,features=signatures, name=NULL)
CD8_subset_anno <- SmoothKNN(CD8_subset_anno, signature.names = names(signatures), reduction="pca")

FeaturePlot(CD8_subset_anno, reduction = "umap", features = "Cytotoxic_kNN")
VlnPlot(CD8_subset_anno, features = "Cytotoxic_kNN", pt.size = 0, split.by = "seurat_clusters")+ NoLegend()

FeaturePlot(CD8_subset_anno, reduction = "umap", features = "Exhaustion_kNN")
VlnPlot(CD8_subset_anno, features = "Exhaustion_kNN", pt.size = 0, split.by = "seurat_clusters")+ NoLegend()


CD8.markers <- FindAllMarkers(CD8_subset_anno, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CD8.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(CD8_subset_anno, features = top10$gene) + NoLegend()


##################################### Stats #########################################

CD8_subset_cancer <- subset(x = CD8_subset_anno, Type == 'Tumor' & Outcome != 'Meta' & Patient_ID_SPC != F)
CD8_subset_cancer_MSI <- subset(x = CD8_subset_cancer, MSI_type == 'MSI')
CD8_subset_cancer_MSS <- subset(x = CD8_subset_cancer, MSI_type == 'MSS')

SCpubr::do_BarPlot(sample = CD8_subset_cancer, 
                   group.by = "seurat_clusters", 
                   legend.position = "none", 
                   plot.title = "",  xlab = "",  ylab = "Number of cells per cluster", font.size = 11)


# 전체
SCpubr::do_BarPlot(sample = CD8_subset_cancer, 
                   group.by = "seurat_clusters", 
                   split.by = "Patient_ID_SPC",
                   position = "Fill",legend.position = "bottom",
                   plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)

# 부분
SCpubr::do_BarPlot(sample = CD8_subset_cancer_MSI, 
                   group.by = "seurat_clusters", 
                   split.by = "Patient_ID_SPC",
                   position = "Fill",legend.position = "bottom",
                   plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)

SCpubr::do_BarPlot(sample = CD8_subset_cancer_MSS, 
                   group.by = "seurat_clusters", 
                   split.by = "Patient_ID_SPC",
                   position = "Fill",legend.position = "bottom",
                   plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)


##################################### Dataframe for Outcome






















################################################################################## Annotation #####
CD8_subset

DimPlot(CD8_subset_outcome, reduction = 'umap', raster=FALSE, label = F) + NoLegend()+NoAxes()
DimPlot(CD8_subset_outcome, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("CD8 T cells GS")
DimPlot(CD8_subset_outcome, reduction = 'umap', split.by = "Outcome", raster=FALSE, label = T) + NoLegend() + ggtitle("CD8 T cells GS")


SCpubr::do_BarPlot(sample = CD8_subset, 
                   group.by = "seurat_clusters", 
                   position = "stack",legend.position = "none",
                   plot.title = "", xlab = "",  ylab = "Number of cells per cluster", font.size = 11)

################################################################################## Save data ----
