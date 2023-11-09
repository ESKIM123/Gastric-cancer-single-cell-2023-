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


##################################### Recur vs Non-recur : in CD8 T cells#########################################
CD8_subset <- subset(x = immune_annotation, clusters == 'CD8 T cells')

CD8_subset <- RunPCA(CD8_subset)
CD8_subset <- FindNeighbors(CD8_subset, reduction = "pca", dims = 1:20)
CD8_subset <- FindClusters(CD8_subset, resolution = 0.25)
CD8_subset <- RunUMAP(CD8_subset, reduction = "pca", dims = 1:20)
# 1:20 - 0.25

DimPlot(CD8_subset, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("CD8 T cells")
DimPlot(CD8_subset, reduction = 'umap', raster=FALSE, label = F) + ggtitle("")+ NoLegend()+NoAxes()
DimPlot(CD8_subset, reduction = 'umap', split.by = "Outcome", raster=FALSE, label = T) + NoLegend() + ggtitle("CD8 T cells")


##################################### Stats #########################################
CD8_subset_outcome <- subset(x = CD8_subset, Type == 'Tumor' & Outcome != 'Meta')
immune_annotation_outcome <- subset(x = immune_annotation, Type == 'Tumor' & Outcome != 'Meta')

##################################### Dataframe for Outcome
md.CD8_subset_outcome<- CD8_subset_outcome@meta.data %>% as.data.table

# Vector of sample numbers
sample_numbers <- c(4, 5, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 51, 52, 53, 54, 56, 59, 60, 62, 64, 65, 66, 68, 71, 73, 75)

# Vector of cluster names
cluster_vec <- unique(md.CD8_subset_outcome$seurat_clusters)
type_vec <- c("Non_recur", "Recur")

unique(md.CD8_subset_outcome$Patient_ID_SPC)

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
      subset_length <- length(which(md.CD8_subset_outcome$Patient_ID_SPC == number & md.CD8_subset_outcome$Outcome == type & md.CD8_subset_outcome$seurat_clusters == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.CD8_subset_outcome$Patient_ID_SPC == number & md.CD8_subset_outcome$Outcome == type))
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

################################################################################## Graph (in CD8 T cells) #### 

# Create an empty list to store the plots
plot_list <- list()

for (cluster in cluster_vec) {
  # Subset the data for each type_vec
  subset_df <- output_df %>% subset(cluster_vec == cluster)
  Non_recur_df <- subset_df %>% subset(type_vec == "Non_recur")
  Recur_df <- subset_df %>% subset(type_vec == "Recur")
  # Perform the unpaired t-test
  ttest <- t.test(Non_recur_df$subset_fraction, Recur_df$subset_fraction)
  wtest <- wilcox.test(Non_recur_df$subset_fraction, Recur_df$subset_fraction)
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
    geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) +
    geom_boxplot(alpha = 0.5, color = "black", outlier.shape = NA)+
    scale_fill_manual(values = c("brown", "grey")) +
    labs(title = paste(cluster), #paste(cluster)
         x = NULL, y = NULL) +
    #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          panel.border = element_rect(color = "black", fill = NA, size = 0.2),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          axis.text.x = element_text(angle = 45, hjust = 1),
          #axis.ticks.length = unit(0.1, "cm"),
          #plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          plot.title = element_text(size = 8, hjust = 0.5), 
          legend.position = "none")+
  
  #stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")
  annotate("text", x = 1.5, y = 0.5, label = paste("tp =", format(ttest$p.value, digits = 3)))+
  annotate("text", x = 1.5, y = 0.2, label = paste("wp =", format(wtest$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[which(cluster_vec == cluster)]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 1) |
  grid.text("CD8 T cells", y = 0.05)

################################################################################## Graph -- clean

# Create an empty list to store the plots
plot_list <- list()

for (cluster in cluster_vec) {
  # Subset the data for each type_vec
  subset_df <- output_df %>% subset(cluster_vec == cluster)
  Non_recur_df <- subset_df %>% subset(type_vec == "Non_recur")
  Recur_df <- subset_df %>% subset(type_vec == "Recur")
  # Perform the unpaired t-test
  ttest <- t.test(Non_recur_df$subset_fraction, Recur_df$subset_fraction)
  wtest <- wilcox.test(Non_recur_df$subset_fraction, Recur_df$subset_fraction)
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
    geom_boxplot(color = "black", outlier.shape = NA)+
    scale_fill_manual(values = c("brown", "grey")) +
    labs(title = paste(cluster), #paste(cluster)
         x = NULL, y = NULL) +
    #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          panel.border = element_rect(color = "black", fill = NA, size = 0.2),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          axis.text.x = element_text(angle = 45, hjust = 1),
          #axis.ticks.length = unit(0.1, "cm"),
          #plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          plot.title = element_text(size = 12, hjust = 0.5), 
          legend.position = "none")
    
    #stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")
    #annotate("text", x = 1.5, y = 0.5, label = paste("p =", format(ttest$p.value, digits = 3)))
    #annotate("text", x = 1.5, y = 0.8, label = paste("p =", format(wtest$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[which(cluster_vec == cluster)]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 1)



################################################################################## Graph (in Total cells) ####
md.immune_annotation_outcome<- immune_annotation_outcome@meta.data %>% as.data.table

# Vector of cluster names
cluster_vec <- unique(md.CD8_subset_outcome$seurat_clusters)
type_vec <- c("Non_recur", "Recur")

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
      subset_length <- length(which(md.CD8_subset_outcome$Patient_ID_SPC == number & md.CD8_subset_outcome$Outcome == type & md.CD8_subset_outcome$seurat_clusters == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.immune_annotation_outcome$Patient_ID_SPC == number & md.immune_annotation_outcome$Outcome == type))
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


# Create an empty list to store the plots
plot_list <- list()

for (cluster in cluster_vec) {
  # Subset the data for each type_vec
  subset_df <- output_df %>% subset(cluster_vec == cluster)
  Non_recur_df <- subset_df %>% subset(type_vec == "Non_recur")
  Recur_df <- subset_df %>% subset(type_vec == "Recur")
  # Perform the unpaired t-test
  ttest <- t.test(Non_recur_df$subset_fraction, Recur_df$subset_fraction)
  wtest <- wilcox.test(Non_recur_df$subset_fraction, Recur_df$subset_fraction)
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
    geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) +
    geom_boxplot(alpha = 0.5, color = "black", outlier.shape = NA)+
    scale_fill_manual(values = c("brown", "grey")) +
    labs(title = paste(cluster), #paste(cluster)
         x = NULL, y = NULL) +
    #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          panel.border = element_rect(color = "black", fill = NA, size = 0.2),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          axis.text.x = element_text(angle = 45, hjust = 1),
          #axis.ticks.length = unit(0.1, "cm"),
          #plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          plot.title = element_text(size = 8, hjust = 0.5), 
          legend.position = "none")+
    
    #stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")
    annotate("text", x = 1.5, y = 0.5, label = paste("tp =", format(ttest$p.value, digits = 3)))+
    annotate("text", x = 1.5, y = 0.2, label = paste("wp =", format(wtest$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[which(cluster_vec == cluster)]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 1) |
  grid.text("CD8 T cells (Total immune)", y = 0.05)

### clean
# Create an empty list to store the plots
plot_list <- list()

for (cluster in cluster_vec) {
  # Subset the data for each type_vec
  subset_df <- output_df %>% subset(cluster_vec == cluster)
  Non_recur_df <- subset_df %>% subset(type_vec == "Non_recur")
  Recur_df <- subset_df %>% subset(type_vec == "Recur")
  # Perform the unpaired t-test
  ttest <- t.test(Non_recur_df$subset_fraction, Recur_df$subset_fraction)
  wtest <- wilcox.test(Non_recur_df$subset_fraction, Recur_df$subset_fraction)
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
    geom_boxplot(color = "black", outlier.shape = NA)+
    scale_fill_manual(values = c("brown", "grey")) +
    labs(title = paste(cluster), #paste(cluster)
         x = NULL, y = NULL) +
    #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          panel.border = element_rect(color = "black", fill = NA, size = 0.2),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          axis.text.x = element_text(angle = 45, hjust = 1),
          #axis.ticks.length = unit(0.1, "cm"),
          #plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          plot.title = element_text(size = 12, hjust = 0.5), 
          legend.position = "none")
  
  #stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")
  #annotate("text", x = 1.5, y = 0.5, label = paste("p =", format(ttest$p.value, digits = 3)))
  #annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[which(cluster_vec == cluster)]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 1)

################################################################################## Cell count #####
CD8_subset

DimPlot(CD8_subset, reduction = 'umap', raster=FALSE, label = F) + NoLegend()+NoAxes()

SCpubr::do_BarPlot(sample = CD8_subset, 
                   group.by = "seurat_clusters", 
                   position = "stack",legend.position = "none",
                   plot.title = "", xlab = "",  ylab = "Number of cells per cluster", font.size = 11)

################################################################################## Save data ----

saveRDS(CD8_subset, file="D:/DATA/Gastric cancer/CD8_subset.rds")
CD8_subset <- readRDS(file="D:/DATA/Gastric cancer/CD8_subset.rds")

################################################################################## Annotation #####

CD8_anno <- CD8_subset

CD8_anno.markers <- FindAllMarkers(CD8_anno, only.pos = TRUE, logfc.threshold = 0.25, test.use = 'roc')
write.csv(CD8_anno.markers, file="D:/DATA/Gastric cancer/CD8_anno_markers.csv", row.names = FALSE)


signatures <- list(
  CD8_naive_like = c("LEF1", "TCF7", "KLF2", "TXK", "BACH2", "LTB", "FLT3LG", "TNFSF8", "CMTM8", "IL23A", "TIMP1", "WNT7A", "CCR7", "IL7R", "IL6R",
                     "IFNGR2", "SELL", "MAL", "EEF1A1", "ACTN1", "TRABD2A", "TPT1", "EEF1B2", "NELL2", "NOSIP", "PABPC1"),
  CD8_memory_like = c("SELL-",  "CCR7-",  "TBX21",  "EOMES",  "IL7R", "KLRG1", "CD27"),
  CD8_effector_like = c("TBX21",  "GZMB",  "PRDM1",  "PRF1", "IFNG"), 
  terminal_Tex_T_cells = c("RBPJ", "ETV1", "TOX", "ZBED2", "TOX2", "CXCL13", "TNFSF4", "FAM3C", "GZMB", "CSF1", "CCL3", "CD70", "IFNG", 
                           "NAMPT", "FASLG",	"IL2RA", "CXCR6", "CD74", "IL2RB", "IL2RG", "TNFRSF9", "LAYN", "ENTPD1", "HAVCR2", "CTLA4",
                           "KRT86", "TNFRSF18", "GEM", "TIGIT", "DUSP4"),
  NK_like_T_cells = c("IKZF2", "EOMES", "NR4A2", "LITAF", "ZNF331", "GZMK", "XCL2", "XCL1", "GZMM", "TNFSF9", "IFNGR1", "CXCR4", "CD74", 
                      "IL2RB", "KIR2DL3", "KIR3DL2", "CD160", "TYROBP", "KLRF1", "KLRD1", "KLRG1", "CMC1", "GCSAM", "DUSP2"),
  ISG_CD8_cells = c("PLSCR1", "STAT1", "IRF7", "STAT2", "SP100", "TNFSF10", "CCR1", "CD74", "IFIT1", "RSAD2", "IFIT3", "IFI44L", 
                      "MX1", "IFI6", "OAS1", "CMPK2", "ISG15", "OAS3"),
  GZMK_TEM_CD8_cells = c( "EOMES", "GZMK", "CXCR4", "CD74", "CXCR5", "CCR4", "CD44", "DUSP2", "CMC1", "CST7", "DKK3", "SH2D1A", "ENC1", "TRAT1", 
                          "DTHD1", "PIK3R1", "CRTAM", "SUB1", "GZMA", "CCL5", "GZMH", "IL32", "CCL4", "CCR5", "CXCR3", "HLA-DRB1", "HLA-DPA1", "COTL1",
                          "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRB5", "ITM2C", "APOBEC3G"),
  TEMRA_CD8_cells = c("ASCL2", "KLF2", "KLF3", "ZEB2", "TBX21", "GZMH", "GZMB", "GZMM", "GZMA", "CX3CR1", "CXCR2", "CMKLR1", "CXCR1", "FGFBP2", "FCGR3A",
                      "S1PR5", "PRSS23", "GNLY", "NKG7", "KLRD1", "FGR", "PLEK", "C1orf21"),
  IL7R_TM_CD8_cells = c("ZFP36L2", "TSC22D3", "IL7R", "CXCR4", "ZFP36", "BTG1", "ANXA1", "LMNA", "CD55", "TNFAIP3", "FTH1", "RGCC", "GPR183", "PABPC1"),
  ZNF683_CXCR6_TM = c("ZNF683", "CCL5", "IL32", "GZMA", "TNF", "CKLF", "ACTB", "CD52", "TRAF3IP3", "SH3BGRL3", "SIT1", "S100A4", "CISH", "MYL12A", "PTPRCAP", "PFN1"),
  TM = c("ZNF683", "CCL5", "IL32", "GZMA", "TNF", "CKLF", "ACTB", "CD52", "TRAF3IP3", "SH3BGRL3", "SIT1", "S100A4", "CISH", "MYL12A", "PTPRCAP", "PFN1",
         "ZFP36L2", "TSC22D3", "IL7R", "CXCR4", "ZFP36", "BTG1", "ANXA1", "LMNA", "CD55", "TNFAIP3", "FTH1", "RGCC", "GPR183", "PABPC1"),
  TCF7_Tex = c("TSHZ2","TCF7","NR3C1","TOX","BATF","CXCL13","EBI3","TNFSF8","CD40LG","CCR7","IL6R","IFNAR2","CCR4",
               "LHFP","CD200","GNG4","TNFRSF4","IGFL2","CPM","NMB","SESN3","BTLA","IGFBP4"),
  Tc17 = c("RORC","ZBTB16","CEBPD","RORA","NR1D1", "CCL20","CD40LG","LTB","TNFSF13B","IL26","IL17A","FLT3LG","TNF","IL23A", "CCR6","IL23R","IL7R","IL17RE","IL18RAP",
           "SLC4A10","KLRB1","TMIGD2","IL4I1","NCR3","LTK","CA2","ME1","AQP3","ADAM12"),
  TCM = c("SELL", "CCR7", "IL7R", "BCL2", "TBX21", "CD27", "CD28", "LEF1", "EOMES", "KLRG1"),
  ZNF683_CXCR6_Trm = c("ZNF683", "HOPX", "ID2", "ZFP36L2", "RBPJ", "CKLF", "IL32", "GZMB", "XCL1", "CCL5", "GZMA", "XCL2", "CXCR6", "CXCR3", "CAPG", 
                       "TMSB4X", "S100A4", "LGALS3", "ACTB", "SH3BGRL3", "CD52", "LGALS1", "ITGA1", "LDLRAD4"),
  Tcm_CD8_cells = c("IFNG", "IL2", "TNF","CCR7", "CD27", "CD28", "PTPRC", "SELL", "IL7R", "CCR7","EOMES", "TBX21")
  )


CD8_anno <- AddModuleScore_UCell(CD8_anno,features=signatures, name=NULL)
CD8_anno <- SmoothKNN(CD8_anno, signature.names = names(signatures), reduction="pca")


FeaturePlot(CD8_anno, reduction = "umap", features = c("CD8_naive_like","CD8_naive_like_kNN"))
VlnPlot(CD8_anno, features = "CD8_naive_like_kNN", pt.size = 0, split.by = "seurat_clusters")

FeaturePlot(CD8_anno, reduction = "umap", features = c("ISG_CD8_cells","ISG_CD8_cells_kNN"))
VlnPlot(CD8_anno, features = "ISG_CD8_cells_kNN", pt.size = 0, split.by = "seurat_clusters")

FeaturePlot(CD8_anno, reduction = "umap", features = c("GZMK_TEM_CD8_cells","GZMK_TEM_CD8_cells_kNN"))
VlnPlot(CD8_anno, features = "GZMK_TEM_CD8_cells_kNN", pt.size = 0, split.by = "seurat_clusters")

FeaturePlot(CD8_anno, reduction = "umap", features = c("TM","TM_kNN"))
VlnPlot(CD8_anno, features = "TM_kNN", pt.size = 0, split.by = "seurat_clusters")

FeaturePlot(CD8_anno, reduction = "umap", features = c("Tcm_CD8_cells","Tcm_CD8_cells_kNN"))
VlnPlot(CD8_anno, features = "Tcm_CD8_cells_kNN", pt.size = 0, split.by = "seurat_clusters")

FeaturePlot(CD8_anno, reduction = "umap", features = c("TEMRA_CD8_cells","TEMRA_CD8_cells_kNN"))
VlnPlot(CD8_anno, features = "TEMRA_CD8_cells_kNN", pt.size = 0, split.by = "seurat_clusters")

FeaturePlot(CD8_anno, reduction = "umap", features = c("IL7R_TM_CD8_cells","IL7R_TM_CD8_cells_kNN"))
VlnPlot(CD8_anno, features = "IL7R_TM_CD8_cells_kNN", pt.size = 0, split.by = "seurat_clusters")

FeaturePlot(CD8_anno, reduction = "umap", features = c("ZNF683_CXCR6_TM","ZNF683_CXCR6_TM_kNN"))
VlnPlot(CD8_anno, features = "ZNF683_CXCR6_TM_kNN", pt.size = 0, split.by = "seurat_clusters")

FeaturePlot(CD8_anno, reduction = "umap", features = c("TCF7_Tex","TCF7_Tex_kNN"))
VlnPlot(CD8_anno, features = "TCF7_Tex_kNN", pt.size = 0, split.by = "seurat_clusters")

FeaturePlot(CD8_anno, reduction = "umap", features = c("Tc17","Tc17_kNN"))
VlnPlot(CD8_anno, features = "Tc17_kNN", pt.size = 0, split.by = "seurat_clusters")

FeaturePlot(CD8_anno, reduction = "umap", features = c("TCM","TCM_kNN"))
VlnPlot(CD8_anno, features = "TCM_kNN", pt.size = 0, split.by = "seurat_clusters")

FeaturePlot(CD8_anno, reduction = "umap", features = c("ZNF683_CXCR6_Trm","ZNF683_CXCR6_Trm_kNN"))
VlnPlot(CD8_anno, features = "ZNF683_CXCR6_Trm_kNN", pt.size = 0, split.by = "seurat_clusters")

##################################################################################
i = "CENPF"
DefaultAssay(CD8_anno) <- "integrated"
FeaturePlot(CD8_anno, features = i, min.cutoff = '0', max.cutoff = '3', raster=FALSE) |
  dittoBoxPlot(CD8_anno, assay = "integrated", i, group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
               boxplot.width = 0.7, max=1) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))
##################################################################################
cluster2.markers <- FindMarkers(Plasma_anno, ident.1 = 2, ident.2 = c(0,1), only.pos = F, test.use = 'roc')
write.csv(cluster2.markers, file="D:/DATA/Gastric cancer/CD8_anno_markers_cluster2.csv", row.names = T)


MALAT1, NEAT1, DNAJA4, AHNAK, KLF6, CCL24

##################################################################################

dittoBoxPlot(CD8_anno, assay = "integrated", "PRF1", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
             boxplot.width = 0.7, max=5) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))

dittoBoxPlot(CD8_anno, assay = "integrated", "GNLY", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
             boxplot.width = 0.7, max=5) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))

dittoBoxPlot(CD8_anno, assay = "integrated", "GZMH", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
             boxplot.width = 0.7, max=5) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))

dittoBoxPlot(CD8_anno, assay = "integrated", "PDCD1", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
             boxplot.width = 0.7, max=5) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))



dittoBoxPlot(CD8_anno, assay = "integrated", "B3GAT1", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F, adjustment = "z-score", 
             boxplot.width = 0.7, min = -1, max = 1) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))

dittoBoxPlot(CD8_anno, assay = "integrated", "HNRNPLL", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
             boxplot.width = 0.7, min=-1, max=5) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))


dittoBoxPlot(CD8_anno, assay = "ADT", "CD45RA", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F, adjustment = "relative.to.max",
             boxplot.width = 0.7) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))

dittoBoxPlot(CD8_anno, assay = "ADT", "CD69", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F, adjustment = "relative.to.max",
             boxplot.width = 0.7) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))


FeaturePlot(CD8_anno, features = 'CD69', min.cutoff = '0', max.cutoff = '3', raster=FALSE)






##################################################################################
b.seu <- CD8_anno
DefaultAssay(b.seu) <- "RNA"
b.seu <- NormalizeData(b.seu)
b.seu <- FindVariableFeatures(b.seu)
b.seu <- ScaleData(b.seu)
b.seu <- RunPCA(b.seu)
b.seu <- FindNeighbors(b.seu, dims = 1:30)
b.seu <- FindClusters(b.seu, resolution = 0.9)
b.seu <- RunUMAP(b.seu, dims = 1:30, n.neighbors = 50)
DefaultAssay(b.seu) <- "integrated"

dittoBoxPlot(b.seu, assay = "RNA", "HNRNPLL", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F, adjustment = "relative.to.max",
             boxplot.width = 0.7) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))


DimPlot(b.seu, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Initial")


FeatureScatter(cbmc, feature1 = "adt_CD19", feature2 = "adt_CD3")






CD8_anno_0 <- subset(x = CD8_anno, seurat_clusters == 0)
CD8_anno_1 <- subset(x = CD8_anno, seurat_clusters == 1)
CD8_anno_2 <- subset(x = CD8_anno, seurat_clusters == 2)
CD8_anno_3 <- subset(x = CD8_anno, seurat_clusters == 3)
CD8_anno_4 <- subset(x = CD8_anno, seurat_clusters == 4)
CD8_anno_5 <- subset(x = CD8_anno, seurat_clusters == 5)
CD8_anno_6 <- subset(x = CD8_anno, seurat_clusters == 6)
CD8_anno_7 <- subset(x = CD8_anno, seurat_clusters == 7)

FeatureScatter(CD8_anno_0, feature1 = "adt_CD4", feature2 = "adt_CD8")
FeatureScatter(CD8_anno_3, feature1 = "adt_CD69", feature2 = "adt_CD103")
FeatureScatter(CD8_anno_0, feature1 = "adt_CD69", feature2 = "adt_CD103")
FeatureScatter(CD8_anno_7, feature1 = "adt_CD69", feature2 = "adt_CD103")
FeatureScatter(CD8_anno_1, feature1 = "adt_CD69", feature2 = "adt_CD103")
FeatureScatter(CD8_anno_2, feature1 = "adt_CD69", feature2 = "adt_CD103")


FeatureScatter(CD8_anno_0, feature1 = "adt_CD197", feature2 = "adt_CD45RA")
FeatureScatter(CD8_anno_1, feature1 = "adt_CD197", feature2 = "adt_CD45RA")
FeatureScatter(CD8_anno_2, feature1 = "adt_CD197", feature2 = "adt_CD45RA")
FeatureScatter(CD8_anno_3, feature1 = "adt_CD197", feature2 = "adt_CD45RA")
FeatureScatter(CD8_anno_4, feature1 = "adt_CD197", feature2 = "adt_CD45RA")
FeatureScatter(CD8_anno_5, feature1 = "adt_CD197", feature2 = "adt_CD45RA")
FeatureScatter(CD8_anno_6, feature1 = "adt_CD197", feature2 = "adt_CD45RA")
FeatureScatter(CD8_anno_7, feature1 = "adt_CD197", feature2 = "adt_CD45RA")




pl1 <- FeatureScatter(CD8_anno_2, feature1 = "adt_CD197", feature2 = "adt_CD45RA")











CD62L^lo, CCR7^lo T-bet Eomes CD127, KLRG1, and CD27

VlnPlot(CD8_anno, features = "SVIL", pt.size = 0, split.by = "seurat_clusters")


genes <- c("CD2","CSF1R")
seurat.object <- SmoothKNN(seurat.object, signature.names=genes,
                           assay="RNA", reduction="pca", k=20, suffix = "_smooth")


DefaultAssay(seurat.object) <- "RNA"
a <- FeaturePlot(seurat.object, reduction = "umap", features = genes)
DefaultAssay(seurat.object) <- "RNA_smooth"
b <- FeaturePlot(seurat.object, reduction = "umap", features = genes)
a / b











signaturesHumanTILs


data.seurat.tcells <- AddModuleScore_UCell(CD8_anno, features = CD8_naive_like, ncores = 4)

featnames <- paste0(names(signaturesHumanTILs), "_UCell")


FeaturePlot(data.seurat.tcells, features = CD8_naive_like, pt.size = 0.1, order = T)

VlnPlot(data.seurat.tcells, features = CD8_naive_like, pt.size = 0, split.by = "seurat_clusters")


CD8_naive_like <- c("CD8A",  "CD8B",  "CCR7",  "IL7R",  "SELL",  "TCF7",  "S1PR1", "LEF1")
FeaturePlot(data.seurat.tcells, features = "CD4", pt.size = 0.1, order = T)




CD8_anno_cancer <- subset(x = CD8_subset, Type == 'Tumor' & Outcome != 'Meta')
CD8_anno_normaltumor <- subset(x = CD8_subset, Type == 'Tumor' & Outcome != 'Meta')

DimPlot(CD8_subset, reduction = 'umap', split.by = "Type", raster=FALSE, label = F) +NoLegend()

DimPlot(CD8_anno_cancer, reduction = 'umap', split.by = "Outcome", raster=FALSE, label = F) +NoLegend()




dittoBoxPlot(CD8_anno, assay = "integrated", "ITGAE", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
             boxplot.width = 0.7, max=5) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))

dittoBoxPlot(CD8_anno, assay = "integrated", "ENTPD1", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
             boxplot.width = 0.7, max=5) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))


dittoBoxPlot(CD8_anno, assay = "integrated", "PDCD1", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
             boxplot.width = 0.7, max=5) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))

dittoBoxPlot(CD8_anno, assay = "integrated", "CTLA4", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
             boxplot.width = 0.7, max=5) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))

dittoBoxPlot(CD8_anno, assay = "integrated", "TIGIT", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
             boxplot.width = 0.7, max=5) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))

dittoBoxPlot(CD8_anno, assay = "integrated", "LAG3", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
             boxplot.width = 0.7, max=8) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))

dittoBoxPlot(CD8_anno, assay = "integrated", "HAVCR2", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
             boxplot.width = 0.7, max=5) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))




dittoBoxPlot(CD8_anno, assay = "integrated", "GZMB", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
             boxplot.width = 0.7, max=10) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))

dittoBoxPlot(CD8_anno, assay = "integrated", "IFNG", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
             boxplot.width = 0.7, max=5) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))

dittoBoxPlot(CD8_anno, assay = "integrated", "TNFRSF18", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
             boxplot.width = 0.7, max=5) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))

################################################################################## Annotation #####

result <- FindMarkers(CD8_anno, ident.1 = "3", ident.2 = "4")
result <- rownames_to_column(result, var = "rowname")


results <- FindMarkers(CD8_anno, ident.1 = "3", ident.2 = "4", min.pct = 0.5, logfc.threshold = log(2), only.pos = T)
results <- rownames_to_column(results, var = "rowname")
results_group <- as.character(results$rowname)
out <- SCpubr::do_FunctionalAnnotationPlot(genes = results_group,
                                           org.db = org.Hs.eg.db,
                                           min.overlap = 15, legend.position = "none")

out$BarPlot






dittoBoxPlot(CD8_subset, "ENTPD1", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
             boxplot.width = 0.7, max=5)


SCpubr::do_BoxPlot(sample = CD8_subset, feature = "ITGAE")
SCpubr::do_BoxPlot(sample = CD8_subset, feature = "CD69")



median.stat <- function(x){
  out <- quantile(x, probs = c(0.5))
  names(out) <- c("ymed")
  return(out) 
}

VlnPlot(CD8_subset, features = "CD69", raster=FALSE, pt.size = 0) + 
  stat_summary(fun.y = median.stat, geom='point', size = 10, colour = "black", shape = 95) 


VlnPlot(CD8_subset, features = "ITGAE", raster=FALSE, pt.size = 0.1)
VlnPlot(CD8_subset, features = "CD69", raster=FALSE, pt.size = 0.1)


VlnPlot(CD8_subset, features = "CTLA4", raster=FALSE, pt.size = 0.1)
VlnPlot(CD8_subset, features = "PDCD1", raster=FALSE, pt.size = 0.1)
VlnPlot(CD8_subset, features = "LAG3", raster=FALSE, pt.size = 0.1)


CD8_subset






CD8_subset.markers <- FindAllMarkers(CD8_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


plot_density(CD8_subset, "CD28")
plot_density(CD8_subset, "TCF7")
plot_density(CD8_subset, "GZMK")
plot_density(CD8_subset, "IL18R1")
plot_density(CD8_subset, "IL23A")
do_NebulosaPlot(CD8_subset, "IL23A", viridis_direction = -1, legend.position="none")





plot_density(CD8_subset, "CCR7")

do_NebulosaPlot(CD8_subset, "CD28", viridis_direction = 1, legend.position="none")


do_NebulosaPlot(CD8_subset, "PDCD1", viridis_direction = 1, pt.size = 1, plot_cell_borders = F, legend.position="none")
do_NebulosaPlot(CD8_subset, "TIGIT", viridis_direction = 1, pt.size = 1, plot_cell_borders = F, legend.position="none")
do_NebulosaPlot(CD8_subset, "LAG3", viridis_direction = 1, pt.size = 1, plot_cell_borders = F, legend.position="none")
do_NebulosaPlot(CD8_subset, "CTLA4", viridis_direction = -1, pt.size = 1, plot_cell_borders = F, legend.position="none")
do_NebulosaPlot(CD8_subset, "ITGAE", viridis_direction = 1, pt.size = 1, plot_cell_borders = F, legend.position="none")
do_NebulosaPlot(CD8_subset, "ENTPD1", viridis_direction = -1, pt.size = 1, plot_cell_borders = F, legend.position="none")

do_NebulosaPlot(CD8_subset, "CCR9", viridis_direction = 1, pt.size = 1, plot_cell_borders = F, legend.position="none")

do_NebulosaPlot(CD8_subset, "IFNG", viridis_direction = 1, pt.size = 1, plot_cell_borders = F, legend.position="none")
do_NebulosaPlot(CD8_subset, "TGF", viridis_direction = 1, pt.size = 1, plot_cell_borders = F, legend.position="none")



do_NebulosaPlot(CD8_subset, "IFNG", viridis_direction = 1, pt.size = 1, plot_cell_borders = F, legend.position="none")
do_NebulosaPlot(CD8_subset, "GZMB", viridis_direction = 1, pt.size = 1, plot_cell_borders = F, legend.position="none")
do_NebulosaPlot(CD8_subset, "GZMA", viridis_direction = 1, pt.size = 1, plot_cell_borders = F, legend.position="none")
do_NebulosaPlot(CD8_subset, "EOMES", viridis_direction = 1, pt.size = 1, plot_cell_borders = F, legend.position="none")
do_NebulosaPlot(CD8_subset, "TBX21", viridis_direction = 1, pt.size = 1, plot_cell_borders = F, legend.position="none")
do_NebulosaPlot(CD8_subset, "IL2", viridis_direction = 1, pt.size = 1, plot_cell_borders = F, legend.position="none")
do_NebulosaPlot(CD8_subset, "MKI67", viridis_direction = 1, pt.size = 1, plot_cell_borders = F, legend.position="none")


plot_density(CD8_subset, "IL2RA")
plot_density(CD8_subset, "GZMK")
plot_density(CD8_subset, "TNF")
plot_density(CD8_subset, "IFNG")




plot_density(CD8_subset, "IL7R")
plot_density(CD8_subset, "KLRG1")
plot_density(CD8_subset, "IL15RA")


plot_density(CD8_subset, "CCR4")



plot_density(CD8_subset, "GFPT2")

plot_density(CD8_subset, "ITGAE")
plot_density(CD8_subset, "ENTPD1")

FeaturePlot(CD8_subset, features = 'CCR7', raster=FALSE)
FeaturePlot(CD8_subset, features = 'CCR9', raster=FALSE)

cluster2.markers <- FindMarkers(CD8_subset, ident.1 = 2, ident.2 = c(0,1,3,4,5,6,7), only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = 'roc')



cluster0.markers <- FindMarkers(CD8_subset, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE, norm.method = NULL)
cluster1.markers <- FindMarkers(CD8_subset, ident.1 = 1, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE, norm.method = NULL)





dittoBoxPlot(CD8_anno, assay = "integrated", "CCR9", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
             boxplot.width = 0.7, min = -1, max=2) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))










do_NebulosaPlot(CD8_subset, "CCR7", viridis_direction = -1, legend.position="none")




CD8_subset.markers %>%
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_log2FC)

CD8_subset.markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) -> top50

CD8_subset_0.markers <- subset(CD8_subset.markers, cluster=="0")
write.csv(CD8_subset_0.markers, file = "D:/DATA/Gastric cancer/CD8_subset_0.csv")

CD8_subset_1.markers <- subset(CD8_subset.markers, cluster=="1")
write.csv(CD8_subset_1.markers, file = "D:/DATA/Gastric cancer/CD8_subset_1.csv")

CD8_subset_2.markers <- subset(CD8_subset.markers, cluster=="2")
write.csv(CD8_subset_2.markers, file = "D:/DATA/Gastric cancer/CD8_subset_2.csv")

CD8_subset_3.markers <- subset(CD8_subset.markers, cluster=="3")
write.csv(CD8_subset_3.markers, file = "D:/DATA/Gastric cancer/CD8_subset_3.csv")

CD8_subset_4.markers <- subset(CD8_subset.markers, cluster=="4")
write.csv(CD8_subset_4.markers, file = "D:/DATA/Gastric cancer/CD8_subset_4.csv")


DoHeatmap(CD8_subset, features = top50$gene, label=F) + NoLegend()


FeaturePlot(CD8_subset, features = 'GZMA', min.cutoff = '0', max.cutoff = '3', raster=FALSE)
FeaturePlot(CD8_subset, features = 'PDCD1', min.cutoff = '0', max.cutoff = '3', raster=FALSE)
FeaturePlot(CD8_subset, features = 'CCR7', min.cutoff = '0', max.cutoff = '1', raster=FALSE)
FeaturePlot(CD8_subset, features = 'CTLA4', min.cutoff = '0', max.cutoff = '3', raster=FALSE)
