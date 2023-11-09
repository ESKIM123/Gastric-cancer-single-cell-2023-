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
immune_annotation <- readRDS(file="D:/DATA/Gastric cancer/immune_cell_annotation_clinical_data_add.rds")
DimPlot(immune_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Initial")


##################################### Recur vs Non-recur : in Plasma cells#########################################
Plasma_subset <- subset(x = immune_annotation, clusters == 'Plasma cells')

Plasma_subset <- RunPCA(Plasma_subset)
Plasma_subset <- FindNeighbors(Plasma_subset, reduction = "pca", dims = 1:20)
Plasma_subset <- FindClusters(Plasma_subset, resolution = 0.2)
Plasma_subset <- RunUMAP(Plasma_subset, reduction = "pca", dims = 1:20)

DimPlot(Plasma_subset, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Plasma cells")
DimPlot(Plasma_subset, reduction = 'umap', split.by = "Outcome", raster=FALSE, label = T) + NoLegend() + ggtitle("Plasma cells") |
  FeaturePlot(Plasma_subset, features = 'HSPB1', raster=FALSE)

# Check outsiders
Plasma_subset <- subset(x = Plasma_subset, seurat_clusters != '5')

Plasma_subset <- RunPCA(Plasma_subset)
Plasma_subset <- FindNeighbors(Plasma_subset, reduction = "pca", dims = 1:20)
Plasma_subset <- FindClusters(Plasma_subset, resolution = 0.2)
Plasma_subset <- RunUMAP(Plasma_subset, reduction = "pca", dims = 1:20)

DimPlot(Plasma_subset, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Plasma cells")

##################################### Stats #########################################
Plasma_subset_outcome <- subset(x = Plasma_subset, Type == 'Tumor' & Outcome != 'Meta')
immune_annotation_outcome <- subset(x = immune_annotation, Type == 'Tumor' & Outcome != 'Meta')

##################################### Dataframe for Outcome
md.Plasma_subset_outcome<- Plasma_subset_outcome@meta.data %>% as.data.table

length(which(md.Plasma_subset_outcome$Patient_ID==68 & md.Plasma_subset_outcome$Outcome=='Recur'& md.Plasma_subset_outcome$seurat_clusters==1))

# Vector of sample numbers
sample_numbers <- c(4, 5, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 51, 52, 53, 54, 56, 59, 60, 62, 64, 65, 66, 68, 71, 73, 75)

# Vector of cluster names
cluster_vec <- unique(md.Plasma_subset_outcome$seurat_clusters)
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
      subset_length <- length(which(md.Plasma_subset_outcome$Patient_ID_SPC == number & md.Plasma_subset_outcome$Outcome == type & md.Plasma_subset_outcome$seurat_clusters == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.Plasma_subset_outcome$Patient_ID_SPC == number & md.Plasma_subset_outcome$Outcome == type))
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

################################################################################## Graph (In Plasma cells) ####

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
  grid.text("Plasma T cells", y = 0.05)

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
  #annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[which(cluster_vec == cluster)]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 1)


################################################################################## Graph (Total cells) ####
md.immune_annotation_outcome<- immune_annotation_outcome@meta.data %>% as.data.table

# Vector of cluster names
cluster_vec <- unique(md.Plasma_subset_outcome$seurat_clusters)
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
      subset_length <- length(which(md.Plasma_subset_outcome$Patient_ID_SPC == number & md.Plasma_subset_outcome$Outcome == type & md.Plasma_subset_outcome$seurat_clusters == cluster))
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
  grid.text("Plasma T cells (Total immune)", y = 0.05)

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
Plasma_subset

DimPlot(Plasma_subset, reduction = 'umap', raster=FALSE, label = F) + NoLegend()+NoAxes()

SCpubr::do_BarPlot(sample = Plasma_subset, 
                   group.by = "seurat_clusters", 
                   position = "stack",legend.position = "none",
                   plot.title = "", xlab = "",  ylab = "Number of cells per cluster", font.size = 11)

################################################################################## Save data ----

saveRDS(Plasma_subset, file="D:/DATA/Gastric cancer/Plasma_subset.rds")
Plasma_subset <- readRDS(file="D:/DATA/Gastric cancer/Plasma_subset.rds")

############################################################################## Annotation #####

Plasma_anno <- Plasma_subset

Plasma_anno.markers <- FindAllMarkers(Plasma_anno, only.pos = TRUE, logfc.threshold = 0.25, test.use = 'roc')
write.csv(Plasma_anno.markers, file="D:/DATA/Gastric cancer/Plasma_anno_markers.csv", row.names = FALSE)


signatures <- list(
  Plasma_naive_like = c("TCF7", "LEF1", "TXK", "KLF2", "SCML1", "CCR7", "IL7R", "EEF1A1", "SELL", "EEF1B2", "PLAC8", "TPT1", "ACTN1", "GIMAP7", "MAL", "LINC00861", "FHIT"),
  ISG_Plasma_cells = c("PLSCR1", "IRF7", "STAT1", "SP100", "IFI16", "TNFSF10", "TNFSF13B", "IFIT3", "IFIT1", "RSAD2", "IFI44L", "MX1", "ISG15", "IFI6", "MX2", "IFIT2", "CMPK2"),
  CCR6_Th17_Plasma_cells = c("BHLHE40", "HOPX", "ZFP36L2", "MYBL1", "RORA", "TNFSF13B", "CCL20", "Plasma0LG", "CCL5", "TIMP1", "FLT3LG", "IL7R", "CXCR4", "CCR6", "IFNGR1", 
                          "GPR35", "ANXA1", "KLRB1", "IL4I1", "VIM", "FKBP11", "CAPG", "ABCB1", "TPT1", "PERP", "DPP4"),
  IL26_Th17_Plasma_cells = c("RORC", "RORA", "CEBPD", "RUNX2", "ID2", "IL26", "IL17F", "IL17A", "TNFSF13B", "CCL20", "Plasma0LG", "GZMA", "CKLF", "IL22", "TNFSF14", "CXCR6", 
                          "CCR6", "IL23R", "IL17RE", "IL7R", "KLRB1", "CTSH", "CA10", "CAPG", "PTPN13", "MGAT4A", "ERN1", "IL4I1", "ANKRD28", "TMIGD2"),
  GZMK_Tem_Plasma_cells = c("EOMES", "LITAF", "NFATC2", "RUNX3", "NR4A2", "GZMK", "CCL4", "GZMA", "CCL5", "GZMH", "TNFSF9", "CCL3L3", "IFNG", "FASLG", "CCL3", "CD74", "CCR5", 
                         "CXCR4", "CXCR3", "Plasma", "NKG7", "CST7", "CRTAM", "CCL4L2", "ENC1", "ITM2C", "SLAMF7", "AOAH", "F2R", "DTHD1"),
  IL21_Tfh_Plasma_cells = c("NR3C1", "TOX2", "TOX", "TSHZ2", "RBPJ", "CXCL13", "IL21", "TNFSF8", "TNFSF11", "TNFSF4", "LIF", "Plasma0LG", "IL6ST", "CXCR5", "IL6R", "GNG4", "CD200",
                         "NMB", "IGFL2", "LHFP", "CPM", "BTLA", "FKBP5", "PDE7B", "ITM2A"),
  IFNG_Tfh_Th1_Plasma_cells = c("ZBED2", "ZEB2", "ID2", "BHLHE40", "RBPJ", "CXCL13", "CCL4", "IFNG", "GZMB", "CCL3", "GZMA", "CCL5", "CSF2", "IL21", "FAM3C", "CXCR6", "CCR5", 
                             "Plasma", "IL2RG", "CD74", "LAG3", "HAVCR2", "PDCD1", "KRT86", "MYO7A", "PTMS", "RDH10", "CCL4L2", "PDE7B", "DUSP4"),
  TEMRA_Plasma_cells = c("GZMA", "XCL2", "GZMM", "CCL3", "FASLG", "GZMK", "CX3CR1", "CMKLR1", "IL5RA", "CXCR2", "IL18RAP", "FGFBP2", "NKG7", "GNLY", "S1PR5", "KLRD1", "C1orf21", 
                      "PRF1", "PLEK", "CTSW", "PRSS23"),
  CREM_Tm_Plasma_cells = c("CREM","ZFP36L2","ZNF331","NR4A2","SKIL", "CCL5","GZMA","TGFB1","AREG","CXCL16","GZMK", "CXCR4","IFNGR1","GPR35","Plasma4","IL7R",
                        "ZFP36","FTH1","SRGN","FAM177A1","TNFAIP3","YPEL5","LMNA","ANXA1","PTGER4","CDKN1A"),
  CREM_Tm_Plasma_cells = c("CREM","ZFP36L2","ZNF331","NR4A2","SKIL", "CCL5","GZMA","TGFB1","AREG","CXCL16","GZMK", "CXCR4","IFNGR1","GPR35","Plasma4","IL7R",
                        "ZFP36","FTH1","SRGN","FAM177A1","TNFAIP3","YPEL5","LMNA","ANXA1","PTGER4","CDKN1A"),
  CCL5_Tm_Plasma_cells = c("ZFP36L2", "HOPX", "ID2", "BHLHE40", "NR4A2", "CCL5", "GZMA", "GZMK", "XCL1", "Plasma0LG", "CCL4", "XCL2", "GZMH", "IFNG",
                        "GZMM", "CXCR4", "IL7R", "IFNGR1", "Plasma4", "IL18RAP", "ANXA1", "ALOX5AP", "GLUL", "ITGA1", "PARP8", "PTGER4", "CD69", "ZFP36", "CD99", "CLU"),
  CAPG_Tm_Plasma_cells = c("HOPX", "BHLHE40", "ID2", "RBPJ", "PPARG", "GZMA", "CCL5", "CKLF", "TIMP1", "IL32", "TGFB1", "Plasma0LG", "HMGB1", "CSF2", "CXCR3", "CD74", 
                        "IL2RG", "CXCR6", "SH3BGRL3", "TMSB4X", "ACTB", "CAPG", "ALOX5AP", "CD52", "PFN1", "S100A11", "S100A4", "COTL1"),
  Th2_Plasma_cells = c("AREG", "IL4", "IL5", "IL13","CXCR4", "CCR4", "CCR8", "PTGDR2", "HAVCR1", "IL17RB", "IL1RL1","BATF", "GATA3", "IRF4", "STAT6"),
  Tfh_Plasma_cells = c("IL21","CXCR3", "CXCR5", "ICOS", "PDCD1","BATF", "BCL6", "MAF", "IRF4", "STAT3")
)


Plasma_anno <- AddModuleScore_UCell(Plasma_anno,features=signatures, name=NULL)
Plasma_anno <- SmoothKNN(Plasma_anno, signature.names = names(signatures), reduction="pca")

FeaturePlot(Plasma_anno, reduction = "umap", features = c("Plasma_naive_like","Plasma_naive_like_kNN"))
VlnPlot(Plasma_anno, features = "Plasma_naive_like_kNN", pt.size = 0, split.by = "seurat_clusters")

FeaturePlot(Plasma_anno, reduction = "umap", features = c("ISG_Plasma_cells","ISG_Plasma_cells_kNN"))
VlnPlot(Plasma_anno, features = "ISG_Plasma_cells_kNN", pt.size = 0, split.by = "seurat_clusters")

FeaturePlot(Plasma_anno, reduction = "umap", features = c("GZMK_Tem_Plasma_cells","GZMK_Tem_Plasma_cells_kNN"))
VlnPlot(Plasma_anno, features = "GZMK_Tem_Plasma_cells_kNN", pt.size = 0, split.by = "seurat_clusters")

FeaturePlot(Plasma_anno, reduction = "umap", features = c("CCR6_Th17_Plasma_cells","CCR6_Th17_Plasma_cells_kNN"))
VlnPlot(Plasma_anno, features = "CCR6_Th17_Plasma_cells_kNN", pt.size = 0, split.by = "seurat_clusters")

FeaturePlot(Plasma_anno, reduction = "umap", features = c("IL21_Tfh_Plasma_cells","IL21_Tfh_Plasma_cells_kNN"))
VlnPlot(Plasma_anno, features = "IL21_Tfh_Plasma_cells_kNN", pt.size = 0, split.by = "seurat_clusters")

FeaturePlot(Plasma_anno, reduction = "umap", features = c("IFNG_Tfh_Th1_Plasma_cells","IFNG_Tfh_Th1_Plasma_cells_kNN"))
VlnPlot(Plasma_anno, features = "IFNG_Tfh_Th1_Plasma_cells_kNN", pt.size = 0, split.by = "seurat_clusters")

FeaturePlot(Plasma_anno, reduction = "umap", features = c("TEMRA_Plasma_cells","TEMRA_Plasma_cells_kNN"))
VlnPlot(Plasma_anno, features = "TEMRA_Plasma_cells_kNN", pt.size = 0, split.by = "seurat_clusters")

FeaturePlot(Plasma_anno, reduction = "umap", features = c("CREM_Tm_Plasma_cells","CREM_Tm_Plasma_cells_kNN"))
VlnPlot(Plasma_anno, features = "CREM_Tm_Plasma_cells_kNN", pt.size = 0, split.by = "seurat_clusters")

FeaturePlot(Plasma_anno, reduction = "umap", features = c("CCL5_Tm_Plasma_cells","CCL5_Tm_Plasma_cells_kNN"))
VlnPlot(Plasma_anno, features = "CCL5_Tm_Plasma_cells_kNN", pt.size = 0, split.by = "seurat_clusters")

FeaturePlot(Plasma_anno, reduction = "umap", features = c("CAPG_Tm_Plasma_cells","CAPG_Tm_Plasma_cells_kNN"))
VlnPlot(Plasma_anno, features = "CAPG_Tm_Plasma_cells_kNN", pt.size = 0, split.by = "seurat_clusters")

FeaturePlot(Plasma_anno, reduction = "umap", features = c("Th2_Plasma_cells","Th2_Plasma_cells_kNN"))
VlnPlot(Plasma_anno, features = "Th2_Plasma_cells_kNN", pt.size = 0, split.by = "seurat_clusters")

FeaturePlot(Plasma_anno, reduction = "umap", features = c("Tfh_Plasma_cells","Tfh_Plasma_cells_kNN"))
VlnPlot(Plasma_anno, features = "Tfh_Plasma_cells_kNN", pt.size = 0, split.by = "seurat_clusters")

##################################################################################
i = "CD27"
FeaturePlot(Plasma_anno, features = i, min.cutoff = '0', max.cutoff = '3', raster=FALSE) |
  dittoBoxPlot(Plasma_anno, assay = "integrated", i, group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
               boxplot.width = 0.7, max=15) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))

i = "IGHA2"
FeaturePlot(Plasma_anno, features = i, min.cutoff = '0', max.cutoff = '3', raster=FALSE) |
  dittoBoxPlot(Plasma_anno, assay = "integrated", i, group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
               boxplot.width = 0.7, max=15) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))

i = "CD38"
FeaturePlot(Plasma_anno, features = i, min.cutoff = '0', max.cutoff = '3', raster=FALSE) |
  dittoBoxPlot(Plasma_anno, assay = "integrated", i, group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
               boxplot.width = 0.7, max=15) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))

i = "SDC1"
FeaturePlot(Plasma_anno, features = i, min.cutoff = '0', max.cutoff = '3', raster=FALSE) |
  dittoBoxPlot(Plasma_anno, assay = "integrated", i, group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
               boxplot.width = 0.7, max=15) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))

i = "IGKC"
FeaturePlot(Plasma_anno, features = i, min.cutoff = '0', max.cutoff = '3', raster=FALSE) |
  dittoBoxPlot(Plasma_anno, assay = "integrated", i, group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
               boxplot.width = 0.7, max=15) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))
 





##################################################################################
cluster2.markers <- FindMarkers(Plasma_anno, ident.1 = 1, ident.2 = 0, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25, test.use = 'roc')
write.csv(cluster2.markers, file="D:/DATA/Gastric cancer/Plasma_anno_markers_cluster2.csv", row.names = T)