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

Plasma_subset <- subset(x = immune_annotation, idents = "Plasma cells")

##################################### Plasma T cells clustering ######################################### OLD

Plasma_subset <- RunPCA(Plasma_subset)
Plasma_subset <- FindNeighbors(Plasma_subset, reduction = "pca", dims = 1:20)
Plasma_subset <- FindClusters(Plasma_subset, resolution = 0.2)
Plasma_subset <- RunUMAP(Plasma_subset, reduction = "pca", dims = 1:20)

DimPlot(Plasma_subset, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Plasma cells")

# Check outsiders
Plasma_subset <- subset(x = Plasma_subset, seurat_clusters != '5')

Plasma_subset <- RunPCA(Plasma_subset)
Plasma_subset <- FindNeighbors(Plasma_subset, reduction = "pca", dims = 1:20)
Plasma_subset <- FindClusters(Plasma_subset, resolution = 0.2)
Plasma_subset <- RunUMAP(Plasma_subset, reduction = "pca", dims = 1:20)

DimPlot(Plasma_subset, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Plasma cells")

Plasma_subset_anno <- Plasma_subset
DimPlot(Plasma_subset_anno, reduction = 'umap', raster=FALSE, label = F) + ggtitle("")+ NoLegend() +NoAxes()


Plasma_subset_anno <- RenameIdents(Plasma_subset_anno, `0`  = 'C0')
Plasma_subset_anno <- RenameIdents(Plasma_subset_anno, `1`  = 'C0')
Plasma_subset_anno <- RenameIdents(Plasma_subset_anno, `2`  = 'C1')
Plasma_subset_anno <- RenameIdents(Plasma_subset_anno, `3`  = 'C2')
Plasma_subset_anno <- RenameIdents(Plasma_subset_anno, `4`  = 'C1')
Plasma_subset_anno <- RenameIdents(Plasma_subset_anno, `5`  = 'C3')

DimPlot(Plasma_subset_anno, reduction = 'umap', raster=FALSE, label = F) + ggtitle("")+ NoLegend() +NoAxes()

view(Plasma_subset_anno@meta.data)
Plasma_subset_anno

Plasma_subset_anno@meta.data$seurat_clusters <- Plasma_subset_anno@active.ident

saveRDS(Plasma_subset_anno, file="D:/DATA/Gastric cancer/GS_subsets/Plasma_subset.rds")
Plasma_subset_anno <- readRDS(file="D:/DATA/Gastric cancer/GS_subsets/Plasma_subset.rds")
DimPlot(Plasma_subset_anno, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Plasma T cells") 
DimPlot(Plasma_subset_anno, reduction = 'umap', raster=FALSE, label = F) + ggtitle("")+NoLegend()+NoAxes()


##################################### Plasma T cells clustering ######################################### NEW
Plasma_subset_filtered_Trex <- RunPCA(Plasma_subset_filtered_Trex)
Plasma_subset_filtered_Trex <- FindNeighbors(Plasma_subset_filtered_Trex, reduction = "pca", dims = 1:50)
Plasma_subset_filtered_Trex <- FindClusters(Plasma_subset_filtered_Trex, resolution = 0.3)
Plasma_subset_filtered_Trex <- RunUMAP(Plasma_subset_filtered_Trex, reduction = "pca", dims = 1:50)

DimPlot(Plasma_subset_filtered_Trex, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Plasma T cells") 



############################################################################## Annotation #####

# Markers
Plasma_markers <- c("CCR7", "TCF7", "SELL", 
                 "ISG15", "IFIT3", "STAT1",
                 "GNLY", "PRF1",
                 "ENTPD1", "ITGAE",
                 "PDCD1", "TIGIT", "CTLA4",
                 "FOXP3", "EOMES", "TBX21", "GZMK", "GZMB"
)


# Dot plots
DotPlot(Plasma_subset_anno, features = Plasma_markers, col.min = 0,  col.max = 2)+ RotatedAxis() + xlab(NULL) + ylab(NULL) + theme(axis.text.x = element_text(face = "italic"))
SCpubr::do_ExpressionHeatmap(sample = Plasma_subset_anno,
                             features = Plasma_markers,
                             flip = F,
                             cluster_cols = F,
                             cluster_rows = F,
                             enforce_symmetry = TRUE,
                             use_viridis = FALSE,
                             max.cutoff = 3,
                             min.cutoff = 0)

# Dot plots
DotPlot(CD8_subset_anno, features = genes)+ RotatedAxis() + xlab(NULL) + ylab(NULL) + theme(axis.text.x = element_text(face = "italic"))
DotPlot(CD8_subset_anno, features = genes)+ RotatedAxis() + xlab(NULL) + ylab(NULL) +NoLegend() + 
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))

do_ExpressionHeatmap(sample = CD8_subset_anno,
                     features = genes,
                     viridis_direction = -1)


FeaturePlot(Plasma_subset_anno, reduction = "umap", features = "Cytotoxic_kNN")
VlnPlot(Plasma_subset_anno, features = "Cytotoxic_kNN", pt.size = 0, split.by = "seurat_clusters")+ NoLegend()

FeaturePlot(Plasma_subset_anno, reduction = "umap", features = "Exhaustion_kNN")
VlnPlot(Plasma_subset_anno, features = "Exhaustion_kNN", pt.size = 0, split.by = "seurat_clusters")+ NoLegend()


Plasma.markers <- FindAllMarkers(Plasma_subset_anno, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Plasma.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(Plasma_subset_anno, features = top10$gene) + NoLegend()


##################################### Stats #########################################

Plasma_subset_cancer <- subset(x = Plasma_subset_anno, Type == 'Tumor' & Outcome != 'Meta' & Patient_ID_SPC != F)
Plasma_subset_cancer_MSI <- subset(x = Plasma_subset_cancer, MSI_type == 'MSI')
Plasma_subset_cancer_MSS <- subset(x = Plasma_subset_cancer, MSI_type == 'MSS')

SCpubr::do_BarPlot(sample = Plasma_subset_cancer, 
                   group.by = "seurat_clusters", 
                   legend.position = "none", 
                   plot.title = "",  xlab = "",  ylab = "Number of cells per cluster", font.size = 11)


# 전체
SCpubr::do_BarPlot(sample = Plasma_subset_cancer, 
                   group.by = "seurat_clusters", 
                   split.by = "Patient_ID_SPC",
                   position = "Fill",legend.position = "bottom",
                   plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)

# 부분
SCpubr::do_BarPlot(sample = Plasma_subset_cancer_MSI, 
                   group.by = "seurat_clusters", 
                   split.by = "Patient_ID_SPC",
                   position = "Fill",legend.position = "bottom",
                   plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)

SCpubr::do_BarPlot(sample = Plasma_subset_cancer_MSS, 
                   group.by = "seurat_clusters", 
                   split.by = "Patient_ID_SPC",
                   position = "Fill",legend.position = "bottom",
                   plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)




####################################################################################################################################################

TCR.subset <- Plasma_subset_anno

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

Plasma_TCRtpye <- TCR.subset

cancer_subset <- subset(x = Plasma_TCRtpye, Type == 'Tumor' & Patient_ID_SPC != FALSE & Stage != 'Stage IV')
##################################### TCR type comparison #########################################

# Show as UMAP

DimPlot(Plasma_TCRtpye, reduction = 'umap', group.by = "TCR_type", raster=FALSE, label = F, pt.size = 1, order=T,
        cols = c('C0_EXP' = '#e69f00', 'C1_EXP' = '#ebe6dd', 'C2_EXP' = '#ebe6dd')) + NoLegend() +ggtitle("Total expanded clones") 
DimPlot(cancer_subset, reduction = 'umap', group.by = "TCR_type", raster=FALSE, label = F, pt.size = 1, order=T,
        cols = c('C0_EXP' = '#e69f00', 'C1_EXP' = '#ebe6dd', 'C2_EXP' = '#ebe6dd')) + NoLegend() +ggtitle("") + NoAxes()

DimPlot(Plasma_TCRtpye, reduction = 'umap', group.by = "TCR_type", raster=FALSE, label = F, pt.size = 1, order=T,
        cols = c('C0_EXP' = '#ebe6dd', 'C1_EXP' = '#56b4e9', 'C2_EXP' = '#ebe6dd')) + NoLegend() +ggtitle("Total expanded clones") 
DimPlot(cancer_subset, reduction = 'umap', group.by = "TCR_type", raster=FALSE, label = F, pt.size = 1, order=T,
        cols = c('C0_EXP' = '#ebe6dd', 'C1_EXP' = '#56b4e9', 'C2_EXP' = '#ebe6dd')) + NoLegend() +ggtitle("") + NoAxes()

DimPlot(Plasma_TCRtpye, reduction = 'umap', group.by = "TCR_type", raster=FALSE, label = F, pt.size = 1, order=T,
        cols = c('C0_EXP' = '#ebe6dd', 'C1_EXP' = '#ebe6dd', 'C2_EXP' = '#009e73')) + NoLegend() +ggtitle("Total expanded clones") 
DimPlot(cancer_subset, reduction = 'umap', group.by = "TCR_type", raster=FALSE, label = F, pt.size = 1, order=T,
        cols = c('C0_EXP' = '#ebe6dd', 'C1_EXP' = '#ebe6dd', 'C2_EXP' = '#009e73')) + NoLegend() +ggtitle("") + NoAxes()



dittoBarPlot(cancer_subset, "seurat_clusters", group.by = "TCR_type", scale = "count", main = "")+ 
  NoLegend()+ xlab(NULL)+ ylab(NULL)+
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))

dittoBarPlot(cancer_subset, "seurat_clusters", group.by = "TCR_type", scale = "percent", main = "")+ 
  NoLegend()+ xlab(NULL)+ ylab(NULL)+
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))


dittoBarPlot(cancer_subset, "seurat_clusters", group.by = "Patient_ID_SPC", scale = "percent", main = "")+ 
  NoLegend()+ xlab(NULL)+ ylab(NULL)+
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))

##################################### Stats #########################################
##################################### Dataframe for Outcome
md.immune_annotation<- cancer_subset@meta.data %>% as.data.table

length(which(md.immune_annotation$Patient_ID==68 & md.immune_annotation$Type=='Tumor'& md.immune_annotation$seurat_clusters==1))

# Vector of sample numbers
sample_numbers <- unique(cancer_subset@meta.data$Patient_ID_SPC)

# Vector of cluster names
cluster_vec <- unique(md.immune_annotation$seurat_clusters)
type_vec <- unique(md.immune_annotation$TCR_type)



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
      subset_length <- length(which(md.immune_annotation$Patient_ID_SPC == number & md.immune_annotation$TCR_type == type & md.immune_annotation$seurat_clusters == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.immune_annotation$Patient_ID_SPC == number & md.immune_annotation$TCR_type == type))
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
  C0_EXP_df <- subset_df %>% subset(type_vec == "C0_EXP")
  C1_EXP_df <- subset_df %>% subset(type_vec == "C1_EXP")
  C2_EXP_df <- subset_df %>% subset(type_vec == "C2_EXP")
  
  # Perform the unpaired t-test
  ttest01 <- t.test(C0_EXP_df$subset_fraction, C1_EXP_df$subset_fraction)
  ttest02 <- t.test(C0_EXP_df$subset_fraction, C2_EXP_df$subset_fraction)
  ttest12 <- t.test(C1_EXP_df$subset_fraction, C2_EXP_df$subset_fraction)
  
  
  #wtest <- wilcox.test(normal_df$subset_fraction, tumor_df$subset_fraction)
  
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
    geom_boxplot(color = "black", outlier.shape = NA)+
    scale_fill_manual(values = c("#e69f00", "#56b4e9", "#009e73")) +
    scale_y_continuous(labels = function(x) paste0(x * 100))+
    labs(title = NULL,
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
          legend.position = "none")
  
  #stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")
  #annotate("text", x = 1.5, y = 0.5, label = paste("p =", format(ttest01$p.value, digits = 3)))+
  #annotate("text", x = 2, y = 0.8, label = paste("p =", format(ttest02$p.value, digits = 3)))+
  #annotate("text", x = 2.5, y = 0.4, label = paste("p =", format(ttest12$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[which(cluster_vec == cluster)]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 1)

############## Stat add

for (cluster in cluster_vec) {
  # Subset the data for each type_vec
  subset_df <- output_df %>% subset(cluster_vec == cluster)
  C0_EXP_df <- subset_df %>% subset(type_vec == "C0_EXP")
  C1_EXP_df <- subset_df %>% subset(type_vec == "C1_EXP")
  C2_EXP_df <- subset_df %>% subset(type_vec == "C2_EXP")
  
  # Perform the unpaired t-test
  ttest01 <- t.test(C0_EXP_df$subset_fraction, C1_EXP_df$subset_fraction)
  ttest02 <- t.test(C0_EXP_df$subset_fraction, C2_EXP_df$subset_fraction)
  ttest12 <- t.test(C1_EXP_df$subset_fraction, C2_EXP_df$subset_fraction)
  
  
  #wtest <- wilcox.test(normal_df$subset_fraction, tumor_df$subset_fraction)
  
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
    geom_boxplot(color = "black", outlier.shape = NA)+
    scale_fill_manual(values = c("#e69f00", "#56b4e9", "#009e73")) +
    scale_y_continuous(labels = function(x) paste0(x * 100))+
    labs(title = cluster,
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
    annotate("text", x = 1.5, y = 0.5, label = paste("p =", format(ttest01$p.value, digits = 3)))+
    annotate("text", x = 2, y = 0.8, label = paste("p =", format(ttest02$p.value, digits = 3)))+
    annotate("text", x = 2.5, y = 0.4, label = paste("p =", format(ttest12$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[which(cluster_vec == cluster)]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 1)
