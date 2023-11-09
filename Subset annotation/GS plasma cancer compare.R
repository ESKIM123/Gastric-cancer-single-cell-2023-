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
Plasma_subset_outcome <- subset(x = Plasma_subset, Type == 'Tumor' & Outcome != 'Meta'& Molecular_type == "GS" & Patient_ID_SPC != "NA")
immune_annotation_outcome <- subset(x = immune_annotation, Type == 'Tumor' & Outcome != 'Meta'& Molecular_type == "GS" & Patient_ID_SPC != "NA")

##################################### Dataframe for Outcome
md.Plasma_subset_outcome<- Plasma_subset_outcome@meta.data %>% as.data.table

length(which(md.Plasma_subset_outcome$Patient_ID==68 & md.Plasma_subset_outcome$Outcome=='Recur'& md.Plasma_subset_outcome$seurat_clusters==1))

# Vector of sample numbers
sample_numbers <- c(59, 66, 44, 41, 18, 54, 4, 56, 52, 68, 29, 7, 36, 28, 8)

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

################################################################################## Graph

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
  grid.text("Plasma cells", y = 0.05)

################################################################################## Graph --- clean

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
          plot.title = element_text(size = 8, hjust = 0.5), 
          legend.position = "none")
  
  #stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")
  #annotate("text", x = 1.5, y = 0.5, label = paste("p =", format(ttest$p.value, digits = 3)))
  #annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[which(cluster_vec == cluster)]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 1)

################################################################################## Graph Total ####
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
  grid.text("Plasma cells (Total immune)", y = 0.05)

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

################################################################################## Annotation #####
Plasma_subset

DimPlot(Plasma_subset, reduction = 'umap', raster=FALSE, label = F) + NoLegend()+NoAxes()

SCpubr::do_BarPlot(sample = Plasma_subset, 
                   group.by = "seurat_clusters", 
                   position = "stack",legend.position = "none",
                   plot.title = "", xlab = "",  ylab = "Number of cells per cluster", font.size = 11)
################################################################################## Annotation #####

VlnPlot(Plasma_subset, features = "CXCR3", raster=FALSE, pt.size = 0.1)
VlnPlot(Plasma_subset, features = "IGHG1", raster=FALSE, pt.size = 0.1)
VlnPlot(Plasma_subset, features = "CD38", raster=FALSE, pt.size = 0.1)
VlnPlot(Plasma_subset, features = "IFNG", raster=FALSE, pt.size = 0.1)

 

dittoBoxPlot(Plasma_subset, assay = "integrated", "KLF2", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
             boxplot.width = 0.7, max=1) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))

dittoBoxPlot(Plasma_subset, assay = "integrated", "IGHG1", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
             boxplot.width = 0.7, max=15) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))

dittoBoxPlot(Plasma_subset, assay = "integrated", "IGHA1", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F,
             boxplot.width = 0.7, max=20) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))

dittoBoxPlot(Plasma_subset, assay = "integrated", "IGHM", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F, 
             boxplot.width = 0.7, max=20) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))

dittoBoxPlot(Plasma_subset, assay = "integrated", "CXCR5", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F, 
             boxplot.width = 0.7, max=0) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))

dittoBoxPlot(Plasma_subset, assay = "integrated", "CXCL13", group.by = "seurat_clusters", plots = "boxplot", boxplot.show.outliers = F, 
             boxplot.width = 0.7, max=0) + NoLegend() +ggtitle("") +
  theme(axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20))
