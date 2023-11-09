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
immune_annotation <- readRDS(file="D:/DATA/Gastric cancer/immune_cell_annotation.rds")
DimPlot(immune_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Initial")


##################################### Recur vs Non-recur : in R_CD4 T cells#########################################
R_CD4_subset <- subset(x = immune_annotation, clusters == 'Regulatory CD4 T cells')

R_CD4_subset <- RunPCA(R_CD4_subset)

R_CD4_subset <- FindNeighbors(R_CD4_subset, reduction = "pca", dims = 1:30)
R_CD4_subset <- FindClusters(R_CD4_subset, resolution = 0.3)
R_CD4_subset <- RunUMAP(R_CD4_subset, reduction = "pca", dims = 1:30)

DimPlot(R_CD4_subset, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("R_CD4 T cells")
DimPlot(R_CD4_subset, reduction = 'umap', split.by = "Outcome", raster=FALSE, label = T) + NoLegend() + ggtitle("R_CD4 T cells")

##################################### Stats #########################################
R_CD4_subset_outcome <- subset(x = R_CD4_subset, Type == 'Tumor' & Outcome != 'Meta')

##################################### Dataframe for Outcome
md.R_CD4_subset_outcome<- R_CD4_subset_outcome@meta.data %>% as.data.table

length(which(md.R_CD4_subset_outcome$Patient_ID==68 & md.R_CD4_subset_outcome$Outcome=='Recur'& md.R_CD4_subset_outcome$seurat_clusters==1))

# Vector of sample numbers
sample_numbers <- c(4, 5, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 51, 52, 53, 54, 56, 59, 60, 62, 64, 65, 66, 68, 71, 73, 75)

# Vector of cluster names
cluster_vec <- unique(md.R_CD4_subset_outcome$seurat_clusters)
type_vec <- c("Non_recur", "Recur")

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
      subset_length <- length(which(md.R_CD4_subset_outcome$Patient_ID_SPC == number & md.R_CD4_subset_outcome$Outcome == type & md.R_CD4_subset_outcome$seurat_clusters == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.R_CD4_subset_outcome$Patient_ID_SPC == number & md.R_CD4_subset_outcome$Outcome == type))
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
          legend.position = "none")+
    
    #stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")
    annotate("text", x = 1.5, y = 0.5, label = paste("p =", format(ttest$p.value, digits = 3)))
  #annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[which(cluster_vec == cluster)]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 1)

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

################################################################################## Annotation #####