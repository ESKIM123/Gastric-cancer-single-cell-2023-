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

##################################### Normal vs Tumors #########################################

DimPlot(immune_annotation, reduction = 'umap', group.by = "Type", raster=FALSE, label = F)+NoAxes()+NoLegend()+ ggtitle("")

SCpubr::do_BarPlot(sample = immune_annotation, 
                   group.by = "clusters", 
                   split.by = "Type",
                   position = "stack",legend.position = "none",
                   plot.title = "", xlab = "",  ylab = "Number of cells per cluster", font.size = 11)

SCpubr::do_BarPlot(sample = immune_annotation, 
                   group.by = "clusters", 
                   split.by = "Type",
                   position = "Fill",legend.position = "none",
                   plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)

SCpubr::do_BarPlot(sample = immune_annotation, 
                   group.by = "clusters", 
                   split.by = "Type",
                   position = "Fill",legend.position = "bottom",
                   plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)


cancer_subset <- subset(x = immune_annotation, Type == 'Tumor' & Patient_ID_SPC != FALSE)
SCpubr::do_BarPlot(sample = cancer_subset, 
                   group.by = "clusters", 
                   split.by = "Patient_ID_SPC",
                   position = "Fill",legend.position = "bottom",
                   plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)

Normal_subset <- subset(x = immune_annotation, Type == 'Normal' & Patient_ID_SPC != FALSE)
SCpubr::do_BarPlot(sample = Normal_subset, 
                   group.by = "clusters", 
                   split.by = "Patient_ID_SPC",
                   position = "Fill",legend.position = "bottom",
                   plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)

##################################### Stats #########################################
##################################### Dataframe for Outcome
md.immune_annotation<- immune_annotation@meta.data %>% as.data.table

length(which(md.immune_annotation$Patient_ID==68 & md.immune_annotation$Type=='Tumor'& md.immune_annotation$seurat_clusters==1))

# Vector of sample numbers
sample_numbers <- c(4, 5, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 51, 52, 53, 54, 56, 59, 60, 62, 64, 65, 66, 68, 71, 73, 75)

# Vector of cluster names
cluster_vec <- unique(md.immune_annotation$clusters)
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
      subset_length <- length(which(md.immune_annotation$Patient_ID_SPC == number & md.immune_annotation$Type == type & md.immune_annotation$clusters == cluster))
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
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
    geom_boxplot(color = "black", outlier.shape = NA)+
    scale_fill_manual(values = c("brown", "grey")) +
    labs(title = NULL,
         x = NULL, y = NULL) +
    #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          panel.border = element_rect(color = "black", fill = NA, size = 0.2),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
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




g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
  geom_point(size = 4, shape =21, color = "black")+
  scale_fill_manual(values = c("brown", "grey")) +
  #geom_line(aes(group = sample_numbers), color = "grey") +
  labs(title = paste(cluster),
       x = "Type", y = "Subset Fraction") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
  theme(panel.background = element_rect(fill = "#FFFFFF"),
        axis.line = element_line(color = "black", size = 0.2),
        axis.ticks = element_line(size = 0.2),
        axis.ticks.length = unit(0.1, "cm"),
        #axis.text.y = element_blank(),
        #axis.text.x = element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        legend.position = "none")+ 
  
  stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
  stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", size = 0.7, width=0.5, color="#03839C", alpha = 0.5)

################################################################################## Voilin plots


NK_subset <- subset(x = immune_annotation, seurat_clusters == '6')

tumor_features <- c("PDCD1", "TIGIT", "LAG3", "CTLA4", "HAVCR2")

immune_annotation.markers <- FindAllMarkers(immune_annotation, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
immune_annotation.markers %>%
  group_by(Type) %>%
  slice_max(n = 2, order_by = avg_log2FC)

immune_annotation.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(immune_annotation, features = top10$gene) + NoLegend()



SCpubr::do_BoxPlot(sample = immune_annotation,
                   feature = "PDCD1",
                   split.by = "Type")

SCpubr::do_BoxPlot(sample = immune_annotation,
                   feature = "TIGIT",
                   split.by = "Type")

SCpubr::do_BoxPlot(sample = immune_annotation,
                   feature = "LAG3",
                   split.by = "Type")

SCpubr::do_BoxPlot(sample = immune_annotation,
                   feature = "CTLA4",
                   split.by = "Type")

SCpubr::do_BoxPlot(sample = immune_annotation,
                   feature = "HAVCR2",
                   split.by = "Type")


##################################### Recur vs Non recur #########################################

cancer_subset <- subset(x = immune_annotation, Type == 'Tumor' & Outcome != 'Meta')

DimPlot(cancer_subset, reduction = 'umap', group.by = "Outcome", raster=FALSE, label = F)

do_BarPlot(sample = cancer_subset, 
           group.by = "Outcome",
           split.by = "clusters",
           plot.title = "Number of cells per sample in each cluster",
           position = "stack")


SCpubr::do_BarPlot(sample = cancer_subset, 
                   group.by = "clusters", 
                   split.by = "Outcome",
                   position = "fill",
                   plot.title = "Fraction of cells by cluster")

SCpubr::do_BarPlot(sample = cancer_subset, 
                   group.by = "clusters", 
                   split.by = "Patient_ID_SPC",
                   position = "fill",
                   plot.title = "Fraction of cells by cluster in cancer cells")

##################################### Stats #########################################
##################################### Dataframe for Outcome
md.cancer_subset <- cancer_subset@meta.data %>% as.data.table

# Vector of sample numbers
sample_numbers <- c(4, 5, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 51, 52, 53, 54, 56, 59, 60, 62, 64, 65, 66, 68, 71, 73, 75)

# Vector of cluster names
cluster_vec <- unique(md.cancer_subset$clusters)
type_vec <- c("Recur", "Non_recur")

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
      subset_length <- length(which(md.cancer_subset$Patient_ID_SPC == number & md.cancer_subset$Outcome == type & md.cancer_subset$clusters == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.cancer_subset$Patient_ID_SPC == number & md.cancer_subset$Outcome == type))
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
  Recur_df <- subset_df %>% subset(type_vec == "Recur")
  Non_recur_df <- subset_df %>% subset(type_vec == "Non_recur")
  # Perform the unpaired t-test
  ttest <- t.test(Recur_df$subset_fraction, Non_recur_df$subset_fraction)
  wtest <- wilcox.test(Recur_df$subset_fraction, Non_recur_df$subset_fraction)
  
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
    geom_point(size = 4, shape =21, color = "black")+
    scale_fill_manual(values = c("brown", "grey")) +
    #geom_line(aes(group = sample_numbers), color = "grey") +
    labs(title = paste(cluster),
         x = "Type", y = "Subset Fraction") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          axis.ticks.length = unit(0.1, "cm"),
          #axis.text.y = element_blank(),
          #axis.text.x = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          legend.position = "none")+ 
    
    stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
    stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", size = 0.7, width=0.5, color="#03839C", alpha = 0.5)+
    annotate("text", x = 1.5, y = 0.9, label = paste("t.p =", format(ttest$p.value, digits = 3)))+
    annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[which(cluster_vec == cluster)]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 3)

################################################################################## Voilin plots

tumor_markers <- c("PDCD1", "TIGIT", "LAG3", "CTLA4", "HAVCR2", "TNFRSF18")

plot_list <- list()

for (markers in tumor_markers) {

p2 <- SCpubr::do_BoxPlot(sample = cancer_subset, legend.position = 'none', font.size = 10,
                   feature = markers,
                   split.by = "Outcome")

# Add the plot to the list
plot_list[[which(tumor_markers == markers)]] <- p2
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 3)


##################################### CD3 T cell subset #########################################

immune_noannotation <- readRDS(file="D:/DATA/Gastric cancer/immune_cell_UMAP.rds")

CD3_subset <- subset(x = immune_noannotation, seurat_clusters == '6' | seurat_clusters == '4' 
                     | seurat_clusters == '0' | seurat_clusters == '1' | seurat_clusters == '7' | seurat_clusters == '9')

CD3_subset <- RunPCA(CD3_subset, verbose = FALSE)

CD3_subset <- FindNeighbors(CD3_subset, reduction = "pca", dims = 1:40)
CD3_subset <- FindClusters(CD3_subset, resolution = 0.3)
CD3_subset <- RunUMAP(CD3_subset, reduction = "pca", dims = 1:40)

DefaultAssay(CD3_subset) <- "ADT"
DimPlot(CD3_subset, reduction = 'umap', raster=FALSE, label = T) + NoLegend() + ggtitle("CD3 T cells") |
  FeaturePlot(CD3_subset, features = 'CD4', min.cutoff = '1', max.cutoff = '3', raster=FALSE) |
  FeaturePlot(CD3_subset, features = 'CD8', min.cutoff = '1', max.cutoff = '3', raster=FALSE) 
DefaultAssay(CD3_subset) <- "integrated"


######################## Change dimensionality and visualize UMAP ####################
for (i in 10:20) {
  CD3_subset <- FindNeighbors(CD3_subset, reduction = "pca", dims = 1:i)
  CD3_subset <- FindClusters(CD3_subset, resolution = 0.3)
  CD3_subset <- RunUMAP(CD3_subset, reduction = "pca", dims = 1:i) 
  
  DefaultAssay(CD3_subset) <- "ADT"
  # plot the results and save as PNG file
  plot_title <- paste("Unsupervised clustering 1:", i, sep = "")
  png_CD3_i <- DimPlot(CD3_subset, reduction = 'umap', raster=FALSE, label = T) + NoLegend() + ggtitle("CD3 T cells") |
    FeaturePlot(CD3_subset, features = 'CD4', min.cutoff = '1', max.cutoff = '3', raster=FALSE) |
    FeaturePlot(CD3_subset, features = 'CD8', min.cutoff = '1', max.cutoff = '3', raster=FALSE) 
  
  DefaultAssay(CD3_subset) <- "integrated"
  ggsave(paste0("png_CD3_", i, ".png"), png_CD3_i, width = 24, height = 7, units = "in")
}












KLRB1

DimPlot(CD3_subset, reduction = 'umap', raster=FALSE, label = T) + NoLegend() + ggtitle("CD3 T cells") |
  Nebulosa::plot_density(CD3_subset, features = "TRAV1-2",  pal = "viridis")

DimPlot(CD3_subset, reduction = 'umap', raster=FALSE, label = T) + NoLegend() + ggtitle("CD3 T cells") |
  Nebulosa::plot_density(CD3_subset, features = "TRDC",  pal = "viridis")

DimPlot(CD3_subset, reduction = 'umap', raster=FALSE, label = T) + NoLegend() + ggtitle("CD3 T cells") |
  Nebulosa::plot_density(CD3_subset, features = "CCR7",  pal = "viridis")


SCpubr::do_NebulosaPlot(CD3_subset, features = c("IL23R", "CCR6", "KLRB1", "RORC"),joint = TRUE, return_only_joint = TRUE)
SCpubr::do_NebulosaPlot(CD3_subset, features = c("CCR7", "IL7R", "TCF7", "SELL"),joint = TRUE, return_only_joint = TRUE)
SCpubr::do_NebulosaPlot(CD3_subset, features = c("CXCR5", "BCL6", "CD200", "BTLA", "CXCL13"),joint = TRUE, return_only_joint = TRUE)
SCpubr::do_NebulosaPlot(CD3_subset, features = c("IFNG", "CXCR3", "PDCD1"),joint = TRUE, return_only_joint = TRUE)
SCpubr::do_NebulosaPlot(CD3_subset, features = c("IL17A", "IL22", "IL23R"),joint = TRUE, return_only_joint = TRUE)
SCpubr::do_NebulosaPlot(CD3_subset, features = c("FOXP3", "IL3RA", "TNFRSF18"),joint = TRUE, return_only_joint = TRUE)
SCpubr::do_NebulosaPlot(CD3_subset, features = c("MKI67", "STMN1"),joint = TRUE, return_only_joint = TRUE)




SCpubr::do_NebulosaPlot(CD3_subset, features = c("HSPA6", "HSPA1B", "HSPA1A", "HSPH1"))



(IL-23R), CCR6, CD161, and the transcription factor retinoic acid-related orphan receptor C (RORC),


FeaturePlot(CD3_subset, features = 'TRAV1-2', min.cutoff = '1', max.cutoff = '3', raster=FALSE)
FeaturePlot(CD3_subset, features = 'TRDC', min.cutoff = '1', max.cutoff = '3', raster=FALSE)
FeaturePlot(CD3_subset, features = 'FOXP3', min.cutoff = '1', max.cutoff = '3', raster=FALSE)
FeaturePlot(CD3_subset, features = 'KLRB1', min.cutoff = '1', max.cutoff = '3', raster=FALSE)
FeaturePlot(CD3_subset, features = 'CCR7', min.cutoff = '1', max.cutoff = '1', raster=FALSE)
FeaturePlot(CD3_subset, features = 'CTLA4', min.cutoff = '1', max.cutoff = '3', raster=FALSE)

##################################### CD3 : recur vs non recur #########################################

cancer_CD3_subset <- subset(x = CD3_subset, Type == 'Tumor' & Outcome != 'Meta')

DimPlot(cancer_CD3_subset, reduction = 'umap', split.by = "Outcome", raster=FALSE, label = T) + NoLegend()+ ggtitle("CD8 T cells")
DimPlot(cancer_CD3_subset, reduction = 'umap', group.by = "Outcome", raster=FALSE, label = F) + ggtitle("CD8 T cells")

do_BarPlot(sample = cancer_CD3_subset, 
           group.by = "Outcome",
           split.by = "seurat_clusters",
           plot.title = "Number of cells per sample in each cluster",
           position = "stack")


SCpubr::do_BarPlot(sample = cancer_CD3_subset, 
                   group.by = "seurat_clusters", 
                   split.by = "Outcome",
                   position = "fill",
                   plot.title = "Fraction of cells by cluster")

##################################### Stats #########################################
##################################### Dataframe for Outcome
md.cancer_CD3_subset <- cancer_CD3_subset@meta.data %>% as.data.table

# Vector of sample numbers
sample_numbers <- c(4, 5, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 51, 52, 53, 54, 56, 59, 60, 62, 64, 65, 66, 68, 71, 73, 75)

# Vector of cluster names
cluster_vec <- unique(md.cancer_CD3_subset$seurat_clusters)
type_vec <- c("Recur", "Non_recur")

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
      subset_length <- length(which(md.cancer_CD3_subset$Patient_ID_SPC == number & md.cancer_CD3_subset$Outcome == type & md.cancer_CD3_subset$seurat_clusters == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.cancer_CD3_subset$Patient_ID_SPC == number & md.cancer_CD3_subset$Outcome == type))
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
  Recur_df <- subset_df %>% subset(type_vec == "Recur")
  Non_recur_df <- subset_df %>% subset(type_vec == "Non_recur")
  # Perform the unpaired t-test
  ttest <- t.test(Recur_df$subset_fraction, Non_recur_df$subset_fraction)
  wtest <- wilcox.test(Recur_df$subset_fraction, Non_recur_df$subset_fraction)
  
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
    geom_point(size = 4, shape =21, color = "black")+
    scale_fill_manual(values = c("brown", "grey")) +
    #geom_line(aes(group = sample_numbers), color = "grey") +
    labs(title = paste(cluster),
         x = "Type", y = "Subset Fraction") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          axis.ticks.length = unit(0.1, "cm"),
          #axis.text.y = element_blank(),
          #axis.text.x = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          legend.position = "none")+ 
    
    stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
    stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", size = 0.7, width=0.5, color="#03839C", alpha = 0.5)+
    annotate("text", x = 1.5, y = 0.9, label = paste("t.p =", format(ttest$p.value, digits = 3)))+
    annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[which(cluster_vec == cluster)]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 3)


##################################### CD8 T cell subset #########################################

immune_noannotation <- readRDS(file="D:/DATA/Gastric cancer/immune_cell_UMAP.rds")

CD8_subset <- subset(x = immune_noannotation, seurat_clusters == '4' | seurat_clusters == '0')

CD8_subset <- RunPCA(CD8_subset, verbose = FALSE)

CD8_subset <- FindNeighbors(CD8_subset, reduction = "pca", dims = 1:50)
CD8_subset <- FindClusters(CD8_subset, resolution = 0.3)
CD8_subset <- RunUMAP(CD8_subset, reduction = "pca", dims = 1:50)

DimPlot(CD8_subset, reduction = 'umap', raster=FALSE, label = T) + NoLegend() + ggtitle("CD8 T cells")
FeaturePlot(CD8_subset, features = 'TRAV1-2', min.cutoff = '0', max.cutoff = '3', raster=FALSE)
FeaturePlot(CD8_subset, features = 'TRDC', min.cutoff = '0', max.cutoff = '3', raster=FALSE)
FeaturePlot(CD8_subset, features = 'GZMA', min.cutoff = '0', max.cutoff = '3', raster=FALSE)
FeaturePlot(CD8_subset, features = 'PDCD1', min.cutoff = '0', max.cutoff = '3', raster=FALSE)
FeaturePlot(CD8_subset, features = 'CCR7', min.cutoff = '0', max.cutoff = '1', raster=FALSE)
FeaturePlot(CD8_subset, features = 'CTLA4', min.cutoff = '0', max.cutoff = '3', raster=FALSE)


sample <- CD8_subset

# Set the identities correctly.
Seurat::Idents(CD8_subset) <- CD8_subset$seurat_clusters
de_genes <- tibble::tibble(Seurat::FindAllMarkers(object = CD8_subset))

p <- SCpubr::do_GroupwiseDEPlot(sample = CD8_subset,
                                de_genes = de_genes,
                                group.by = c("seurat_clusters", 
                                             "Outcome"),
                                row_title_expression = c("",
                                                         "Outcome"))

p

##################################### CD8 : recur vs non recur #########################################

cancer_CD8_subset <- subset(x = CD8_subset, Type == 'Tumor' & Outcome != 'Meta')

DimPlot(cancer_CD8_subset, reduction = 'umap', split.by = "Outcome", raster=FALSE, label = T) + NoLegend()+ ggtitle("CD8 T cells")
DimPlot(cancer_CD8_subset, reduction = 'umap', group.by = "Outcome", raster=FALSE, label = F) + ggtitle("CD8 T cells")

do_BarPlot(sample = cancer_CD8_subset, 
           group.by = "Outcome",
           split.by = "seurat_clusters",
           plot.title = "Number of cells per sample in each cluster",
           position = "stack")


SCpubr::do_BarPlot(sample = cancer_CD8_subset, 
                   group.by = "seurat_clusters", 
                   split.by = "Outcome",
                   position = "fill",
                   plot.title = "Fraction of cells by cluster")

##################################### Stats #########################################
##################################### Dataframe for Outcome
md.cancer_CD8_subset <- cancer_CD8_subset@meta.data %>% as.data.table

# Vector of sample numbers
sample_numbers <- c(4, 5, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 51, 52, 53, 54, 56, 59, 60, 62, 64, 65, 66, 68, 71, 73, 75)

# Vector of cluster names
cluster_vec <- unique(md.cancer_CD8_subset$seurat_clusters)
type_vec <- c("Recur", "Non_recur")

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
      subset_length <- length(which(md.cancer_CD8_subset$Patient_ID_SPC == number & md.cancer_CD8_subset$Outcome == type & md.cancer_CD8_subset$seurat_clusters == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.cancer_CD8_subset$Patient_ID_SPC == number & md.cancer_CD8_subset$Outcome == type))
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
  Recur_df <- subset_df %>% subset(type_vec == "Recur")
  Non_recur_df <- subset_df %>% subset(type_vec == "Non_recur")
  # Perform the unpaired t-test
  ttest <- t.test(Recur_df$subset_fraction, Non_recur_df$subset_fraction)
  wtest <- wilcox.test(Recur_df$subset_fraction, Non_recur_df$subset_fraction)
  
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
    geom_point(size = 4, shape =21, color = "black")+
    scale_fill_manual(values = c("brown", "grey")) +
    #geom_line(aes(group = sample_numbers), color = "grey") +
    labs(title = paste(cluster),
         x = "Type", y = "Subset Fraction") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          axis.ticks.length = unit(0.1, "cm"),
          #axis.text.y = element_blank(),
          #axis.text.x = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          legend.position = "none")+ 
    
    stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
    stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", size = 0.7, width=0.5, color="#03839C", alpha = 0.5)+
    annotate("text", x = 1.5, y = 0.9, label = paste("t.p =", format(ttest$p.value, digits = 3)))+
    annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[which(cluster_vec == cluster)]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 3)


##################################### CD4 T cell subset #########################################

CD4_subset <- subset(x = immune_noannotation, seurat_clusters == '1' | seurat_clusters == '7')

CD4_subset <- RunPCA(CD4_subset, verbose = FALSE)

CD4_subset <- FindNeighbors(CD4_subset, reduction = "pca", dims = 1:50)
CD4_subset <- FindClusters(CD4_subset, resolution = 0.3)
CD4_subset <- RunUMAP(CD4_subset, reduction = "pca", dims = 1:50)

DimPlot(CD4_subset, reduction = 'umap', raster=FALSE, label = T) + NoLegend() + ggtitle("CD4 T cells")
FeaturePlot(CD4_subset, features = 'FOXP3', min.cutoff = '0', max.cutoff = '3', raster=FALSE)
FeaturePlot(CD4_subset, features = 'TRAV1-2', min.cutoff = '0', max.cutoff = '3', raster=FALSE)
FeaturePlot(CD4_subset, features = 'TRDC', min.cutoff = '0', max.cutoff = '3', raster=FALSE)
FeaturePlot(CD4_subset, features = 'GZMB', min.cutoff = '0', max.cutoff = '3', raster=FALSE)
FeaturePlot(CD4_subset, features = 'PDCD1', min.cutoff = '0', max.cutoff = '3', raster=FALSE)
FeaturePlot(CD4_subset, features = 'CCR7', min.cutoff = '0', max.cutoff = '1', raster=FALSE)
FeaturePlot(CD4_subset, features = 'CTLA4', min.cutoff = '0', max.cutoff = '3', raster=FALSE)
FeaturePlot(CD4_subset, features = 'ANLN', min.cutoff = '0', max.cutoff = '3', raster=FALSE)

result <- FindMarkers(CD4_subset, ident.1 = "5", ident.2 = "1", min.pct = 0.5, logfc.threshold = log(2))



SCpubr::do_NebulosaPlot(CD4_subset, features = c("CCR7", "IL7R", "TCF7", "SELL"))
SCpubr::do_NebulosaPlot(CD4_subset, features = c("CXCR5", "BCL6", "CD200", "BTLA", "CXCL13"))
SCpubr::do_NebulosaPlot(CD4_subset, features = c("IFNG", "CXCR3", "PDCD1"))
SCpubr::do_NebulosaPlot(CD4_subset, features = c("IL17A", "IL17F", "IL22", "IL23R"))
SCpubr::do_NebulosaPlot(CD4_subset, features = c("FOXP3", "IL3RA", "TNFRSF18"))
SCpubr::do_NebulosaPlot(CD4_subset, features = c("MKI67", "STMN1"))
SCpubr::do_NebulosaPlot(CD4_subset, features = c("HSPA6", "HSPA1B", "HSPA1A", "HSPH1"))













##################################### CD4 : recur vs non recur #########################################

cancer_CD4_subset <- subset(x = CD4_subset, Type == 'Tumor' & Outcome != 'Meta')

DimPlot(cancer_CD4_subset, reduction = 'umap', group.by = "Outcome", raster=FALSE, label = F) + ggtitle("CD4 T cells")
DimPlot(cancer_CD4_subset, reduction = 'umap', split.by = "Outcome", raster=FALSE, label = T) + NoLegend()+ ggtitle("CD4 T cells")


do_BarPlot(sample = cancer_CD4_subset, 
           group.by = "Outcome",
           split.by = "seurat_clusters",
           plot.title = "Number of cells per sample in each cluster",
           position = "stack")


SCpubr::do_BarPlot(sample = cancer_CD4_subset, 
                   group.by = "seurat_clusters", 
                   split.by = "Outcome",
                   position = "fill",
                   plot.title = "Fraction of cells by cluster")

##################################### Stats #########################################
##################################### Dataframe for Outcome
md.cancer_CD4_subset <- cancer_CD4_subset@meta.data %>% as.data.table

# Vector of sample numbers
sample_numbers <- c(4, 5, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 51, 52, 53, 54, 56, 59, 60, 62, 64, 65, 66, 68, 71, 73, 75)

# Vector of cluster names
cluster_vec <- unique(md.cancer_CD4_subset$seurat_clusters)
type_vec <- c("Recur", "Non_recur")

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
      subset_length <- length(which(md.cancer_CD4_subset$Patient_ID_SPC == number & md.cancer_CD4_subset$Outcome == type & md.cancer_CD4_subset$seurat_clusters == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.cancer_CD4_subset$Patient_ID_SPC == number & md.cancer_CD4_subset$Outcome == type))
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
  Recur_df <- subset_df %>% subset(type_vec == "Recur")
  Non_recur_df <- subset_df %>% subset(type_vec == "Non_recur")
  # Perform the unpaired t-test
  ttest <- t.test(Recur_df$subset_fraction, Non_recur_df$subset_fraction)
  wtest <- wilcox.test(Recur_df$subset_fraction, Non_recur_df$subset_fraction)
  
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
    geom_point(size = 4, shape =21, color = "black")+
    scale_fill_manual(values = c("brown", "grey")) +
    #geom_line(aes(group = sample_numbers), color = "grey") +
    labs(title = paste(cluster),
         x = "Type", y = "Subset Fraction") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          axis.ticks.length = unit(0.1, "cm"),
          #axis.text.y = element_blank(),
          #axis.text.x = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          legend.position = "none")+ 
    
    stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
    stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", size = 0.7, width=0.5, color="#03839C", alpha = 0.5)+
    annotate("text", x = 1.5, y = 0.9, label = paste("t.p =", format(ttest$p.value, digits = 3)))+
    annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[which(cluster_vec == cluster)]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 3)



##################################### NK cell subset #########################################

NK_subset <- subset(x = immune_noannotation, seurat_clusters == '6')

NK_subset <- RunPCA(NK_subset, verbose = FALSE)

NK_subset <- FindNeighbors(NK_subset, reduction = "pca", dims = 1:50)
NK_subset <- FindClusters(NK_subset, resolution = 0.3)
NK_subset <- RunUMAP(NK_subset, reduction = "pca", dims = 1:50)

DimPlot(NK_subset, reduction = 'umap', raster=FALSE, label = T) + NoLegend() + ggtitle("NK cells")
FeaturePlot(NK_subset, features = 'CD8A', min.cutoff = '0', max.cutoff = '3', raster=FALSE)
FeaturePlot(NK_subset, features = 'NCAM1', min.cutoff = '0', max.cutoff = '3', raster=FALSE)

##################################### NK : recur vs non recur #########################################

cancer_NK_subset <- subset(x = NK_subset, Type == 'Tumor' & Outcome != 'Meta')

DimPlot(cancer_NK_subset, reduction = 'umap', group.by = "Outcome", raster=FALSE, label = F) + ggtitle("NK cells")
DimPlot(cancer_NK_subset, reduction = 'umap', split.by = "Outcome", raster=FALSE, label = T) + NoLegend()+ ggtitle("NK cells")


do_BarPlot(sample = cancer_NK_subset, 
           group.by = "Outcome",
           split.by = "seurat_clusters",
           plot.title = "Number of cells per sample in each cluster",
           position = "stack")


SCpubr::do_BarPlot(sample = cancer_NK_subset, 
                   group.by = "seurat_clusters", 
                   split.by = "Outcome",
                   position = "fill",
                   plot.title = "Fraction of cells by cluster")

##################################### Stats #########################################
##################################### Dataframe for Outcome
md.cancer_NK_subset <- cancer_NK_subset@meta.data %>% as.data.table

# Vector of sample numbers
sample_numbers <- c(4, 5, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 51, 52, 53, 54, 56, 59, 60, 62, 64, 65, 66, 68, 71, 73, 75)

# Vector of cluster names
cluster_vec <- unique(md.cancer_NK_subset$seurat_clusters)
type_vec <- c("Recur", "Non_recur")

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
      subset_length <- length(which(md.cancer_NK_subset$Patient_ID_SPC == number & md.cancer_NK_subset$Outcome == type & md.cancer_NK_subset$seurat_clusters == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.cancer_NK_subset$Patient_ID_SPC == number & md.cancer_NK_subset$Outcome == type))
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
  Recur_df <- subset_df %>% subset(type_vec == "Recur")
  Non_recur_df <- subset_df %>% subset(type_vec == "Non_recur")
  # Perform the unpaired t-test
  ttest <- t.test(Recur_df$subset_fraction, Non_recur_df$subset_fraction)
  wtest <- wilcox.test(Recur_df$subset_fraction, Non_recur_df$subset_fraction)
  
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
    geom_point(size = 4, shape =21, color = "black")+
    scale_fill_manual(values = c("brown", "grey")) +
    #geom_line(aes(group = sample_numbers), color = "grey") +
    labs(title = paste(cluster),
         x = "Type", y = "Subset Fraction") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          axis.ticks.length = unit(0.1, "cm"),
          #axis.text.y = element_blank(),
          #axis.text.x = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          legend.position = "none")+ 
    
    stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
    stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", size = 0.7, width=0.5, color="#03839C", alpha = 0.5)+
    annotate("text", x = 1.5, y = 0.9, label = paste("t.p =", format(ttest$p.value, digits = 3)))+
    annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[which(cluster_vec == cluster)]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 3)


##################################### B cell subset #########################################

B_subset <- subset(x = immune_noannotation, seurat_clusters == '2')

B_subset <- RunPCA(B_subset, verbose = FALSE)

B_subset <- FindNeighbors(B_subset, reduction = "pca", dims = 1:50)
B_subset <- FindClusters(B_subset, resolution = 0.3)
B_subset <- RunUMAP(B_subset, reduction = "pca", dims = 1:50)

DimPlot(B_subset, reduction = 'umap', raster=FALSE, label = T) + NoLegend() + ggtitle("B cells")
FeaturePlot(B_subset, features = 'CD19', min.cutoff = '0', max.cutoff = '3', raster=FALSE)
FeaturePlot(B_subset, features = 'CD69', min.cutoff = '0', max.cutoff = '3', raster=FALSE)
FeaturePlot(B_subset, features = 'CD80', min.cutoff = '0', max.cutoff = '3', raster=FALSE)
FeaturePlot(B_subset, features = 'IL2RA', min.cutoff = '0', max.cutoff = '3', raster=FALSE)
FeaturePlot(B_subset, features = 'SDC1', min.cutoff = '0', max.cutoff = '3', raster=FALSE)
FeaturePlot(B_subset, features = 'HSPB1', min.cutoff = '0', max.cutoff = '3', raster=FALSE)

##################################### B : recur vs non recur #########################################

cancer_B_subset <- subset(x = B_subset, Type == 'Tumor' & Outcome != 'Meta')

DimPlot(cancer_B_subset, reduction = 'umap', group.by = "Outcome", raster=FALSE, label = F) + ggtitle("B cells")
DimPlot(cancer_B_subset, reduction = 'umap', split.by = "Outcome", raster=FALSE, label = T) + NoLegend()+ ggtitle("B cells")


do_BarPlot(sample = cancer_B_subset, 
           group.by = "Outcome",
           split.by = "seurat_clusters",
           plot.title = "Number of cells per sample in each cluster",
           position = "stack")


SCpubr::do_BarPlot(sample = cancer_B_subset, 
                   group.by = "seurat_clusters", 
                   split.by = "Outcome",
                   position = "fill",
                   plot.title = "Fraction of cells by cluster")

##################################### Stats #########################################
##################################### Dataframe for Outcome
md.cancer_B_subset <- cancer_B_subset@meta.data %>% as.data.table

# Vector of sample numbers
sample_numbers <- c(4, 5, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 51, 52, 53, 54, 56, 59, 60, 62, 64, 65, 66, 68, 71, 73, 75)

# Vector of cluster names
cluster_vec <- unique(md.cancer_B_subset$seurat_clusters)
type_vec <- c("Recur", "Non_recur")

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
      subset_length <- length(which(md.cancer_B_subset$Patient_ID_SPC == number & md.cancer_B_subset$Outcome == type & md.cancer_B_subset$seurat_clusters == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.cancer_B_subset$Patient_ID_SPC == number & md.cancer_B_subset$Outcome == type))
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
  Recur_df <- subset_df %>% subset(type_vec == "Recur")
  Non_recur_df <- subset_df %>% subset(type_vec == "Non_recur")
  # Perform the unpaired t-test
  ttest <- t.test(Recur_df$subset_fraction, Non_recur_df$subset_fraction)
  wtest <- wilcox.test(Recur_df$subset_fraction, Non_recur_df$subset_fraction)
  
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
    geom_point(size = 4, shape =21, color = "black")+
    scale_fill_manual(values = c("brown", "grey")) +
    #geom_line(aes(group = sample_numbers), color = "grey") +
    labs(title = paste(cluster),
         x = "Type", y = "Subset Fraction") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          axis.ticks.length = unit(0.1, "cm"),
          #axis.text.y = element_blank(),
          #axis.text.x = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          legend.position = "none")+ 
    
    stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
    stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", size = 0.7, width=0.5, color="#03839C", alpha = 0.5)+
    annotate("text", x = 1.5, y = 0.9, label = paste("t.p =", format(ttest$p.value, digits = 3)))+
    annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[which(cluster_vec == cluster)]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 3)


##################################### Plasma cell subset #########################################

Plasma_subset <- subset(x = immune_noannotation, seurat_clusters == '5')

Plasma_subset <- RunPCA(Plasma_subset, verbose = FALSE)

Plasma_subset <- FindNeighbors(Plasma_subset, reduction = "pca", dims = 1:50)
Plasma_subset <- FindClusters(Plasma_subset, resolution = 0.3)
Plasma_subset <- RunUMAP(Plasma_subset, reduction = "pca", dims = 1:50)

DimPlot(Plasma_subset, reduction = 'umap', raster=FALSE, label = T) + NoLegend() + ggtitle("Plasma cells")
FeaturePlot(Plasma_subset, features = 'CD19', min.cutoff = '0', max.cutoff = '3', raster=FALSE)
FeaturePlot(Plasma_subset, features = 'GITR', min.cutoff = '0', max.cutoff = '3', raster=FALSE)
FeaturePlot(Plasma_subset, features = 'HSPB1', min.cutoff = '0', max.cutoff = '3', raster=FALSE)

##################################### Plasma : recur vs non recur #########################################

cancer_Plasma_subset <- subset(x = Plasma_subset, Type == 'Tumor' & Outcome != 'Meta')

DimPlot(cancer_Plasma_subset, reduction = 'umap', group.by = "Outcome", raster=FALSE, label = F) + ggtitle("Plasma cells")
DimPlot(cancer_Plasma_subset, reduction = 'umap', split.by = "Outcome", raster=FALSE, label = T) + NoLegend()+ ggtitle("Plasma cells")


do_BarPlot(sample = cancer_Plasma_subset, 
           group.by = "Outcome",
           split.by = "seurat_clusters",
           plot.title = "Number of cells per sample in each cluster",
           position = "stack")


SCpubr::do_BarPlot(sample = cancer_Plasma_subset, 
                   group.by = "seurat_clusters", 
                   split.by = "Outcome",
                   position = "fill",
                   plot.title = "Fraction of cells by cluster")

##################################### Stats #########################################
##################################### Dataframe for Outcome
md.cancer_Plasma_subset <- cancer_Plasma_subset@meta.data %>% as.data.table

# Vector of sample numbers
sample_numbers <- c(4, 5, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 51, 52, 53, 54, 56, 59, 60, 62, 64, 65, 66, 68, 71, 73, 75)

# Vector of cluster names
cluster_vec <- unique(md.cancer_Plasma_subset$seurat_clusters)
type_vec <- c("Recur", "Non_recur")

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
      subset_length <- length(which(md.cancer_Plasma_subset$Patient_ID_SPC == number & md.cancer_Plasma_subset$Outcome == type & md.cancer_Plasma_subset$seurat_clusters == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.cancer_Plasma_subset$Patient_ID_SPC == number & md.cancer_Plasma_subset$Outcome == type))
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
  Recur_df <- subset_df %>% subset(type_vec == "Recur")
  Non_recur_df <- subset_df %>% subset(type_vec == "Non_recur")
  # Perform the unpaired t-test
  ttest <- t.test(Recur_df$subset_fraction, Non_recur_df$subset_fraction)
  wtest <- wilcox.test(Recur_df$subset_fraction, Non_recur_df$subset_fraction)
  
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
    geom_point(size = 4, shape =21, color = "black")+
    scale_fill_manual(values = c("brown", "grey")) +
    #geom_line(aes(group = sample_numbers), color = "grey") +
    labs(title = paste(cluster),
         x = "Type", y = "Subset Fraction") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          axis.ticks.length = unit(0.1, "cm"),
          #axis.text.y = element_blank(),
          #axis.text.x = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          legend.position = "none")+ 
    
    stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
    stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", size = 0.7, width=0.5, color="#03839C", alpha = 0.5)+
    annotate("text", x = 1.5, y = 0.9, label = paste("t.p =", format(ttest$p.value, digits = 3)))+
    annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[which(cluster_vec == cluster)]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 3)


##################################### Myeloid cell subset #########################################

Myeloid_subset <- subset(x = immune_noannotation, seurat_clusters == '3')

Myeloid_subset <- RunPCA(Myeloid_subset, verbose = FALSE)

Myeloid_subset <- FindNeighbors(Myeloid_subset, reduction = "pca", dims = 1:50)
Myeloid_subset <- FindClusters(Myeloid_subset, resolution = 0.3)
Myeloid_subset <- RunUMAP(Myeloid_subset, reduction = "pca", dims = 1:50)

DimPlot(Myeloid_subset, reduction = 'umap', raster=FALSE, label = T) + NoLegend() + ggtitle("Myeloid cells")
FeaturePlot(Myeloid_subset, features = 'CD14', min.cutoff = '0', max.cutoff = '3', raster=FALSE)
FeaturePlot(Myeloid_subset, features = 'FCGR3A', min.cutoff = '0', max.cutoff = '3', raster=FALSE)
FeaturePlot(Myeloid_subset, features = 'CD86', min.cutoff = '0', max.cutoff = '3', raster=FALSE)
FeaturePlot(Myeloid_subset, features = 'MKI67', min.cutoff = '0', max.cutoff = '3', raster=FALSE)
FeaturePlot(Myeloid_subset, features = 'HSPB1', min.cutoff = '0', max.cutoff = '3', raster=FALSE)

##################################### Myeloid : recur vs non recur #########################################

cancer_Myeloid_subset <- subset(x = Myeloid_subset, Type == 'Tumor' & Outcome != 'Meta')

DimPlot(cancer_Myeloid_subset, reduction = 'umap', group.by = "Outcome", raster=FALSE, label = F) + ggtitle("Myeloid cells")
DimPlot(cancer_Myeloid_subset, reduction = 'umap', split.by = "Outcome", raster=FALSE, label = T) + NoLegend()+ ggtitle("Myeloid cells")


do_BarPlot(sample = cancer_Myeloid_subset, 
           group.by = "Outcome",
           split.by = "seurat_clusters",
           plot.title = "Number of cells per sample in each cluster",
           position = "stack")


SCpubr::do_BarPlot(sample = cancer_Myeloid_subset, 
                   group.by = "seurat_clusters", 
                   split.by = "Outcome",
                   position = "fill",
                   plot.title = "Fraction of cells by cluster")

##################################### Stats #########################################
##################################### Dataframe for Outcome
md.cancer_Myeloid_subset <- cancer_Myeloid_subset@meta.data %>% as.data.table

# Vector of sample numbers
sample_numbers <- c(4, 5, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 51, 52, 53, 54, 56, 59, 60, 62, 64, 65, 66, 68, 71, 73, 75)

# Vector of cluster names
cluster_vec <- unique(md.cancer_Myeloid_subset$seurat_clusters)
type_vec <- c("Recur", "Non_recur")

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
      subset_length <- length(which(md.cancer_Myeloid_subset$Patient_ID_SPC == number & md.cancer_Myeloid_subset$Outcome == type & md.cancer_Myeloid_subset$seurat_clusters == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.cancer_Myeloid_subset$Patient_ID_SPC == number & md.cancer_Myeloid_subset$Outcome == type))
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
  Recur_df <- subset_df %>% subset(type_vec == "Recur")
  Non_recur_df <- subset_df %>% subset(type_vec == "Non_recur")
  # Perform the unpaired t-test
  ttest <- t.test(Recur_df$subset_fraction, Non_recur_df$subset_fraction)
  wtest <- wilcox.test(Recur_df$subset_fraction, Non_recur_df$subset_fraction)
  
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
    geom_point(size = 4, shape =21, color = "black")+
    scale_fill_manual(values = c("brown", "grey")) +
    #geom_line(aes(group = sample_numbers), color = "grey") +
    labs(title = paste(cluster),
         x = "Type", y = "Subset Fraction") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          axis.ticks.length = unit(0.1, "cm"),
          #axis.text.y = element_blank(),
          #axis.text.x = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          legend.position = "none")+ 
    
    stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
    stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", size = 0.7, width=0.5, color="#03839C", alpha = 0.5)+
    annotate("text", x = 1.5, y = 0.9, label = paste("t.p =", format(ttest$p.value, digits = 3)))+
    annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[which(cluster_vec == cluster)]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 3)












