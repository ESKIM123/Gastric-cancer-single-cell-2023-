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
CD8_subset <- readRDS(file="D:/DATA/Gastric cancer/GS_subsets/CD8_subset.rds")
DimPlot(CD8_subset, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Initial")



##################################### wtp53 vs mtp53 #########################################

view(cancer_MSS_subset@meta.data)

cancer_subset_MSIMSS <- subset(x = CD8_subset, Type == 'Tumor' & Patient_ID_SPC != FALSE)
cancer_subset_MSIMSS <- subset(x = cancer_subset_MSIMSS, Stage != 'Stage IV')
cancer_MSS_subset <- subset(x = cancer_subset_MSIMSS, MSI_type == 'MSS')

cancer_MSS_subset[["p53_type"]]  <- cancer_MSS_subset[["Molecular_type"]]
cancer_MSS_subset@meta.data$clusters <- cancer_MSS_subset$seurat_clusters
cancer_MSS_subset$clusters <- cancer_MSS_subset$seurat_clusters


cancer_MSS_subset@meta.data$p53_type[cancer_MSS_subset@meta.data$p53_type == "GS"] <- "wtTP53"
cancer_MSS_subset@meta.data$p53_type[cancer_MSS_subset@meta.data$p53_type == "CIN"] <- "mtTP53"

SCpubr::do_BarPlot(sample = cancer_MSS_subset, 
                   group.by = "clusters", 
                   split.by = "p53_type",
                   position = "Fill",legend.position = "bottom",
                   plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)


cancer_wtTP53_subset <- subset(x = cancer_MSS_subset, p53_type == 'wtTP53')
SCpubr::do_BarPlot(sample = cancer_wtTP53_subset, 
                   group.by = "clusters", 
                   split.by = "Patient_ID_SPC",
                   position = "Fill",legend.position = "bottom",
                   plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)

cancer_mtTP53_subset <- subset(x = cancer_MSS_subset, p53_type == 'mtTP53')
SCpubr::do_BarPlot(sample = cancer_mtTP53_subset, 
                   group.by = "clusters", 
                   split.by = "Patient_ID_SPC",
                   position = "Fill",legend.position = "bottom",
                   plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)


##################################### Stats #########################################
##################################### Dataframe for Outcome
md.immune_annotation<- cancer_MSS_subset@meta.data %>% as.data.table

length(which(md.immune_annotation$Patient_ID==68 & md.immune_annotation$Type=='Tumor'& md.immune_annotation$seurat_clusters==1))

# Vector of sample numbers
sample_numbers <- c(4, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 36, 40, 41, 44, 50, 52, 53, 54, 56, 59, 60, 65, 66, 68, 73, 75)

# Vector of cluster names
cluster_vec <- unique(md.immune_annotation$clusters)
type_vec <- c("wtTP53", "mtTP53")

unique(cancer_MSS_subset@meta.data$Patient_ID_SPC)



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
      subset_length <- length(which(md.immune_annotation$Patient_ID_SPC == number & md.immune_annotation$p53_type == type & md.immune_annotation$clusters == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.immune_annotation$Patient_ID_SPC == number & md.immune_annotation$p53_type == type))
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
  tumor_df <- subset_df %>% subset(type_vec == "mtTP53") 
  normal_df <- subset_df %>% subset(type_vec == "wtTP53")
  # Perform the unpaired t-test
  ttest <- t.test(normal_df$subset_fraction, tumor_df$subset_fraction)
  wtest <- wilcox.test(normal_df$subset_fraction, tumor_df$subset_fraction)
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
    geom_boxplot(color = "black", outlier.shape = NA)+
    scale_fill_manual(values = c("peru", "moccasin")) +
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
grid.arrange(grobs = plot_list, nrow = 2)

############## Stat add

for (cluster in cluster_vec) {
  # Subset the data for each type_vec
  subset_df <- output_df %>% subset(cluster_vec == cluster)
  tumor_df <- subset_df %>% subset(type_vec == "mtTP53") 
  normal_df <- subset_df %>% subset(type_vec == "wtTP53")
  # Perform the unpaired t-test
  ttest <- t.test(normal_df$subset_fraction, tumor_df$subset_fraction)
  wtest <- wilcox.test(normal_df$subset_fraction, tumor_df$subset_fraction)
  
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
    geom_boxplot(color = "black", outlier.shape = NA)+
    scale_fill_manual(values = c("peru", "moccasin")) +
    labs(title = cluster,
         x = NULL, y = NULL) +
    #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          panel.border = element_rect(color = "black", fill = NA, size = 0.2),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          #axis.ticks.length = unit(0.1, "cm"),
          #plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          plot.title = element_text(size = 8, hjust = 0.5), 
          legend.position = "none")+
    
    #stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
    annotate("text", x = 1.5, y = 0.7, label = paste("p =", format(ttest$p.value, digits = 3)))+
    annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[which(cluster_vec == cluster)]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 2)

ggsave("output_graph_MSI_unpaired.png", plot = grid.arrange(grobs = plot_list, nrow = 1), width = 13, height = 5, units = "in")