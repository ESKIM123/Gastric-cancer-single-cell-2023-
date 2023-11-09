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
Immune_subset_anno <- readRDS(file="D:/DATA/Gastric cancer/immune_cell_annotation_final_OLD_TCRtype.rds")
DimPlot(Immune_subset_anno, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Initial")

Immune_TCRtpye <- Immune_subset_anno

cancer_subset <- subset(x = Immune_TCRtpye, Type == 'Tumor' & Patient_ID_SPC != FALSE & Stage != 'Stage IV')
##################################### TCR type comparison #########################################

# Show as UMAP

DimPlot(Immune_TCRtpye, reduction = 'umap', group.by = "TCR_type", raster=FALSE, label = F, pt.size = 1, order=T,
        cols = c('C0_EXP' = '#e69f00', 'C1_EXP' = '#ebe6dd', 'C2_EXP' = '#ebe6dd')) + NoLegend() +ggtitle("Total expanded clones") 
DimPlot(cancer_subset, reduction = 'umap', group.by = "TCR_type", raster=FALSE, label = F, pt.size = 1, order=T,
        cols = c('C0_EXP' = '#e69f00', 'C1_EXP' = '#ebe6dd', 'C2_EXP' = '#ebe6dd')) + NoLegend() +ggtitle("") + NoAxes()

DimPlot(Immune_TCRtpye, reduction = 'umap', group.by = "TCR_type", raster=FALSE, label = F, pt.size = 1, order=T,
        cols = c('C0_EXP' = '#ebe6dd', 'C1_EXP' = '#56b4e9', 'C2_EXP' = '#ebe6dd')) + NoLegend() +ggtitle("Total expanded clones") 
DimPlot(cancer_subset, reduction = 'umap', group.by = "TCR_type", raster=FALSE, label = F, pt.size = 1, order=T,
        cols = c('C0_EXP' = '#ebe6dd', 'C1_EXP' = '#56b4e9', 'C2_EXP' = '#ebe6dd')) + NoLegend() +ggtitle("") + NoAxes()

DimPlot(Immune_TCRtpye, reduction = 'umap', group.by = "TCR_type", raster=FALSE, label = F, pt.size = 1, order=T,
        cols = c('C0_EXP' = '#ebe6dd', 'C1_EXP' = '#ebe6dd', 'C2_EXP' = '#009e73')) + NoLegend() +ggtitle("Total expanded clones") 
DimPlot(cancer_subset, reduction = 'umap', group.by = "TCR_type", raster=FALSE, label = F, pt.size = 1, order=T,
        cols = c('C0_EXP' = '#ebe6dd', 'C1_EXP' = '#ebe6dd', 'C2_EXP' = '#009e73')) + NoLegend() +ggtitle("") + NoAxes()



dittoBarPlot(cancer_subset, "clusters", group.by = "TCR_type", scale = "count", main = "")+ 
  NoLegend()+ xlab(NULL)+ ylab(NULL)+
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))

dittoBarPlot(cancer_subset, "clusters", group.by = "TCR_type", scale = "percent", main = "")+ 
  NoLegend()+ xlab(NULL)+ ylab(NULL)+
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))


dittoBarPlot(cancer_subset, "clusters", group.by = "Patient_ID_SPC", scale = "percent", main = "")+ 
  NoLegend()+ xlab(NULL)+ ylab(NULL)+
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))

##################################### Stats #########################################
##################################### Dataframe for Outcome
md.immune_annotation<- cancer_subset@meta.data %>% as.data.table

length(which(md.immune_annotation$Patient_ID==68 & md.immune_annotation$Type=='Tumor'& md.immune_annotation$clusters==1))

# Vector of sample numbers
sample_numbers <- unique(cancer_subset@meta.data$Patient_ID_SPC)

# Vector of cluster names
cluster_vec <- unique(md.immune_annotation$clusters)
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
      subset_length <- length(which(md.immune_annotation$Patient_ID_SPC == number & md.immune_annotation$TCR_type == type & md.immune_annotation$clusters == cluster))
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
          #panel.border = element_rect(color = "black", fill = NA, size = 0.2),
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
grid.arrange(grobs = plot_list, nrow = 2)

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
grid.arrange(grobs = plot_list, nrow = 2)



############## do_AlluvialPlot ####



do_AlluvialPlot(sample = cancer_subset,
                first_group = "TCR_type",
                last_group = "clusters")
