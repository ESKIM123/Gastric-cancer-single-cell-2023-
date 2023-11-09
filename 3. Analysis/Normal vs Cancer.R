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
cancer.subset <- readRDS(file="D:/DATA/Gastric cancer/immune_cell_clinical_data_modify.rds")

immune_annotation <- subset(x = cancer.subset, Patient_ID_SPC != FALSE & Patient_ID_SPC != 8 & Patient_ID_SPC != 18 & Patient_ID_SPC != 52)

DimPlot(immune_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Initial")

##################################### Normal vs Tumors #########################################
DimPlot(immune_annotation, reduction = 'umap', group.by = "Type", raster=FALSE, label = F, cols = c("#6495ED", "lightgrey"))+
  NoAxes()+NoLegend()+ ggtitle("")



SCpubr::do_BarPlot(sample = immune_annotation, 
                   group.by = "clusters", 
                   split.by = "Type",
                   position = "stack",legend.position = "bottom",
                   plot.title = "", xlab = "",  ylab = "Number of cells per cluster", font.size = 11)

SCpubr::do_BarPlot(sample = immune_annotation, 
                   group.by = "clusters", 
                   split.by = "Type",
                   position = "Fill",legend.position = "bottom",
                   plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)

dittoBarPlot(immune_annotation, "clusters", group.by = "Type", scale = "percent", main = "")+ 
  NoLegend()+ xlab(NULL)+ ylab(NULL)+
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))


dittoBarPlot(immune_annotation, "clusters", group.by = "Type", scale = "count", main = "")+ 
  NoLegend()+ xlab(NULL)+ ylab(NULL)+
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))


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


dittoBarPlot(cancer_subset, "clusters", group.by = "Patient_ID_SPC", scale = "percent", main = "")+ 
  NoLegend()+ xlab(NULL)+ ylab(NULL)+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))


dittoBarPlot(Normal_subset, "clusters", group.by = "Patient_ID_SPC", scale = "percent", main = "")+ 
  NoLegend()+ xlab(NULL)+ ylab(NULL)+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))




##################################### Stats PAIRED ★ ######################################### 
##################################### Dataframe for Outcome
md.immune_annotation<- immune_annotation@meta.data %>% as.data.table

length(which(md.immune_annotation$Patient_ID==68 & md.immune_annotation$Type=='Tumor'& md.immune_annotation$seurat_clusters==1))

# Vector of sample numbers
sample_numbers <- c(15, 22, 28, 29, 36, 41, 5, 54, 56, 59, 64, 65, 66, 68, 7, 71, 73)

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
grid.arrange(grobs = plot_list, nrow = 3)






############## Stat add

for (cluster in cluster_vec) {
  # Subset the data for each type_vec
  subset_df <- output_df %>% subset(cluster_vec == cluster)
  normal_df <- subset_df %>% subset(type_vec == "Normal")
  tumor_df <- subset_df %>% subset(type_vec == "Tumor")
  # Perform the unpaired t-test
  ttest <- t.test(normal_df$subset_fraction, tumor_df$subset_fraction)
  wtest <- wilcox.test(normal_df$subset_fraction, tumor_df$subset_fraction)
  
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction)) +
    geom_boxplot(color = "black", outlier.shape = NA)+
    scale_fill_manual(values = c("grey", "brown")) +
    geom_point(aes(color = type_vec, shape = factor(type_vec))) +
    geom_line(aes(group = sample_numbers))+
    ggtitle(cluster) +
    xlab("") + ylab("") +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          axis.ticks.length = unit(0.1, "cm"),
          axis.text.y = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          legend.position = "none")+
    
    #stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
    annotate("text", x = 1.5, y = 0.7, label = paste("p =", format(ttest$p.value, digits = 3)))+
    annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[which(cluster_vec == cluster)]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 3)




##################################### Stats MSI ★ ######################################### 
##################################### Dataframe for Outcome
md.immune_annotation<- immune_annotation@meta.data %>% as.data.table

length(which(md.immune_annotation$Patient_ID==68 & md.immune_annotation$Type=='Tumor'& md.immune_annotation$seurat_clusters==1))

# Vector of sample numbers
sample_numbers <- c(71)

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
  # ttest <- t.test(normal_df$subset_fraction, tumor_df$subset_fraction)
  # wtest <- wilcox.test(normal_df$subset_fraction, tumor_df$subset_fraction)
  
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
grid.arrange(grobs = plot_list, nrow = 3)


# ############## Stat add
# 
# for (cluster in cluster_vec) {
#   # Subset the data for each type_vec
#   subset_df <- output_df %>% subset(cluster_vec == cluster)
#   normal_df <- subset_df %>% subset(type_vec == "Normal")
#   tumor_df <- subset_df %>% subset(type_vec == "Tumor")
#   # Perform the unpaired t-test
#   ttest <- t.test(normal_df$subset_fraction, tumor_df$subset_fraction)
#   wtest <- wilcox.test(normal_df$subset_fraction, tumor_df$subset_fraction)
#   
#   g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction)) +
#     geom_boxplot(color = "black", outlier.shape = NA)+
#     scale_fill_manual(values = c("grey", "brown")) +
#     geom_point(aes(color = type_vec, shape = factor(type_vec))) +
#     geom_line(aes(group = sample_numbers))+
#     ggtitle(cluster) +
#     xlab("") + ylab("") +
#     theme(panel.background = element_rect(fill = "#FFFFFF"),
#           axis.line = element_line(color = "black", size = 0.2),
#           axis.ticks = element_line(size = 0.2),
#           axis.ticks.length = unit(0.1, "cm"),
#           axis.text.y = element_blank(),
#           plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
#           legend.position = "none")+
#     
#     #stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
#     annotate("text", x = 1.5, y = 0.7, label = paste("p =", format(ttest$p.value, digits = 3)))+
#     annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
#   
#   # Add the plot to the list
#   plot_list[[which(cluster_vec == cluster)]] <- g1
# }
# 
# # Arrange the plots in a grid
# grid.arrange(grobs = plot_list, nrow = 3)



##################################### Stats CIN ★ ######################################### 
##################################### Dataframe for Outcome
md.immune_annotation<- immune_annotation@meta.data %>% as.data.table

length(which(md.immune_annotation$Patient_ID==68 & md.immune_annotation$Type=='Tumor'& md.immune_annotation$seurat_clusters==1))

# Vector of sample numbers
sample_numbers <- c(15, 22, 29, 36, 65, 66, 73)

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
grid.arrange(grobs = plot_list, nrow = 3)



############## Stat add

for (cluster in cluster_vec) {
  # Subset the data for each type_vec
  subset_df <- output_df %>% subset(cluster_vec == cluster)
  normal_df <- subset_df %>% subset(type_vec == "Normal")
  tumor_df <- subset_df %>% subset(type_vec == "Tumor")
  # Perform the unpaired t-test
  ttest <- t.test(normal_df$subset_fraction, tumor_df$subset_fraction)
  wtest <- wilcox.test(normal_df$subset_fraction, tumor_df$subset_fraction)

  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction)) +
    geom_boxplot(color = "black", outlier.shape = NA)+
    scale_fill_manual(values = c("grey", "brown")) +
    geom_point(aes(color = type_vec, shape = factor(type_vec))) +
    geom_line(aes(group = sample_numbers))+
    ggtitle(cluster) +
    xlab("") + ylab("") +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          axis.ticks.length = unit(0.1, "cm"),
          axis.text.y = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          legend.position = "none")+

    #stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
    annotate("text", x = 1.5, y = 0.7, label = paste("p =", format(ttest$p.value, digits = 3)))+
    annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))

  # Add the plot to the list
  plot_list[[which(cluster_vec == cluster)]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 3)




##################################### Stats GS ★ ######################################### 
##################################### Dataframe for Outcome
md.immune_annotation<- immune_annotation@meta.data %>% as.data.table

length(which(md.immune_annotation$Patient_ID==68 & md.immune_annotation$Type=='Tumor'& md.immune_annotation$seurat_clusters==1))

# Vector of sample numbers
sample_numbers <- c(5, 7, 28, 41, 54, 56, 59, 64, 68)

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
grid.arrange(grobs = plot_list, nrow = 3)



############## Stat add

for (cluster in cluster_vec) {
  # Subset the data for each type_vec
  subset_df <- output_df %>% subset(cluster_vec == cluster)
  normal_df <- subset_df %>% subset(type_vec == "Normal")
  tumor_df <- subset_df %>% subset(type_vec == "Tumor")
  # Perform the unpaired t-test
  ttest <- t.test(normal_df$subset_fraction, tumor_df$subset_fraction)
  wtest <- wilcox.test(normal_df$subset_fraction, tumor_df$subset_fraction)
  
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction)) +
    geom_boxplot(color = "black", outlier.shape = NA)+
    scale_fill_manual(values = c("grey", "brown")) +
    geom_point(aes(color = type_vec, shape = factor(type_vec))) +
    geom_line(aes(group = sample_numbers))+
    ggtitle(cluster) +
    xlab("") + ylab("") +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          axis.ticks.length = unit(0.1, "cm"),
          axis.text.y = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          legend.position = "none")+
    
    #stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
    annotate("text", x = 1.5, y = 0.7, label = paste("p =", format(ttest$p.value, digits = 3)))+
    annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
  
  # Add the plot to the list
  plot_list[[which(cluster_vec == cluster)]] <- g1
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 3)















# 
# ggsave("output_graph_paired.png", plot = grid.arrange(grobs = plot_list, nrow = 3), width = 13, height = 5, units = "in")
# 
# 
# 
# 
# 
# 
# ############## Stat add
# 
# for (cluster in cluster_vec) {
#   # Subset the data for each type_vec
#   subset_df <- output_df %>% subset(cluster_vec == cluster)
#   normal_df <- subset_df %>% subset(type_vec == "Normal")
#   tumor_df <- subset_df %>% subset(type_vec == "Tumor")
#   # Perform the unpaired t-test
#   ttest <- t.test(normal_df$subset_fraction, tumor_df$subset_fraction)
#   wtest <- wilcox.test(normal_df$subset_fraction, tumor_df$subset_fraction)
#   
#   g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction)) +
#     geom_boxplot(color = "black", outlier.shape = NA, fill = c("grey", "brown"), alpha=0.5)+
#     geom_point(shape = 21, color = "black", alpha=0.8, aes(fill=type_vec)) + 
#     scale_fill_manual(values = c("grey", "brown")) +
#     geom_line(aes(group = sample_numbers), alpha=0.5)+
#     ggtitle(cluster) + 
#     xlab("") + ylab("") +
#     theme(panel.background = element_rect(fill = "#FFFFFF"),
#           axis.line = element_line(color = "black", size = 0.2),
#           axis.ticks = element_line(size = 0.2),
#           axis.ticks.length = unit(0.1, "cm"),
#           plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
#           legend.position = "none")+
#     
#     #stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
#     annotate("text", x = 1.5, y = 0.7, label = paste("p =", format(ttest$p.value, digits = 3)))+
#     annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
#   
#   # Add the plot to the list
#   plot_list[[which(cluster_vec == cluster)]] <- g1
# }
# 
# # Arrange the plots in a grid
# grid.arrange(grobs = plot_list, nrow = 3)
# 
# 
# # 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# immune_annotation_notMETA <- subset(x = immune_annotation, Stage == 'Stage II' | Stage == 'Stage III')
# view(immune_annotation_notMETA@meta.data)
# 
# SCpubr::do_BarPlot(sample = immune_annotation_notMETA, 
#                    group.by = "clusters", 
#                    split.by = "Type",
#                    position = "stack",legend.position = "bottom",
#                    plot.title = "", xlab = "",  ylab = "Number of cells per cluster", font.size = 11)
# 
# SCpubr::do_BarPlot(sample = immune_annotation_notMETA, 
#                    group.by = "clusters", 
#                    split.by = "Type",
#                    position = "Fill",legend.position = "bottom",
#                    plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)
# 
# 
# cancer_subset_stage <- subset(x = immune_annotation_notMETA, Type == 'Tumor' & Patient_ID_SPC != FALSE)
# SCpubr::do_BarPlot(sample = cancer_subset_stage, 
#                    group.by = "clusters", 
#                    split.by = "Patient_ID_SPC",
#                    position = "Fill",legend.position = "bottom",
#                    plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)
# 
# 
# dittoBarPlot(cancer_subset_stage, "clusters", group.by = "Patient_ID_SPC", scale = "percent", main = "")+ 
#   NoLegend()+ xlab(NULL)+ ylab(NULL)+
#   scale_y_continuous(expand = c(0, 0))+
#   theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))
# 
# 
# Normal_subset_stage <- subset(x = immune_annotation_notMETA, Type == 'Normal' & Patient_ID_SPC != FALSE)
# SCpubr::do_BarPlot(sample = Normal_subset_stage, 
#                    group.by = "clusters", 
#                    split.by = "Patient_ID_SPC",
#                    position = "Fill",legend.position = "bottom",
#                    plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)
# 
# dittoBarPlot(Normal_subset_stage, "clusters", group.by = "Patient_ID_SPC", scale = "percent", main = "")+ 
#   NoLegend()+ xlab(NULL)+ ylab(NULL)+
#   scale_y_continuous(expand = c(0, 0))+
#   theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))
# 
# dittoBarPlot(Normal_subset_stage, "clusters", group.by = "Patient_ID_SPC", scale = "percent", main = "")
# 
# ##################################### Stats #########################################
# ##################################### Dataframe for Outcome
# md.immune_annotation<- immune_annotation@meta.data %>% as.data.table
# 
# length(which(md.immune_annotation$Patient_ID==68 & md.immune_annotation$Type=='Tumor'& md.immune_annotation$seurat_clusters==1))
# 
# # Vector of sample numbers
# sample_numbers <- c(4, 5, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 51, 52, 53, 54, 56, 59, 60, 62, 64, 65, 66, 68, 71, 73, 75)
# 
# # Vector of cluster names
# cluster_vec <- unique(md.immune_annotation$clusters)
# type_vec <- c("Tumor", "Normal")
# 
# unique(cancer_subset@meta.data$Patient_ID_SPC)
# 
# 
# 
# # Create an empty data frame to store the results
# output_df <- data.frame(sample_numbers = numeric(),
#                         cluster_vec = character(),
#                         type_vec = character(),
#                         count = numeric(),
#                         total_cell = numeric(),
#                         subset_fraction = numeric(),
#                         stringsAsFactors = FALSE)
# 
# # Loop over the variables and add rows to the data frame
# for (type in type_vec) {
#   for (cluster in cluster_vec) {
#     for (number in sample_numbers) {
#       # Calculate the length of the subset
#       subset_length <- length(which(md.immune_annotation$Patient_ID_SPC == number & md.immune_annotation$Type == type & md.immune_annotation$clusters == cluster))
#       # Calculate the subset fraction
#       total_count <- length(which(md.immune_annotation$Patient_ID_SPC == number & md.immune_annotation$Type == type))
#       subset_fraction <- subset_length / total_count
#       # Add a row to the data frame
#       new_row <- data.frame(sample_numbers = number,
#                             cluster_vec = cluster,
#                             type_vec = type,
#                             count = subset_length,
#                             total_cell = total_count,
#                             subset_fraction = subset_fraction,
#                             stringsAsFactors = FALSE)
#       output_df <- rbind(output_df, new_row)
#     }
#   }
# }
# 
# 
# output_df
# 
# ################################################################################## Graph
# 
# # Create an empty list to store the plots
# plot_list <- list()
# 
# for (cluster in cluster_vec) {
#   # Subset the data for each type_vec
#   subset_df <- output_df %>% subset(cluster_vec == cluster)
#   normal_df <- subset_df %>% subset(type_vec == "Normal")
#   tumor_df <- subset_df %>% subset(type_vec == "Tumor")
#   # Perform the unpaired t-test
#   ttest <- t.test(normal_df$subset_fraction, tumor_df$subset_fraction)
#   wtest <- wilcox.test(normal_df$subset_fraction, tumor_df$subset_fraction)
#   g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
#     geom_boxplot(color = "black", outlier.shape = NA)+
#     scale_fill_manual(values = c("grey", "brown")) +
#     labs(title = NULL,
#          x = NULL, y = NULL) +
#     #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
#     theme(panel.background = element_rect(fill = "#FFFFFF"),
#           panel.border = element_rect(color = "black", fill = NA, size = 0.2),
#           axis.line = element_line(color = "black", size = 0.2),
#           axis.ticks = element_line(size = 0.2),
#           #axis.ticks.length = unit(0.1, "cm"),
#           #plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
#           plot.title = element_text(size = 8, hjust = 0.5), 
#           legend.position = "none")
#   
#   #stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")
#   #annotate("text", x = 1.5, y = 0.5, label = paste("p =", format(ttest$p.value, digits = 3)))
#   #annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
#   
#   # Add the plot to the list
#   plot_list[[which(cluster_vec == cluster)]] <- g1
# }
# 
# # Arrange the plots in a grid
# grid.arrange(grobs = plot_list, nrow = 1)
# 
# ############## Stat add
# 
# for (cluster in cluster_vec) {
#   # Subset the data for each type_vec
#   subset_df <- output_df %>% subset(cluster_vec == cluster)
#   normal_df <- subset_df %>% subset(type_vec == "Normal")
#   tumor_df <- subset_df %>% subset(type_vec == "Tumor")
#   # Perform the unpaired t-test
#   ttest <- t.test(normal_df$subset_fraction, tumor_df$subset_fraction)
#   wtest <- wilcox.test(normal_df$subset_fraction, tumor_df$subset_fraction)
#   
#   g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
#     geom_boxplot(color = "black", outlier.shape = NA)+
#     scale_fill_manual(values = c("grey", "brown")) +
#     labs(title = cluster,
#          x = NULL, y = NULL) +
#     #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
#     theme(panel.background = element_rect(fill = "#FFFFFF"),
#           panel.border = element_rect(color = "black", fill = NA, size = 0.2),
#           axis.line = element_line(color = "black", size = 0.2),
#           axis.ticks = element_line(size = 0.2),
#           #axis.ticks.length = unit(0.1, "cm"),
#           #plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
#           plot.title = element_text(size = 8, hjust = 0.5), 
#           legend.position = "none")+
#   
#   #stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
#   annotate("text", x = 1.5, y = 0.7, label = paste("p =", format(ttest$p.value, digits = 3)))+
#   annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
#   
#   # Add the plot to the list
#   plot_list[[which(cluster_vec == cluster)]] <- g1
# }
# 
# # Arrange the plots in a grid
# grid.arrange(grobs = plot_list, nrow = 1)
# 
# ggsave("output_graph_unpaired.png", plot = grid.arrange(grobs = plot_list, nrow = 1), width = 13, height = 5, units = "in")
# 
# 
# 
# ##################################### Stats PAIRED ★ ######################################### 
# ##################################### Dataframe for Outcome
# md.immune_annotation<- immune_annotation@meta.data %>% as.data.table
# 
# length(which(md.immune_annotation$Patient_ID==68 & md.immune_annotation$Type=='Tumor'& md.immune_annotation$seurat_clusters==1))
# 
# # Vector of sample numbers
# sample_numbers <- c(15, 18, 22, 28, 29, 36, 41, 5, 54, 56, 59, 64, 65, 66, 68, 7, 71, 73, 8)
# 
# # Vector of cluster names
# cluster_vec <- unique(md.immune_annotation$clusters)
# type_vec <- c("Tumor", "Normal")
# 
# unique(cancer_subset@meta.data$Patient_ID_SPC)
# 
# 
# 
# # Create an empty data frame to store the results
# output_df <- data.frame(sample_numbers = numeric(),
#                         cluster_vec = character(),
#                         type_vec = character(),
#                         count = numeric(),
#                         total_cell = numeric(),
#                         subset_fraction = numeric(),
#                         stringsAsFactors = FALSE)
# 
# # Loop over the variables and add rows to the data frame
# for (type in type_vec) {
#   for (cluster in cluster_vec) {
#     for (number in sample_numbers) {
#       # Calculate the length of the subset
#       subset_length <- length(which(md.immune_annotation$Patient_ID_SPC == number & md.immune_annotation$Type == type & md.immune_annotation$clusters == cluster))
#       # Calculate the subset fraction
#       total_count <- length(which(md.immune_annotation$Patient_ID_SPC == number & md.immune_annotation$Type == type))
#       subset_fraction <- subset_length / total_count
#       # Add a row to the data frame
#       new_row <- data.frame(sample_numbers = number,
#                             cluster_vec = cluster,
#                             type_vec = type,
#                             count = subset_length,
#                             total_cell = total_count,
#                             subset_fraction = subset_fraction,
#                             stringsAsFactors = FALSE)
#       output_df <- rbind(output_df, new_row)
#     }
#   }
# }
# 
# 
# output_df
# 
# ################################################################################## Graph
# 
# # Create an empty list to store the plots
# plot_list <- list()
# 
# for (cluster in cluster_vec) {
#   # Subset the data for each type_vec
#   subset_df <- output_df %>% subset(cluster_vec == cluster)
#   normal_df <- subset_df %>% subset(type_vec == "Normal")
#   tumor_df <- subset_df %>% subset(type_vec == "Tumor")
#   # Perform the unpaired t-test
#   ttest <- t.test(normal_df$subset_fraction, tumor_df$subset_fraction)
#   wtest <- wilcox.test(normal_df$subset_fraction, tumor_df$subset_fraction)
#   
#   g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction)) +
#     geom_boxplot(color = "black", outlier.shape = NA, fill = c("#E1D4BB", "#537188"), alpha=0.8)+
#     geom_point(shape = 21, color = "black", alpha=0.8, aes(fill=type_vec)) + 
#     scale_fill_manual(values = c("grey", "brown")) +
#     scale_y_continuous(labels = function(x) paste0(x * 100))+
#     geom_line(aes(group = sample_numbers), alpha=0.5)+
#     ggtitle("") + 
#     xlab("") + ylab("") +
#     theme(panel.background = element_rect(fill = "#FFFFFF"),
#           panel.border = element_rect(color = "black", fill = NA, size = 0.2),
#           axis.line = element_line(color = "black", size = 0.2),
#           axis.ticks = element_line(size = 0.2),
#           axis.ticks.length = unit(0.1, "cm"),
#           plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
#           legend.position = "none")
#     
#     #stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
#     #annotate("text", x = 1.5, y = 0.7, label = paste("p =", format(ttest$p.value, digits = 3)))+
#     #annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
#   
#   # Add the plot to the list
#   plot_list[[which(cluster_vec == cluster)]] <- g1
# }
# 
# # Arrange the plots in a grid
# grid.arrange(grobs = plot_list, nrow = 3)
# 
# 
# ############## Stat add
# 
# for (cluster in cluster_vec) {
#   # Subset the data for each type_vec
#   subset_df <- output_df %>% subset(cluster_vec == cluster)
#   normal_df <- subset_df %>% subset(type_vec == "Normal")
#   tumor_df <- subset_df %>% subset(type_vec == "Tumor")
#   # Perform the unpaired t-test
#   ttest <- t.test(normal_df$subset_fraction, tumor_df$subset_fraction)
#   wtest <- wilcox.test(normal_df$subset_fraction, tumor_df$subset_fraction)
#   
#   g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction)) +
#     geom_boxplot(color = "black", outlier.shape = NA, fill = c("grey", "brown"), alpha=0.5)+
#     geom_point(shape = 21, color = "black", alpha=0.8, aes(fill=type_vec)) + 
#     scale_fill_manual(values = c("grey", "brown")) +
#     geom_line(aes(group = sample_numbers), alpha=0.5)+
#     ggtitle(cluster) + 
#     xlab("") + ylab("") +
#     theme(panel.background = element_rect(fill = "#FFFFFF"),
#           axis.line = element_line(color = "black", size = 0.2),
#           axis.ticks = element_line(size = 0.2),
#           axis.ticks.length = unit(0.1, "cm"),
#           plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
#           legend.position = "none")+
#   
#   #stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
#   annotate("text", x = 1.5, y = 0.7, label = paste("p =", format(ttest$p.value, digits = 3)))+
#   annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
#   
#   # Add the plot to the list
#   plot_list[[which(cluster_vec == cluster)]] <- g1
# }
# 
# # Arrange the plots in a grid
# grid.arrange(grobs = plot_list, nrow = 3)
# 
# 
# 
# ############## Stat add
# 
# for (cluster in cluster_vec) {
#   # Subset the data for each type_vec
#   subset_df <- output_df %>% subset(cluster_vec == cluster)
#   normal_df <- subset_df %>% subset(type_vec == "Normal")
#   tumor_df <- subset_df %>% subset(type_vec == "Tumor")
#   # Perform the unpaired t-test
#   ttest <- t.test(normal_df$subset_fraction, tumor_df$subset_fraction)
#   wtest <- wilcox.test(normal_df$subset_fraction, tumor_df$subset_fraction)
#   
#   g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction)) +
#     geom_boxplot(color = "black", outlier.shape = NA)+
#     scale_fill_manual(values = c("grey", "brown")) +
#     geom_point(aes(color = type_vec, shape = factor(type_vec))) + 
#     geom_line(aes(group = sample_numbers))+
#     ggtitle(cluster) + 
#     xlab("") + ylab("") +
#     theme(panel.background = element_rect(fill = "#FFFFFF"),
#           axis.line = element_line(color = "black", size = 0.2),
#           axis.ticks = element_line(size = 0.2),
#           axis.ticks.length = unit(0.1, "cm"),
#           axis.text.y = element_blank(),
#           plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
#           legend.position = "none")+
#     
#     #stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
#     annotate("text", x = 1.5, y = 0.7, label = paste("p =", format(ttest$p.value, digits = 3)))+
#     annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
#   
#   # Add the plot to the list
#   plot_list[[which(cluster_vec == cluster)]] <- g1
# }
# 
# # Arrange the plots in a grid
# grid.arrange(grobs = plot_list, nrow = 1)
# 
# 
# 
# ggsave("output_graph_paired.png", plot = grid.arrange(grobs = plot_list, nrow = 1), width = 13, height = 5, units = "in")
# 
# 
# ##################################### Stats: Stage #########################################
# ##################################### Dataframe for Outcome
# md.immune_annotation<- immune_annotation@meta.data %>% as.data.table
# 
# length(which(md.immune_annotation$Patient_ID==68 & md.immune_annotation$Type=='Tumor'& md.immune_annotation$seurat_clusters==1))
# 
# # Vector of sample numbers
# sample_numbers <- c(4, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 52, 53, 54, 56, 59, 60, 65, 66, 68, 71, 73, 75)
# 
# # Vector of cluster names
# cluster_vec <- unique(md.immune_annotation$clusters)
# type_vec <- c("Tumor", "Normal")
# 
# unique(cancer_subset@meta.data$Patient_ID_SPC)
# 
# 
# 
# # Create an empty data frame to store the results
# output_df <- data.frame(sample_numbers = numeric(),
#                         cluster_vec = character(),
#                         type_vec = character(),
#                         count = numeric(),
#                         total_cell = numeric(),
#                         subset_fraction = numeric(),
#                         stringsAsFactors = FALSE)
# 
# # Loop over the variables and add rows to the data frame
# for (type in type_vec) {
#   for (cluster in cluster_vec) {
#     for (number in sample_numbers) {
#       # Calculate the length of the subset
#       subset_length <- length(which(md.immune_annotation$Patient_ID_SPC == number & md.immune_annotation$Type == type & md.immune_annotation$clusters == cluster))
#       # Calculate the subset fraction
#       total_count <- length(which(md.immune_annotation$Patient_ID_SPC == number & md.immune_annotation$Type == type))
#       subset_fraction <- subset_length / total_count
#       # Add a row to the data frame
#       new_row <- data.frame(sample_numbers = number,
#                             cluster_vec = cluster,
#                             type_vec = type,
#                             count = subset_length,
#                             total_cell = total_count,
#                             subset_fraction = subset_fraction,
#                             stringsAsFactors = FALSE)
#       output_df <- rbind(output_df, new_row)
#     }
#   }
# }
# 
# 
# output_df
# 
# ################################################################################## Graph
# 
# # Create an empty list to store the plots
# plot_list <- list()
# 
# for (cluster in cluster_vec) {
#   # Subset the data for each type_vec
#   subset_df <- output_df %>% subset(cluster_vec == cluster)
#   normal_df <- subset_df %>% subset(type_vec == "Normal")
#   tumor_df <- subset_df %>% subset(type_vec == "Tumor")
#   # Perform the unpaired t-test
#   ttest <- t.test(normal_df$subset_fraction, tumor_df$subset_fraction)
#   wtest <- wilcox.test(normal_df$subset_fraction, tumor_df$subset_fraction)
#   g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
#     geom_boxplot(color = "black", outlier.shape = NA)+
#     scale_fill_manual(values = c("brown", "grey")) +
#     labs(title = NULL,
#          x = NULL, y = NULL) +
#     #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
#     theme(panel.background = element_rect(fill = "#FFFFFF"),
#           panel.border = element_rect(color = "black", fill = NA, size = 0.2),
#           axis.line = element_line(color = "black", size = 0.2),
#           axis.ticks = element_line(size = 0.2),
#           #axis.ticks.length = unit(0.1, "cm"),
#           #plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
#           plot.title = element_text(size = 8, hjust = 0.5), 
#           legend.position = "none")
#   
#   #stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")
#   #annotate("text", x = 1.5, y = 0.5, label = paste("p =", format(ttest$p.value, digits = 3)))
#   #annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
#   
#   # Add the plot to the list
#   plot_list[[which(cluster_vec == cluster)]] <- g1
# }
# 
# # Arrange the plots in a grid
# grid.arrange(grobs = plot_list, nrow = 1)
# 
# ############## Stat add
# 
# for (cluster in cluster_vec) {
#   # Subset the data for each type_vec
#   subset_df <- output_df %>% subset(cluster_vec == cluster)
#   normal_df <- subset_df %>% subset(type_vec == "Normal")
#   tumor_df <- subset_df %>% subset(type_vec == "Tumor")
#   # Perform the unpaired t-test
#   ttest <- t.test(normal_df$subset_fraction, tumor_df$subset_fraction)
#   wtest <- wilcox.test(normal_df$subset_fraction, tumor_df$subset_fraction)
#   
#   g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
#     geom_boxplot(color = "black", outlier.shape = NA)+
#     scale_fill_manual(values = c("brown", "grey")) +
#     labs(title = cluster,
#          x = NULL, y = NULL) +
#     #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
#     theme(panel.background = element_rect(fill = "#FFFFFF"),
#           panel.border = element_rect(color = "black", fill = NA, size = 0.2),
#           axis.line = element_line(color = "black", size = 0.2),
#           axis.ticks = element_line(size = 0.2),
#           #axis.ticks.length = unit(0.1, "cm"),
#           #plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
#           plot.title = element_text(size = 8, hjust = 0.5), 
#           legend.position = "none")+
#     
#     #stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
#     annotate("text", x = 1.5, y = 0.7, label = paste("p =", format(ttest$p.value, digits = 3)))+
#     annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
#   
#   # Add the plot to the list
#   plot_list[[which(cluster_vec == cluster)]] <- g1
# }
# 
# # Arrange the plots in a grid
# grid.arrange(grobs = plot_list, nrow = 1)
# 
# ggsave("output_graph_unpaired.png", plot = grid.arrange(grobs = plot_list, nrow = 1), width = 13, height = 5, units = "in")
# 
# 
# 
# ##################################### Stats ######################################### PAIRED
# ##################################### Dataframe for Outcome
# md.immune_annotation<- immune_annotation@meta.data %>% as.data.table
# 
# length(which(md.immune_annotation$Patient_ID==68 & md.immune_annotation$Type=='Tumor'& md.immune_annotation$seurat_clusters==1))
# 
# # Vector of sample numbers
# sample_numbers <- c(15, 18, 22, 28, 29, 36, 41, 54, 56, 59, 65, 66, 68, 7, 71, 73, 8)
# 
# # Vector of cluster names
# cluster_vec <- unique(md.immune_annotation$clusters)
# type_vec <- c("Tumor", "Normal")
# 
# unique(cancer_subset@meta.data$Patient_ID_SPC)
# 
# 
# 
# # Create an empty data frame to store the results
# output_df <- data.frame(sample_numbers = numeric(),
#                         cluster_vec = character(),
#                         type_vec = character(),
#                         count = numeric(),
#                         total_cell = numeric(),
#                         subset_fraction = numeric(),
#                         stringsAsFactors = FALSE)
# 
# # Loop over the variables and add rows to the data frame
# for (type in type_vec) {
#   for (cluster in cluster_vec) {
#     for (number in sample_numbers) {
#       # Calculate the length of the subset
#       subset_length <- length(which(md.immune_annotation$Patient_ID_SPC == number & md.immune_annotation$Type == type & md.immune_annotation$clusters == cluster))
#       # Calculate the subset fraction
#       total_count <- length(which(md.immune_annotation$Patient_ID_SPC == number & md.immune_annotation$Type == type))
#       subset_fraction <- subset_length / total_count
#       # Add a row to the data frame
#       new_row <- data.frame(sample_numbers = number,
#                             cluster_vec = cluster,
#                             type_vec = type,
#                             count = subset_length,
#                             total_cell = total_count,
#                             subset_fraction = subset_fraction,
#                             stringsAsFactors = FALSE)
#       output_df <- rbind(output_df, new_row)
#     }
#   }
# }
# 
# 
# output_df
# 
# ################################################################################## Graph
# 
# # Create an empty list to store the plots
# plot_list <- list()
# 
# for (cluster in cluster_vec) {
#   # Subset the data for each type_vec
#   subset_df <- output_df %>% subset(cluster_vec == cluster)
#   normal_df <- subset_df %>% subset(type_vec == "Normal")
#   tumor_df <- subset_df %>% subset(type_vec == "Tumor")
#   # Perform the unpaired t-test
#   ttest <- t.test(normal_df$subset_fraction, tumor_df$subset_fraction, paired = TRUE)
#   wtest <- wilcox.test(normal_df$subset_fraction, tumor_df$subset_fraction, paired = TRUE)
#   g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
#     geom_boxplot(color = "black", outlier.shape = NA)+
#     scale_fill_manual(values = c("brown", "grey")) +
#     labs(title = NULL,
#          x = NULL, y = NULL) +
#     #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
#     theme(panel.background = element_rect(fill = "#FFFFFF"),
#           panel.border = element_rect(color = "black", fill = NA, size = 0.2),
#           axis.line = element_line(color = "black", size = 0.2),
#           axis.ticks = element_line(size = 0.2),
#           #axis.ticks.length = unit(0.1, "cm"),
#           #plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
#           plot.title = element_text(size = 8, hjust = 0.5), 
#           legend.position = "none")
#   
#   #stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")
#   #annotate("text", x = 1.5, y = 0.5, label = paste("p =", format(ttest$p.value, digits = 3)))
#   #annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
#   
#   # Add the plot to the list
#   plot_list[[which(cluster_vec == cluster)]] <- g1
# }
# 
# # Arrange the plots in a grid
# grid.arrange(grobs = plot_list, nrow = 1)
# 
# ############## Stat add
# 
# for (cluster in cluster_vec) {
#   # Subset the data for each type_vec
#   subset_df <- output_df %>% subset(cluster_vec == cluster)
#   normal_df <- subset_df %>% subset(type_vec == "Normal")
#   tumor_df <- subset_df %>% subset(type_vec == "Tumor")
#   # Perform the unpaired t-test
#   ttest <- t.test(normal_df$subset_fraction, tumor_df$subset_fraction)
#   wtest <- wilcox.test(normal_df$subset_fraction, tumor_df$subset_fraction)
#   
#   g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
#     geom_boxplot(color = "black", outlier.shape = NA)+
#     scale_fill_manual(values = c("brown", "grey")) +
#     labs(title = cluster,
#          x = NULL, y = NULL) +
#     #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
#     theme(panel.background = element_rect(fill = "#FFFFFF"),
#           panel.border = element_rect(color = "black", fill = NA, size = 0.2),
#           axis.line = element_line(color = "black", size = 0.2),
#           axis.ticks = element_line(size = 0.2),
#           #axis.ticks.length = unit(0.1, "cm"),
#           #plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
#           plot.title = element_text(size = 8, hjust = 0.5), 
#           legend.position = "none")+
#     
#     #stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
#     annotate("text", x = 1.5, y = 0.7, label = paste("p =", format(ttest$p.value, digits = 3)))+
#     annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
#   
#   # Add the plot to the list
#   plot_list[[which(cluster_vec == cluster)]] <- g1
# }
# 
# # Arrange the plots in a grid
# grid.arrange(grobs = plot_list, nrow = 1)
# 
# ggsave("output_graph_paired.png", plot = grid.arrange(grobs = plot_list, nrow = 1), width = 13, height = 5, units = "in")
# 
# 
# 
# 
# ##################################### Stats: Stage cancer#########################################
# ##################################### Dataframe for Outcome
# 
# 
# cancer_subset <- subset(x = immune_annotation, Type == 'Tumor' & Patient_ID_SPC != FALSE)
# 
# md.immune_annotation<- cancer_subset@meta.data %>% as.data.table
# 
# length(which(md.immune_annotation$Patient_ID==68 & md.immune_annotation$Type=='Tumor'& md.immune_annotation$seurat_clusters==1))
# 
# # Vector of sample numbers
# sample_numbers <- c(4, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 52, 53, 54, 56, 59, 60, 65, 66, 68, 71, 73, 75)
# 
# # Vector of cluster names
# cluster_vec <- unique(md.immune_annotation$clusters)
# type_vec <- c("Stage II", "Stage III")
# 
# unique(cancer_subset@meta.data$Patient_ID_SPC)
# 
# 
# 
# # Create an empty data frame to store the results
# output_df <- data.frame(sample_numbers = numeric(),
#                         cluster_vec = character(),
#                         type_vec = character(),
#                         count = numeric(),
#                         total_cell = numeric(),
#                         subset_fraction = numeric(),
#                         stringsAsFactors = FALSE)
# 
# # Loop over the variables and add rows to the data frame
# for (type in type_vec) {
#   for (cluster in cluster_vec) {
#     for (number in sample_numbers) {
#       # Calculate the length of the subset
#       subset_length <- length(which(md.immune_annotation$Patient_ID_SPC == number & md.immune_annotation$Stage == type & md.immune_annotation$clusters == cluster))
#       # Calculate the subset fraction
#       total_count <- length(which(md.immune_annotation$Patient_ID_SPC == number & md.immune_annotation$Stage == type))
#       subset_fraction <- subset_length / total_count
#       # Add a row to the data frame
#       new_row <- data.frame(sample_numbers = number,
#                             cluster_vec = cluster,
#                             type_vec = type,
#                             count = subset_length,
#                             total_cell = total_count,
#                             subset_fraction = subset_fraction,
#                             stringsAsFactors = FALSE)
#       output_df <- rbind(output_df, new_row)
#     }
#   }
# }
# 
# 
# output_df
# 
# ################################################################################## Graph
# 
# # Create an empty list to store the plots
# plot_list <- list()
# 
# for (cluster in cluster_vec) {
#   # Subset the data for each type_vec
#   subset_df <- output_df %>% subset(cluster_vec == cluster)
#   stageII_df <- subset_df %>% subset(type_vec == "Stage II")
#   stageIII_df <- subset_df %>% subset(type_vec == "Stage III")
#   # Perform the unpaired t-test
#   ttest <- t.test(stageII_df$subset_fraction, stageIII_df$subset_fraction)
#   wtest <- wilcox.test(stageII_df$subset_fraction, stageIII_df$subset_fraction)
#   g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
#     geom_boxplot(color = "black", outlier.shape = NA)+
#     scale_fill_manual(values = c("lightgreen", "darkgreen")) +
#     labs(title = NULL,
#          x = NULL, y = NULL) +
#     #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
#     theme(panel.background = element_rect(fill = "#FFFFFF"),
#           panel.border = element_rect(color = "black", fill = NA, size = 0.2),
#           axis.line = element_line(color = "black", size = 0.2),
#           axis.ticks = element_line(size = 0.2),
#           #axis.ticks.length = unit(0.1, "cm"),
#           #plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
#           plot.title = element_text(size = 8, hjust = 0.5), 
#           legend.position = "none")
#   
#   #stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")
#   #annotate("text", x = 1.5, y = 0.5, label = paste("p =", format(ttest$p.value, digits = 3)))
#   #annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
#   
#   # Add the plot to the list
#   plot_list[[which(cluster_vec == cluster)]] <- g1
# }
# 
# # Arrange the plots in a grid
# grid.arrange(grobs = plot_list, nrow = 1)
# 
# ############## Stat add
# 
# for (cluster in cluster_vec) {
#   # Subset the data for each type_vec
#   subset_df <- output_df %>% subset(cluster_vec == cluster)
#   stageII_df <- subset_df %>% subset(type_vec == "Stage II")
#   stageIII_df <- subset_df %>% subset(type_vec == "Stage III")
#   # Perform the unpaired t-test
#   ttest <- t.test(stageII_df$subset_fraction, stageIII_df$subset_fraction)
#   wtest <- wilcox.test(stageII_df$subset_fraction, stageIII_df$subset_fraction)
#   
#   g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
#     geom_boxplot(color = "black", outlier.shape = NA)+
#     scale_fill_manual(values = c("brown", "grey")) +
#     labs(title = cluster,
#          x = NULL, y = NULL) +
#     #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
#     theme(panel.background = element_rect(fill = "#FFFFFF"),
#           panel.border = element_rect(color = "black", fill = NA, size = 0.2),
#           axis.line = element_line(color = "black", size = 0.2),
#           axis.ticks = element_line(size = 0.2),
#           #axis.ticks.length = unit(0.1, "cm"),
#           #plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
#           plot.title = element_text(size = 8, hjust = 0.5), 
#           legend.position = "none")+
#     
#     #stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")+
#     annotate("text", x = 1.5, y = 0.7, label = paste("p =", format(ttest$p.value, digits = 3)))+
#     annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
#   
#   # Add the plot to the list
#   plot_list[[which(cluster_vec == cluster)]] <- g1
# }
# 
# # Arrange the plots in a grid
# grid.arrange(grobs = plot_list, nrow = 1)
# 
# ggsave("output_graph_unpaired.png", plot = grid.arrange(grobs = plot_list, nrow = 1), width = 13, height = 5, units = "in")
# 
# 
# 
# ##################################### Stats: Stage cancer IV include #########################################
# ##################################### Dataframe for Outcome
# 
# 
# cancer_subset <- subset(x = immune_annotation, Type == 'Tumor' & Patient_ID_SPC != FALSE)
# 
# md.immune_annotation<- cancer_subset@meta.data %>% as.data.table
# 
# length(which(md.immune_annotation$Patient_ID==68 & md.immune_annotation$Type=='Tumor'& md.immune_annotation$seurat_clusters==1))
# 
# # Vector of sample numbers
# sample_numbers <- c(4, 5, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 51, 52, 53, 54, 56, 59, 60, 62, 64, 65, 66, 68, 71, 73, 75)
# 
# # Vector of cluster names
# cluster_vec <- unique(md.immune_annotation$clusters)
# type_vec <- c("Stage II", "Stage III", "Stage IV")
# 
# unique(cancer_subset@meta.data$Patient_ID_SPC)
# 
# 
# 
# # Create an empty data frame to store the results
# output_df <- data.frame(sample_numbers = numeric(),
#                         cluster_vec = character(),
#                         type_vec = character(),
#                         count = numeric(),
#                         total_cell = numeric(),
#                         subset_fraction = numeric(),
#                         stringsAsFactors = FALSE)
# 
# # Loop over the variables and add rows to the data frame
# for (type in type_vec) {
#   for (cluster in cluster_vec) {
#     for (number in sample_numbers) {
#       # Calculate the length of the subset
#       subset_length <- length(which(md.immune_annotation$Patient_ID_SPC == number & md.immune_annotation$Stage == type & md.immune_annotation$clusters == cluster))
#       # Calculate the subset fraction
#       total_count <- length(which(md.immune_annotation$Patient_ID_SPC == number & md.immune_annotation$Stage == type))
#       subset_fraction <- subset_length / total_count
#       # Add a row to the data frame
#       new_row <- data.frame(sample_numbers = number,
#                             cluster_vec = cluster,
#                             type_vec = type,
#                             count = subset_length,
#                             total_cell = total_count,
#                             subset_fraction = subset_fraction,
#                             stringsAsFactors = FALSE)
#       output_df <- rbind(output_df, new_row)
#     }
#   }
# }
# 
# 
# output_df
# 
# ################################################################################## Graph
# 
# # Create an empty list to store the plots
# plot_list <- list()
# 
# for (cluster in cluster_vec) {
#   # Subset the data for each type_vec
#   subset_df <- output_df %>% subset(cluster_vec == cluster)
#   stageII_df <- subset_df %>% subset(type_vec == "Stage II")
#   stageIII_df <- subset_df %>% subset(type_vec == "Stage III")
#   stageIV_df <- subset_df %>% subset(type_vec == "Stage IV")
#   # Perform the unpaired t-test
#   #ttest <- t.test(stageII_df$subset_fraction, stageIII_df$subset_fraction)
#   #wtest <- wilcox.test(stageII_df$subset_fraction, stageIII_df$subset_fraction)
#   g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction, fill = type_vec)) +
#     geom_boxplot(color = "black", outlier.shape = NA)+
#     scale_fill_manual(values = c("lightgreen", "darkgreen", "brown")) +
#     labs(title = NULL,
#          x = NULL, y = NULL) +
#     #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0, 0, 0)) +
#     theme(panel.background = element_rect(fill = "#FFFFFF"),
#           panel.border = element_rect(color = "black", fill = NA, size = 0.2),
#           axis.line = element_line(color = "black", size = 0.2),
#           axis.ticks = element_line(size = 0.2),
#           axis.text.x = element_text(angle = 45, hjust = 1),
#           #axis.ticks.length = unit(0.1, "cm"),
#           #plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
#           plot.title = element_text(size = 8, hjust = 0.5), 
#           legend.position = "none")
#   
#   #stat_compare_means(aes(x = type_vec, y = subset_fraction, group = sample_numbers), label = "p.format")
#   #annotate("text", x = 1.5, y = 0.5, label = paste("p =", format(ttest$p.value, digits = 3)))
#   #annotate("text", x = 1.5, y = 0.8, label = paste("w.p =", format(wtest$p.value, digits = 3)))
#   
#   # Add the plot to the list
#   plot_list[[which(cluster_vec == cluster)]] <- g1
# }
# 
# # Arrange the plots in a grid
# grid.arrange(grobs = plot_list, nrow = 1)
