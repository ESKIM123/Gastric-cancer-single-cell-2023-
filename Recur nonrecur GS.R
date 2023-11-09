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


##################################### Recur vs Non-recur #########################################

cancer_subset <- subset(x = immune_annotation, Type == 'Tumor' & Outcome != 'Meta' & Molecular_type == "GS" & Patient_ID_SPC != "NA")
cancer_subset_outcome <- cancer_subset

SCpubr::do_BarPlot(sample = cancer_subset_outcome, 
                   group.by = "clusters", 
                   split.by = "Outcome",
                   position = "stack",legend.position = "none",
                   plot.title = "", xlab = "",  ylab = "Number of cells per cluster", font.size = 11)

SCpubr::do_BarPlot(sample = cancer_subset_outcome, 
                   group.by = "clusters", 
                   split.by = "Outcome",
                   position = "Fill",legend.position = "none",
                   plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)


Nonrecur_subset <- subset(x = cancer_subset_outcome, Outcome == 'Non_recur' & Patient_ID_SPC != FALSE)
SCpubr::do_BarPlot(sample = Nonrecur_subset, 
                   group.by = "clusters", 
                   split.by = "Patient_ID_SPC",
                   position = "Fill",legend.position = "none",
                   plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)

Recur_subset <- subset(x = cancer_subset_outcome, Outcome == 'Recur' & Patient_ID_SPC != FALSE)
SCpubr::do_BarPlot(sample = Recur_subset, 
                   group.by = "clusters", 
                   split.by = "Patient_ID_SPC",
                   position = "Fill",legend.position = "none",
                   plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)


##################################### Stats #########################################
##################################### Dataframe for Outcome
md.cancer_subset_outcome<- cancer_subset_outcome@meta.data %>% as.data.table

length(which(md.cancer_subset_outcome$Patient_ID==68 & md.cancer_subset_outcome$Outcome=='Recur'& md.cancer_subset_outcome$seurat_clusters==1))

# Vector of sample numbers
sample_numbers <- c(59, 66, 44, 41, 18, 54, 4, 56, 52, 68, 29, 7, 36, 28, 8)

# Vector of cluster names
cluster_vec <- unique(md.cancer_subset_outcome$clusters)
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
      subset_length <- length(which(md.cancer_subset_outcome$Patient_ID_SPC == number & md.cancer_subset_outcome$Outcome == type & md.cancer_subset_outcome$clusters == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.cancer_subset_outcome$Patient_ID_SPC == number & md.cancer_subset_outcome$Outcome == type))
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
  grid.text("Total cells - GS", y = 0.05)

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




##################################### Save data #########################################

saveRDS(cancer_subset_outcome, file="D:/DATA/Gastric cancer/Cancer_recur_notmeta.rds")

