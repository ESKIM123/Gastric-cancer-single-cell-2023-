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
CD4_annotation <- readRDS(file="D:/DATA/Gastric cancer/GS_subsets/CD4_subset.rds")
DimPlot(CD4_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Initial")

##################################### wtp53 vs mtp53 #########################################

CD4_annotation[["clusters"]]  <- CD4_annotation[["seurat_clusters"]] 
CD4_annotation[["TRM_NTRM"]]  <- NA

cancer_subset_TRM_NTRM <- subset(x = CD4_annotation, Type == 'Tumor' & Patient_ID_SPC != FALSE)
cancer_subset_TRM_NTRM <- subset(x = cancer_subset_TRM_NTRM, Stage != 'Stage IV')

patient_id_TRM = c(15,52,35,50,75,7,40)

patient_id_NTRM = c(28,16,65,66,53,18,56,44,4,41,59,8,22,73,27,68,23,36,20,54,29)


for (i in patient_id_TRM){
  cancer_subset_TRM_NTRM@meta.data$TRM_NTRM[cancer_subset_TRM_NTRM@meta.data$Patient_ID_SPC == i] <- "TRM"
}

for (j in patient_id_NTRM){
  cancer_subset_TRM_NTRM@meta.data$TRM_NTRM[cancer_subset_TRM_NTRM@meta.data$Patient_ID_SPC == j] <- "NTRM"
}

cancer_subset_TRM_NTRM <- subset(x = cancer_subset_TRM_NTRM, TRM_NTRM != "NA")


cancer_TRM_subset <- subset(x = cancer_subset_TRM_NTRM, TRM_NTRM == 'TRM')
cancer_NTRM_subset <- subset(x = cancer_subset_TRM_NTRM, TRM_NTRM == 'NTRM')


SCpubr::do_BarPlot(sample = cancer_subset_TRM_NTRM, 
                   group.by = "clusters", 
                   split.by = "TRM_NTRM",
                   position = "Fill",legend.position = "bottom",
                   plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)


cancer_NTRM_subset <- subset(x = cancer_TRM_subset, TRM_NTRM == 'NTRM')
SCpubr::do_BarPlot(sample = cancer_NTRM_subset, 
                   group.by = "clusters", 
                   split.by = "Patient_ID_SPC",
                   position = "Fill",legend.position = "bottom",
                   plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)

cancer_TRM_subset <- subset(x = cancer_TRM_subset, TRM_NTRM == 'TRM')
SCpubr::do_BarPlot(sample = cancer_TRM_subset, 
                   group.by = "clusters", 
                   split.by = "Patient_ID_SPC",
                   position = "Fill",legend.position = "bottom",
                   plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)


##################################### Stats #########################################
##################################### Dataframe for Outcome
md.CD4_annotation<- cancer_subset_TRM_NTRM@meta.data %>% as.data.table

length(which(md.CD4_annotation$Patient_ID==68 & md.CD4_annotation$Type=='Tumor'& md.CD4_annotation$seurat_clusters==1))

# Vector of sample numbers
sample_numbers <- c(15,52,35,50,75,7,40,28,16,65,66,53,18,56,44,4,41,59,8,22,73,27,68,23,36,20,54,29)

# Vector of cluster names
cluster_vec <- unique(md.CD4_annotation$clusters)
type_vec <- c("NTRM", "TRM")

unique(cancer_TRM_subset@meta.data$Patient_ID_SPC)



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
      subset_length <- length(which(md.CD4_annotation$Patient_ID_SPC == number & md.CD4_annotation$TRM_NTRM == type & md.CD4_annotation$clusters == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.CD4_annotation$Patient_ID_SPC == number & md.CD4_annotation$TRM_NTRM == type))
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
  tumor_df <- subset_df %>% subset(type_vec == "TRM") 
  normal_df <- subset_df %>% subset(type_vec == "NTRM")
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
  tumor_df <- subset_df %>% subset(type_vec == "TRM") 
  normal_df <- subset_df %>% subset(type_vec == "NTRM")
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