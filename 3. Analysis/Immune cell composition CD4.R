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
library(fgsea)
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyverse)
library(data.table)
library(presto)
library(ggrepel)
library(gplots)
library(RColorBrewer)


# Get data ----
CD8_subset_anno <- readRDS(file="D:/DATA/Gastric cancer/GS_subsets/CD8_subset.rds")
DimPlot(CD8_subset_anno, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Initial")

view(CD8_subset_anno@meta.data)

# Subset 
immune_annotation <- subset(x = CD8_subset_anno, Type == 'Tumor' & Patient_ID_SPC != FALSE &
                              Patient_ID_SPC != 8 & Patient_ID_SPC != 18 & Patient_ID_SPC != 52)

CD8_TCGAtype <- immune_annotation

view(CD8_TCGAtype@meta.data)

# Show as UMAP

DimPlot(CD8_TCGAtype, reduction = 'umap', group.by = "TCGA", raster=FALSE, label = F, pt.size = 1, order=T,
        cols = c('MSI' = '#e69f00', 'CIN' = '#ebe6dd', 'GS' = '#ebe6dd')) + NoLegend() +ggtitle("Total expanded clones") 

DimPlot(CD8_TCGAtype, reduction = 'umap', group.by = "TCGA", raster=FALSE, label = F, pt.size = 1, order=T,
                cols = c('MSI' = '#e69f00', 'CIN' = '#ebe6dd', 'GS' = '#ebe6dd')) + NoLegend() +ggtitle("") + NoAxes()

DimPlot(CD8_TCGAtype, reduction = 'umap', group.by = "TCGA", raster=FALSE, label = F, pt.size = 1, order=T,
        cols = c('MSI' = '#ebe6dd', 'CIN' = '#56b4e9', 'GS' = '#ebe6dd')) + NoLegend() +ggtitle("Total expanded clones")

DimPlot(CD8_TCGAtype, reduction = 'umap', group.by = "TCGA", raster=FALSE, label = F, pt.size = 1, order=T,
        cols = c('MSI' = '#ebe6dd', 'CIN' = '#56b4e9', 'GS' = '#ebe6dd')) + NoLegend() +ggtitle("") + NoAxes()

DimPlot(CD8_TCGAtype, reduction = 'umap', group.by = "TCGA", raster=FALSE, label = F, pt.size = 1, order=T,
        cols = c('MSI' = '#ebe6dd', 'CIN' = '#ebe6dd', 'GS' = '#009e73')) + NoLegend() +ggtitle("Total expanded clones") 

DimPlot(CD8_TCGAtype, reduction = 'umap', group.by = "TCGA", raster=FALSE, label = F, pt.size = 1, order=T,
        cols = c('MSI' = '#ebe6dd', 'CIN' = '#ebe6dd', 'GS' = '#009e73')) + NoLegend() +ggtitle("") + NoAxes()



dittoBarPlot(CD8_TCGAtype, "seurat_clusters", group.by = "TCGA", scale = "count", main = "")+ 
  NoLegend()+ xlab(NULL)+ ylab(NULL)+
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))

dittoBarPlot(CD8_TCGAtype, "seurat_clusters", group.by = "TCGA", scale = "percent", main = "")+ 
  NoLegend()+ xlab(NULL)+ ylab(NULL)+
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))

dittoBarPlot(CD8_TCGAtype, "seurat_clusters", group.by = "TCGA", scale = "percent", main = "")+ 
  xlab(NULL)+ ylab(NULL)+
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))


dittoBarPlot(CD8_TCGAtype, "seurat_clusters", group.by = "Patient_ID_SPC", scale = "percent", main = "")+ 
  NoLegend()+ xlab(NULL)+ ylab(NULL)+
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))

##################################### Stats #########################################
##################################### Dataframe for Outcome
md.immune_annotation<- CD8_TCGAtype@meta.data %>% as.data.table

length(which(md.immune_annotation$Patient_ID==68 & md.immune_annotation$Type=='Tumor'& md.immune_annotation$clusters==1))

# Vector of sample numbers
sample_numbers <- unique(CD8_TCGAtype@meta.data$Patient_ID_SPC)

# Vector of cluster names
cluster_vec <- unique(md.immune_annotation$seurat_clusters)
cluster_vec
type_vec <- unique(md.immune_annotation$TCGA)
type_vec


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
      subset_length <- length(which(md.immune_annotation$Patient_ID_SPC == number & md.immune_annotation$TCGA == type & md.immune_annotation$seurat_clusters == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.immune_annotation$Patient_ID_SPC == number & md.immune_annotation$TCGA == type))
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
  MSI_df <- subset_df %>% subset(type_vec == "MSI")
  CIN_df <- subset_df %>% subset(type_vec == "CIN")
  GS_df <- subset_df %>% subset(type_vec == "GS")
  
  # Perform the unpaired t-test
  ttest01 <- t.test(MSI_df$subset_fraction, CIN_df$subset_fraction)
  ttest02 <- t.test(MSI_df$subset_fraction, GS_df$subset_fraction)
  ttest12 <- t.test(CIN_df$subset_fraction, GS_df$subset_fraction)
  
  
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
grid.arrange(grobs = plot_list, nrow = 3)

############## Stat add

for (cluster in cluster_vec) {
  # Subset the data for each type_vec
  subset_df <- output_df %>% subset(cluster_vec == cluster)
  MSI_df <- subset_df %>% subset(type_vec == "MSI")
  CIN_df <- subset_df %>% subset(type_vec == "CIN")
  GS_df <- subset_df %>% subset(type_vec == "GS")
  
  # Perform the unpaired t-test
  ttest01 <- t.test(MSI_df$subset_fraction, CIN_df$subset_fraction)
  ttest02 <- t.test(MSI_df$subset_fraction, GS_df$subset_fraction)
  ttest12 <- t.test(CIN_df$subset_fraction, GS_df$subset_fraction)
  
  
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



# Heatmap ------------------



##################################### Dataframe for Outcome
md.immune_annotation<- immune_annotation@meta.data %>% as.data.table

length(which(md.immune_annotation$Patient_ID==68 & md.immune_annotation$Type=='Tumor'& md.immune_annotation$seurat_clusters==1))

# Vector of sample numbers
sample_numbers <- unique(immune_annotation$Patient_ID_SPC)

# Vector of cluster names
cluster_vec <- unique(md.immune_annotation$seurat_clusters)
type_vec <- c("Tumor", "Normal")


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
      subset_length <- length(which(md.immune_annotation$Patient_ID_SPC == number & md.immune_annotation$Type == type & md.immune_annotation$seurat_clusters == cluster))
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


output_df_tumor <- subset(output_df, type_vec %in% "Tumor")

# Get unique values of sample_numbers and cluster_vec
unique_sample_numbers <- unique(output_df_tumor$sample_numbers)
unique_cluster_vec <- unique(output_df_tumor$cluster_vec)

# Create an empty dataframe with rows and columns based on unique values
new_df <- data.frame(
  sample_numbers = unique_sample_numbers,
  stringsAsFactors = FALSE
)

# Add columns for unique cluster_vec values
for (cluster in unique_cluster_vec) {
  new_df[, cluster] <- 0  # Initialize all values to 0
}

# Fill in the values from output_df_tumor
for (i in 1:nrow(output_df_tumor)) {
  row_index <- match(output_df_tumor$sample_numbers[i], new_df$sample_numbers)
  col_name <- output_df_tumor$cluster_vec[i]
  new_df[row_index, col_name] <- output_df_tumor$subset_fraction[i]
}



##################################### Clustering



colnames(new_df)
 new2_df <- new_df
# new2_df <- new_df[, !(names(new_df) %in% c("3"))]
unique_cluster_vec2 <- c("0","1","2", "3", "4", "5", "6")


# un-scaling
new_df_unscaled <- new2_df


# Set row names to be the sample_numbers
row.names(new_df_unscaled) <- new_df_unscaled$sample_numbers

new_df_unscaled$sample_numbers <- NULL

# Create a custom pastel color palette
pastel_palette <- brewer.pal(n = 9, name = "YlGnBu")

# Create the heatmap with a custom color palette
heatmap.2(as.matrix(new_df_unscaled), dendrogram = "row", scale = "none", key = TRUE, col = pastel_palette, margins=c(11,11))



# Min-max scaling
new_df_scaled <- new2_df

for (col in unique_cluster_vec2) {
  col_data <- new_df_scaled[, col]
  new_df_scaled[, col] <- (col_data - min(col_data)) / (max(col_data) - min(col_data))
}


# Set row names to be the sample_numbers
row.names(new_df_scaled) <- new_df_scaled$sample_numbers
new_df_scaled$sample_numbers <- NULL

# Create a custom pastel color palette
pastel_palette <- brewer.pal(n = 9, name = "YlGnBu")

# Create the heatmap with a custom color palette
heatmap.2(as.matrix(new_df_scaled), dendrogram = "row", scale = "none", key = TRUE, col = pastel_palette, margins=c(15,3),
          density.info="none",  trace = "none",
          cexRow = 1,
          srtCol = 45, 
          sepCol = 0.5)



# z-score scaling
for (col in unique_cluster_vec2) {
  col_data <- new_df_scaled[, col]
  new_df_scaled[, col] <- (col_data - mean(col_data)) / sd(col_data)
}

# Create a custom color palette function
custom_color_palette <- colorRampPalette(c("blue", "white", "red"))

# Create the heatmap with the custom color palette
heatmap.2(as.matrix(new_df_scaled), dendrogram = "row", scale = "none", key = TRUE, col = custom_color_palette(100),
          margins = c(15, 3), density.info = "none", trace = "none",
          cexRow = 1, srtCol = 45, sepCol = 0.5)