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
library(immunarch)
library(tibble)
library(vegan)
library(reshape2)


# Get data ----
CD8_subset <- readRDS(file="D:/DATA/Gastric cancer/GS_subsets/CD8_subset.rds")
DimPlot(CD8_subset, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Initial")

########################### TCR ####################################
CD8_subset@meta.data$rownames <- row.names(CD8_subset@meta.data)

colnames(CD8_subset@meta.data)[colnames(CD8_subset@meta.data) == "alpha"] <- "TCR_alpha"
colnames(CD8_subset@meta.data)[colnames(CD8_subset@meta.data) == "beta"] <- "TCR_beta"

# create a new column that concatenates values from columns A and B
CD8_subset@meta.data$concatenated <- ifelse(!is.na(CD8_subset@meta.data$TCR_alpha),
                                                  paste(CD8_subset@meta.data$TCR_alpha, CD8_subset@meta.data$TCR_beta, sep = ""), "")


view(CD8_subset@meta.data)


######################################## By patient ####
# 환자별로 나눠서 Rank 값들을 넣어줍니다. 


CD8_rank_total <- CD8_subset
md.CD8_rank_total<- CD8_rank_total@meta.data 

patient_id = c(4, 5, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 51, 52, 53, 54, 56, 59, 60, 62, 64, 65, 66, 68, 71, 73, 75)


for (i in patient_id){
  CD8_rank <- subset(x = md.CD8_rank_total, Type == 'Tumor' & Patient_ID_SPC == i)
  
  
  CD8_rank$concatenated <- ifelse(!is.na(CD8_rank$TCR_alpha),
                                  paste(CD8_rank$TCR_alpha, CD8_subset$TCR_beta, sep = ""), "")
  
  # create a new column with the rank numbers
  freq <- table(CD8_rank$concatenated)
  rank <- rank(-freq, ties.method = "min")
  CD8_rank$p_rank <- rank[CD8_rank$concatenated]
  
  CD8_rank$p_top50 <- ifelse(CD8_rank$p_rank <= 50, TRUE, FALSE)
  CD8_rank$p_top20 <- ifelse(CD8_rank$p_rank <= 20, TRUE, FALSE)
  CD8_rank$p_top10 <- ifelse(CD8_rank$p_rank <= 10, TRUE, FALSE)
  CD8_rank$p_top5 <- ifelse(CD8_rank$p_rank <= 5, TRUE, FALSE)
  
  CD8_rank_name <- paste("CD8_rank_", i, sep="")
  assign(CD8_rank_name, CD8_rank)
}

CD8_rank_total@meta.data$p_rank <- NA
CD8_rank_total@meta.data$p_top50 <- NA 
CD8_rank_total@meta.data$p_top20 <- NA
CD8_rank_total@meta.data$p_top10 <- NA
CD8_rank_total@meta.data$p_top5 <- NA

for (j in patient_id) {
  df_name <- paste0("CD8_rank_",j)
  df <- get(df_name)
  CD8_rank_total@meta.data$concatenated[CD8_rank_total@meta.data$rownames %in% df$rownames] <- df$concatenated[df$rownames %in% CD8_rank_total@meta.data$rownames]
  CD8_rank_total@meta.data$p_rank[CD8_rank_total@meta.data$rownames %in% df$rownames] <- df$p_rank[df$rownames %in% CD8_rank_total@meta.data$rownames]
  CD8_rank_total@meta.data$p_top50[CD8_rank_total@meta.data$rownames %in% df$rownames] <- df$p_top50[df$rownames %in% CD8_rank_total@meta.data$rownames]
  CD8_rank_total@meta.data$p_top20[CD8_rank_total@meta.data$rownames %in% df$rownames] <- df$p_top20[df$rownames %in% CD8_rank_total@meta.data$rownames]
  CD8_rank_total@meta.data$p_top10[CD8_rank_total@meta.data$rownames %in% df$rownames] <- df$p_top10[df$rownames %in% CD8_rank_total@meta.data$rownames]
  CD8_rank_total@meta.data$p_top5[CD8_rank_total@meta.data$rownames %in% df$rownames] <- df$p_top5[df$rownames %in% CD8_rank_total@meta.data$rownames]
}


# 각 환자별 top5, top 10,20,50을 표시합니다. 
CD8_rank_total_total <- subset(x = CD8_rank_total, Type == 'Tumor')
DimPlot(CD8_rank_total_total, reduction = 'umap', split.by = "p_top10", raster=FALSE, label = TRUE) + ggtitle("Total")
DimPlot(CD8_rank_total_total, reduction = 'umap', split.by = "p_top5", raster=FALSE, label = TRUE) + ggtitle("Total")

CD8_rank_total_MSI <- subset(x = CD8_rank_total, Type == 'Tumor' & MSI_type == 'MSI')
DimPlot(CD8_rank_total_MSI, reduction = 'umap', split.by = "p_top10", raster=FALSE, label = TRUE) + ggtitle("MSI")
DimPlot(CD8_rank_total_MSI, reduction = 'umap', split.by = "p_top5", raster=FALSE, label = TRUE) + ggtitle("MSI")

CD8_rank_total_MSS <- subset(x = CD8_rank_total, Type == 'Tumor' & MSI_type == 'MSS')
DimPlot(CD8_rank_total_MSS, reduction = 'umap', split.by = "p_top10", raster=FALSE, label = TRUE) + ggtitle("MSS")
DimPlot(CD8_rank_total_MSS, reduction = 'umap', split.by = "p_top5", raster=FALSE, label = TRUE) + ggtitle("MSS")

######################################## View expansion ####

CD8_expansion <- subset(x = CD8_rank_total, Type == 'Tumor' & Outcome != 'Meta' & Patient_ID_SPC != F)

CD8_expansion@meta.data$expansion <- ""
expansion_counts <- table(CD8_expansion@meta.data$concatenated)

for (i in 1:length(CD8_expansion@meta.data$concatenated)) {
  if (!is.na(expansion_counts[CD8_expansion@meta.data$concatenated[i]])) {
    if (expansion_counts[CD8_expansion@meta.data$concatenated[i]] >= 10) {
      CD8_expansion@meta.data$expansion[i] <- "10"
    } else if (expansion_counts[CD8_expansion@meta.data$concatenated[i]] >= 5) {
      CD8_expansion@meta.data$expansion[i] <- "5"
    } else if (expansion_counts[CD8_expansion@meta.data$concatenated[i]] >= 2) {
      CD8_expansion@meta.data$expansion[i] <- "2"
    } else {
      CD8_expansion@meta.data$expansion[i] <- "0"
    }
  } else {
    CD8_expansion@meta.data$expansion[i] <- "NA"
  }
}


CD8_expansion_mod <- CD8_expansion
CD8_expansion_mod@meta.data$expansion[CD8_expansion_mod@meta.data$expansion == "NA"] <- "0"
mapping <- c("0" = 0, "2" = 0.25, "5" = 0.5, "10" = 1)
CD8_expansion_mod@meta.data$expansion_continuous <- mapping[as.character(CD8_expansion_mod@meta.data$expansion)]


# Show as UMAP
DimPlot(CD8_expansion_mod, reduction = 'umap', group.by = "expansion", raster=FALSE, label = F, pt.size = 1, order=T,
        cols = c('0' = '#ebe6dd', '2' = '#b59d69', '5' = '#7d6228', '10' = '#402e08')) + NoLegend() +ggtitle("Total expanded clones") 

FeaturePlot(CD8_expansion_mod, features = "expansion_continuous", pt.size = 1, order=T)

SCpubr::do_NebulosaPlot(sample = CD8_expansion_mod, features = "expansion_continuous",
                        pt.size = 1,plot.title = "",border.size = 0, legend.position = "none")

view(CD8_expansion_mod@meta.data)

# Show as UMAP

CD8_expansion_mod_MSI <- subset(x = CD8_expansion_mod, MSI_type == 'MSI')
DimPlot(CD8_expansion_mod_MSI, reduction = 'umap', group.by = "expansion", raster=FALSE, label = F, pt.size = 1, order=T,
        cols = c('0' = '#ebe6dd', '2' = '#b59d69', '5' = '#7d6228', '10' = '#402e08')) + NoLegend() +ggtitle("") 

CD8_expansion_mod_MSS <- subset(x = CD8_expansion_mod, MSI_type == 'MSS')
DimPlot(CD8_expansion_mod_MSS, reduction = 'umap', group.by = "expansion", raster=FALSE, label = F, pt.size = 1, order=T,
        cols = c('0' = '#ebe6dd', '2' = '#b59d69', '5' = '#7d6228', '10' = '#402e08')) + NoLegend() +ggtitle("") 


############################################################################## STATS ########################

# Count TCR counts
sum(CD8_expansion@meta.data$expansion == "10", na.rm = TRUE)
sum(CD8_expansion@meta.data$expansion == "5", na.rm = TRUE)
sum(CD8_expansion@meta.data$expansion == "2", na.rm = TRUE)
sum(CD8_expansion@meta.data$expansion == "0", na.rm = TRUE)
sum(CD8_expansion@meta.data$expansion == "NA", na.rm = TRUE)

sum(CD8_expansion@meta.data$expansion == "10", na.rm = TRUE) +
  sum(CD8_expansion@meta.data$expansion == "5", na.rm = TRUE) +
  sum(CD8_expansion@meta.data$expansion == "2", na.rm = TRUE) +
  sum(CD8_expansion@meta.data$expansion == "0", na.rm = TRUE) 

CD8_expansion
sum(CD8_expansion@meta.data$expansion_2 == "2", na.rm = TRUE)




# Create example data
df_total <- data.frame(x = c("Total", "TCR"),
                 y = c(17046, 12029))
x_order = c("Total", "TCR")
df_total$x <- factor(df_total$x, levels = x_order)

colors <- c('Total' = 'grey', 'TCR' = 'red')

# Create the bar graph
ggplot(df_total, aes(x = x, y = y, fill = x)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  labs(x = "", y = "", title = "")+
  theme(panel.background = element_rect(fill = "#FFFFFF"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.2),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))


# Create example data
df <- data.frame(x = c("N =1", "N ≥2", "N ≥5", "N ≥10"),
                 y = c(9935, 1220, 389, 485))
x_order = c("N =1", "N ≥2", "N ≥5", "N ≥10")
df$x <- factor(df$x, levels = x_order)

colors <- c('N =1' = '#ebe6dd', 'N ≥2' = '#b59d69', 'N ≥5' = '#7d6228', 'N ≥10' = '#402e08')

# Create the bar graph
ggplot(df, aes(x = x, y = y, fill = x)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  labs(x = "", y = "", title = "")+
  theme(panel.background = element_rect(fill = "#FFFFFF"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.2),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))



# Fraction of expanded clones for each clusters ####
CD8_expansion@meta.data$expansion_2 <- CD8_expansion@meta.data$expansion
CD8_expansion@meta.data$expansion_2[CD8_expansion@meta.data$expansion_2 == "5"] <- "2"
CD8_expansion@meta.data$expansion_2[CD8_expansion@meta.data$expansion_2 == "10"] <- "2"


md.CD8_expansion<- CD8_expansion@meta.data %>% as.data.table

unique_values_2 <- unique(md.CD8_expansion$concatenated[md.CD8_expansion$expansion_2 == 2])
length(unique_values_2)

unique_values_0 <- unique(md.CD8_expansion$concatenated[md.CD8_expansion$expansion_2 == 0])
length(unique_values_0)


# Create example data
df_unique <- data.frame(x = c("Non-expanded", "Expanded"),
                       y = c(9935, 598))
x_order = c("Non-expanded", "Expanded")
df_unique$x <- factor(df_unique$x, levels = x_order)

colors <- c('Non-expanded' = '#ebe6dd', 'Expanded' = 'brown')

# Create the bar graph
ggplot(df_unique, aes(x = x, y = y, fill = x)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  labs(x = "", y = "", title = "")+
  theme(panel.background = element_rect(fill = "#FFFFFF"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.2),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))








# Fraction of expanded clones for each clusters ####
CD8_expansion@meta.data$expansion_2 <- CD8_expansion@meta.data$expansion
CD8_expansion@meta.data$expansion_2[CD8_expansion@meta.data$expansion_2 == "5"] <- "2"
CD8_expansion@meta.data$expansion_2[CD8_expansion@meta.data$expansion_2 == "10"] <- "2"

CD8_expansion@meta.data$expansion_5 <- CD8_expansion@meta.data$expansion
CD8_expansion@meta.data$expansion_5[CD8_expansion@meta.data$expansion_5 == "2"] <- "0"
CD8_expansion@meta.data$expansion_5[CD8_expansion@meta.data$expansion_5 == "5"] <- "5"
CD8_expansion@meta.data$expansion_5[CD8_expansion@meta.data$expansion_5 == "10"] <- "5"

CD8_expansion@meta.data$expansion_10 <- CD8_expansion@meta.data$expansion
CD8_expansion@meta.data$expansion_10[CD8_expansion@meta.data$expansion_10 == "2"] <- "0"
CD8_expansion@meta.data$expansion_10[CD8_expansion@meta.data$expansion_10 == "5"] <- "0"
CD8_expansion@meta.data$expansion_10[CD8_expansion@meta.data$expansion_10 == "10"] <- "10"



patient_id_nonmeta = c(4, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 52, 53, 54, 56, 59, 60, 65, 66, 68, 71, 73, 75)
cluster_numbers = c("C0","C1","C2","C3","C4")
#######################################################

# 클러스터 1에 Expanded clone / 전체 T cell 수 (TCR이 확인된)



TCR2_output_df <- data.frame(sample_ID = numeric(),
                            cluster_vec = character(),
                            Fraction_TCR = numeric(),
                            Total_TCR = numeric(),
                            stringsAsFactors = FALSE)

for (i in patient_id_nonmeta) {
  for (j in cluster_numbers){
  Expanded <- sum(CD8_expansion@meta.data$Patient_ID_SPC == i & CD8_expansion@meta.data$expansion_2 == "2" &  CD8_expansion@meta.data$seurat_cluster == j, na.rm = TRUE)
  
  Total <- sum(CD8_expansion@meta.data$Patient_ID_SPC == i & CD8_expansion@meta.data$expansion_2 == "2", na.rm = TRUE) + sum(CD8_expansion@meta.data$Patient_ID_SPC == i & CD8_expansion@meta.data$expansion_2 == "0", na.rm = TRUE)
  
  
  Fraction <- Expanded / Total * 100
  
  new_row <- data.frame(sample_ID = i,
                        cluster_vec = j,
                        Fraction_TCR = Fraction,
                        Total_TCR = Total,
                        stringsAsFactors = FALSE)
  TCR2_output_df <- rbind(TCR2_output_df, new_row)
  }
}

TCR2_output_df
write.csv(TCR2_output_df, file = "D:/DATA/Gastric cancer/GS_subsets/CD8_subset_TCR2.csv", row.names = FALSE)


#######################################################
TCR5_output_df <- data.frame(sample_ID = character(),
                            cluster_vec = character(),
                            Fraction_TCR = numeric(),
                            Total_TCR = numeric(),
                            stringsAsFactors = FALSE)

for (i in patient_id_nonmeta) {
  for (j in cluster_numbers){
    Expanded <- sum(CD8_expansion@meta.data$Patient_ID_SPC == i & CD8_expansion@meta.data$expansion_5 == "5" &  CD8_expansion@meta.data$seurat_cluster == j, na.rm = TRUE)
    
    Total <- sum(CD8_expansion@meta.data$Patient_ID_SPC == i & CD8_expansion@meta.data$expansion_5 == "5", na.rm = TRUE) + sum(CD8_expansion@meta.data$Patient_ID_SPC == i & CD8_expansion@meta.data$expansion_5 == "0", na.rm = TRUE)
    
    
    Fraction <- Expanded / Total * 100
    
    new_row <- data.frame(sample_ID = i,
                          cluster_vec = j,
                          Fraction_TCR = Fraction,
                          Total_TCR = Total,
                          stringsAsFactors = FALSE)
    TCR5_output_df <- rbind(TCR5_output_df, new_row)
  }
}

TCR5_output_df
write.csv(TCR5_output_df, file = "D:/DATA/Gastric cancer/GS_subsets/CD8_subset_TCR5.csv", row.names = FALSE)

######################################################
TCR10_output_df <- data.frame(sample_ID = numeric(),
                             cluster_vec = character(),
                             Fraction_TCR = numeric(),
                             Total_TCR = numeric(),
                             stringsAsFactors = FALSE)

for (i in patient_id_nonmeta) {
  for (j in cluster_numbers){
    Expanded <- sum(CD8_expansion@meta.data$Patient_ID_SPC == i & CD8_expansion@meta.data$expansion_10 == "10" &  CD8_expansion@meta.data$seurat_cluster == j, na.rm = TRUE)
    
    Total <- sum(CD8_expansion@meta.data$Patient_ID_SPC == i & CD8_expansion@meta.data$expansion_10 == "10", na.rm = TRUE) + sum(CD8_expansion@meta.data$Patient_ID_SPC == i & CD8_expansion@meta.data$expansion_10 == "0", na.rm = TRUE)
    
    
    Fraction <- Expanded / Total * 100
    
    new_row <- data.frame(sample_ID = i,
                          cluster_vec = j,
                          Fraction_TCR = Fraction,
                          Total_TCR = Total,
                          stringsAsFactors = FALSE)
    TCR10_output_df <- rbind(TCR10_output_df, new_row)
  }
}

TCR10_output_df
write.csv(TCR10_output_df, file = "D:/DATA/Gastric cancer/GS_subsets/CD8_subset_TCR10.csv", row.names = FALSE)



######################################################

# 클러스터 1에 Expanded clone / 전체 Expanded clone 수 (TCR이 확인된)

TCR2ex_output_df <- data.frame(sample_ID = numeric(),
                              cluster_vec = character(),
                              Fraction_TCR = numeric(),
                              Total_TCR = numeric(),
                              stringsAsFactors = FALSE)

for (i in patient_id_nonmeta) {
  for (j in cluster_numbers){
    Expanded <- sum(CD8_expansion@meta.data$Patient_ID_SPC == i & CD8_expansion@meta.data$expansion_2 == "2" &  CD8_expansion@meta.data$seurat_cluster == j, na.rm = TRUE)
    
    Total <- sum(CD8_expansion@meta.data$Patient_ID_SPC == i & CD8_expansion@meta.data$expansion_2 == "2", na.rm = TRUE)
    
    Fraction <- Expanded / Total * 100
    
    new_row <- data.frame(sample_ID = i,
                          cluster_vec = j,
                          Fraction_TCR = Fraction,
                          Total_TCR = Total,
                          stringsAsFactors = FALSE)
    TCR2ex_output_df <- rbind(TCR2ex_output_df, new_row)
  }
}

TCR2ex_output_df
write.csv(TCR2ex_output_df, file = "D:/DATA/Gastric cancer/GS_subsets/CD8_subset_TCR2ex.csv", row.names = FALSE)


cluster_vec <- c("C0", "C1", "C2", "C3", "C4")


ggplot(TCR2ex_output_df, aes(x = cluster_vec, y = Fraction_TCR))+
  geom_boxplot(color = "black", outlier.shape = NA, alpha=0.8, fill = c("#e69f00", "#56b4e9", "#009e73", "#f0e442", "#0072b2"))+
  geom_point(shape = 21, color = "black", alpha=0.8, fill = "#4A3F4E")+
  geom_jitter(width = 0.1)+
    ggtitle("") + 
  xlab("") + ylab("") +
  theme(panel.background = element_rect(fill = "#FFFFFF"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.2),
        axis.line = element_line(color = "black", size = 0.2),
        axis.ticks = element_line(size = 0.2),
        axis.ticks.length = unit(0.1, "cm"),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        legend.position = "none")


ggplot(TCR2ex_output_df, aes(x = cluster_vec, y = Fraction_TCR))+
  geom_boxplot(color = "black", outlier.shape = NA, alpha=0.8, fill = c("#e69f00", "#56b4e9", "#009e73", "#f0e442", "#0072b2"))+
  geom_point(shape = 21, color = "black", alpha=0.8, fill = "grey")+
  geom_line(aes(group = sample_ID), alpha=0.5)+
  ggtitle("") + 
  xlab("") + ylab("") +
  theme(panel.background = element_rect(fill = "#FFFFFF"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.2),
        axis.line = element_line(color = "black", size = 0.2),
        axis.ticks = element_line(size = 0.2),
        axis.ticks.length = unit(0.1, "cm"),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        legend.position = "none")

########################################################## Draw heatmap ####


# Get unique values of sample_ID and cluster_vec
row_names <- unique(TCR2ex_output_df$sample_ID)
col_names <- unique(TCR2ex_output_df$cluster_vec)

# Create an empty dataframe
new_df <- data.frame(matrix(NA, nrow = length(row_names), ncol = length(col_names)))
rownames(new_df) <- row_names
colnames(new_df) <- col_names

# Populate the new dataframe with Fraction_TCR data
for (i in row_names) {
  for (j in col_names){
  row_name <- TCR2ex_output_df$sample_ID[i]
  col_name <- TCR2ex_output_df$cluster_vec[j]
  value <- subset(TCR2ex_output_df, sample_ID == i & cluster_vec == j)$Fraction_TCR
  
  new_df[row_names==i, col_names==j] <- value
  }
  }

heatmap(as.matrix(new_df, ), Colv = NA, cexCol = 1)


pheatmap(as.matrix(new_df, ), Colv = NA, cexCol = 1, scale = "row", cluster_cols = F, clustering_method = "ward", cellwidth=25)















