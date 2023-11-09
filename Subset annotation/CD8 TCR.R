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
CD8_subset <- readRDS(file="D:/DATA/Gastric cancer/CD8_subset.rds")

CD8_subset@meta.data$rownames <- row.names(CD8_subset@meta.data)

colnames(CD8_subset@meta.data)[colnames(CD8_subset@meta.data) == "alpha"] <- "TCR_alpha"
colnames(CD8_subset@meta.data)[colnames(CD8_subset@meta.data) == "beta"] <- "TCR_beta"

view(CD8_subset@meta.data)

DimPlot(CD8_subset, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Initial")

################################### Divide CD8 by type #########################

CD8_subset_recur <- subset(x = CD8_subset, Type == 'Tumor' & Outcome == 'Recur')
CD8_subset_nonrecur <- subset(x = CD8_subset, Type == 'Tumor' & Outcome == 'Non_recur')

######################################## Recur

# create a new column that concatenates values from columns A and B
CD8_subset_recur@meta.data$concatenated <- ifelse(!is.na(CD8_subset_recur@meta.data$TCR_alpha),
                                            paste(CD8_subset_recur@meta.data$TCR_alpha, CD8_subset@meta.data$TCR_beta, sep = ""), "")

# create a new column with the rank numbers
freq <- table(CD8_subset_recur@meta.data$concatenated)
rank <- rank(-freq, ties.method = "min")
CD8_subset_recur@meta.data$rank <- rank[CD8_subset_recur@meta.data$concatenated]

CD8_subset_recur@meta.data$top50 <- ifelse(CD8_subset_recur@meta.data$rank <= 50, TRUE, FALSE)
CD8_subset_recur@meta.data$top20 <- ifelse(CD8_subset_recur@meta.data$rank <= 20, TRUE, FALSE)
CD8_subset_recur@meta.data$top10 <- ifelse(CD8_subset_recur@meta.data$rank <= 10, TRUE, FALSE)

DimPlot(CD8_subset_recur, reduction = 'umap', split.by = "top20", raster=FALSE, label = TRUE) + ggtitle("Recur")

######################################## Non recur

# create a new column that concatenates values from columns A and B
CD8_subset_nonrecur@meta.data$concatenated <- ifelse(!is.na(CD8_subset_nonrecur@meta.data$TCR_alpha),
                                                  paste(CD8_subset_nonrecur@meta.data$TCR_alpha, CD8_subset@meta.data$TCR_beta, sep = ""), "")

# create a new column with the rank numbers
freq <- table(CD8_subset_nonrecur@meta.data$concatenated)
rank <- rank(-freq, ties.method = "min")
CD8_subset_nonrecur@meta.data$rank <- rank[CD8_subset_nonrecur@meta.data$concatenated]

CD8_subset_nonrecur@meta.data$top50 <- ifelse(CD8_subset_nonrecur@meta.data$rank <= 50, TRUE, FALSE)
CD8_subset_nonrecur@meta.data$top20 <- ifelse(CD8_subset_nonrecur@meta.data$rank <= 20, TRUE, FALSE)
CD8_subset_nonrecur@meta.data$top10 <- ifelse(CD8_subset_nonrecur@meta.data$rank <= 10, TRUE, FALSE)

DimPlot(CD8_subset_nonrecur, reduction = 'umap', split.by = "top20", raster=FALSE, label = TRUE) + ggtitle("Nonrecur")



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
CD8_rank_total_recur <- subset(x = CD8_rank_total, Type == 'Tumor' & Outcome == 'Recur')
DimPlot(CD8_rank_total_recur, reduction = 'umap', split.by = "p_top10", raster=FALSE, label = TRUE) + ggtitle("Recur")
DimPlot(CD8_rank_total_recur, reduction = 'umap', split.by = "p_top5", raster=FALSE, label = TRUE) + ggtitle("Recur")

CD8_rank_total_recur <- subset(x = CD8_rank_total, Type == 'Tumor' & Outcome == 'Recur'& p_top5 == T)
DimPlot(CD8_rank_total_recur, reduction = 'umap', raster=FALSE, label = F) + ggtitle("")+NoLegend()+NoAxes()

CD8_rank_total_nonrecur <- subset(x = CD8_rank_total, Type == 'Tumor' & Outcome == 'Non_recur')
DimPlot(CD8_rank_total_nonrecur, reduction = 'umap', split.by = "p_top10", raster=FALSE, label = TRUE) + ggtitle("Non recur")
DimPlot(CD8_rank_total_nonrecur, reduction = 'umap', split.by = "p_top5", raster=FALSE, label = TRUE) + ggtitle("Non recur")

CD8_rank_total_nonrecur <- subset(x = CD8_rank_total, Type == 'Tumor' & Outcome == 'Non_recur'& p_top5 == T)
DimPlot(CD8_rank_total_nonrecur, reduction = 'umap', raster=FALSE, label = F) + ggtitle("")+NoLegend()+NoAxes()









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
        cols = c('0' = '#ebe6dd', '2' = '#b59d69', '5' = '#7d6228', '10' = '#402e08')) + NoLegend() +ggtitle("")

FeaturePlot(CD8_expansion_mod, features = "expansion_continuous", pt.size = 1, order=T)

SCpubr::do_NebulosaPlot(sample = CD8_expansion_mod, features = "expansion_continuous",
                        pt.size = 1,plot.title = "",border.size = 0, legend.position = "none")


# Show as UMAP
CD8_expansion_mod_recur <- subset(x = CD8_expansion_mod, Outcome == 'Recur')
p1 <- DimPlot(CD8_expansion_mod_recur, reduction = 'umap', group.by = "expansion", raster=FALSE, label = F, pt.size = 1, order=T,
        cols = c('0' = '#ebe6dd', '2' = '#b59d69', '5' = '#7d6228', '10' = '#402e08')) + NoLegend() +ggtitle("") +NoAxes()

CD8_expansion_mod_nonrecur <- subset(x = CD8_expansion_mod, Outcome == 'Non_recur')
p2 <- DimPlot(CD8_expansion_mod_nonrecur, reduction = 'umap', group.by = "expansion", raster=FALSE, label = F, pt.size = 1, order=T,
        cols = c('0' = '#ebe6dd', '2' = '#b59d69', '5' = '#7d6228', '10' = '#402e08')) + NoLegend() +ggtitle("") +NoAxes()


p3 <- SCpubr::do_NebulosaPlot(sample = CD8_expansion_mod_recur, features = "expansion_continuous",
                        pt.size = 1,plot.title = "",border.size = 0, legend.position = "none")

p4 <- SCpubr::do_NebulosaPlot(sample = CD8_expansion_mod, features = "expansion_continuous",
                        pt.size = 1,plot.title = "",border.size = 0, legend.position = "none")


grid.arrange(p1, p2, p3, p4, ncol = 2)



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




# Fraction of expanded clones for each clusters ####
CD8_expansion@meta.data$expansion_2 <- CD8_expansion@meta.data$expansion
CD8_expansion@meta.data$expansion_2[CD8_expansion@meta.data$expansion_2 == "5"] <- "2"
CD8_expansion@meta.data$expansion_2[CD8_expansion@meta.data$expansion_2 == "10"] <- "2"


fraction_list <- list()
for (i in 0:7) {
  Expanded <- sum(CD8_expansion@meta.data$expansion_2 == "2" & CD8_expansion@meta.data$seurat_cluster == i, na.rm = TRUE)
  Total <- sum(CD8_expansion@meta.data$expansion_2 == "2" & CD8_expansion@meta.data$seurat_cluster == i, na.rm = TRUE) + 
    sum(CD8_expansion@meta.data$expansion_2 == "0" & CD8_expansion@meta.data$seurat_cluster == i, na.rm = TRUE) 
  Fraction <- Expanded / Total * 100
  fraction_list[[i+1]] <- Fraction
}

df <- data.frame(Cluster = 0:7, Fraction = unlist(fraction_list))

# Create the barplot
ggplot(df, aes(x = Cluster, y = Fraction)) +
  geom_bar(stat = "identity", fill = "#b59d69", color = "black") +
  xlab("Cluster") + ylab("Fraction of expanded clones (≥2) (%)") + theme_classic() + scale_x_continuous(breaks = 0:7) + 
  scale_y_continuous(expand = c(0, 0))








# Fraction of expaned clones for each clusters ####
CD8_expansion@meta.data$expansion_5 <- CD8_expansion@meta.data$expansion
CD8_expansion@meta.data$expansion_5[CD8_expansion@meta.data$expansion_5 == "2"] <- "0"
CD8_expansion@meta.data$expansion_5[CD8_expansion@meta.data$expansion_5 == "5"] <- "2"
CD8_expansion@meta.data$expansion_5[CD8_expansion@meta.data$expansion_5 == "10"] <- "2"


fraction_list <- list()
for (i in 0:7) {
  Expanded <- sum(CD8_expansion@meta.data$expansion_5 == "2" & CD8_expansion@meta.data$seurat_cluster == i, na.rm = TRUE)
  Total <- sum(CD8_expansion@meta.data$expansion_5 == "2" & CD8_expansion@meta.data$seurat_cluster == i, na.rm = TRUE) + 
    sum(CD8_expansion@meta.data$expansion_5 == "0" & CD8_expansion@meta.data$seurat_cluster == i, na.rm = TRUE) 
  Fraction <- Expanded / Total * 100
  fraction_list[[i+1]] <- Fraction
}

df <- data.frame(Cluster = 0:7, Fraction = unlist(fraction_list))

# Create the barplot
ggplot(df, aes(x = Cluster, y = Fraction)) +
  geom_bar(stat = "identity", fill = "#7d6228", color = "black") +
  xlab("Cluster") + ylab("Fraction of Expanded clones (≥5) (%)") + theme_classic() + scale_x_continuous(breaks = 0:7) + scale_y_continuous(expand = c(0, 0))

# Fraction of expaned clones for each clusters ####
CD8_expansion@meta.data$expansion_10 <- CD8_expansion@meta.data$expansion
CD8_expansion@meta.data$expansion_10[CD8_expansion@meta.data$expansion_10 == "2"] <- "0"
CD8_expansion@meta.data$expansion_10[CD8_expansion@meta.data$expansion_10 == "5"] <- "0"
CD8_expansion@meta.data$expansion_10[CD8_expansion@meta.data$expansion_10 == "10"] <- "2"


fraction_list <- list()
for (i in 0:7) {
  Expanded <- sum(CD8_expansion@meta.data$expansion_10 == "2" & CD8_expansion@meta.data$seurat_cluster == i, na.rm = TRUE)
  Total <- sum(CD8_expansion@meta.data$expansion_10 == "2" & CD8_expansion@meta.data$seurat_cluster == i, na.rm = TRUE) + 
    sum(CD8_expansion@meta.data$expansion_10 == "0" & CD8_expansion@meta.data$seurat_cluster == i, na.rm = TRUE) 
  Fraction <- Expanded / Total * 100
  fraction_list[[i+1]] <- Fraction
}

df <- data.frame(Cluster = 0:7, Fraction = unlist(fraction_list))

# Create the barplot
ggplot(df, aes(x = Cluster, y = Fraction)) +
  geom_bar(stat = "identity", fill = "#402e08", color = "black") +
  xlab("Cluster") + ylab("Fraction of Expanded clones (≥10) (%)") + theme_classic() + scale_x_continuous(breaks = 0:7) + scale_y_continuous(expand = c(0, 0))

################################################################################ 각 클러스터별 평균 expanded cells ####
#CD8_expansion <- subset(x = CD8_rank_total, Type == 'Tumor' & Outcome != 'Meta' & Patient_ID_SPC != F)

patient_id_nometa = c(4, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 52, 53, 54, 56, 59, 60, 65, 66, 68, 71, 73, 75)

# Create an empty dataframe to store the results
results_df <- data.frame(Patient_ID_SPC = numeric(),
                         Cluster = numeric(),
                         Fraction = numeric())

for (i in patient_id_nometa){
  CD8_expansion_id <- subset(x = CD8_expansion, Patient_ID_SPC == i)
  CD8_expansion_id@meta.data$expansion_2 <- CD8_expansion_id@meta.data$expansion
  CD8_expansion_id@meta.data$expansion_2[CD8_expansion_id@meta.data$expansion_2 == "2"] <- "2"
  CD8_expansion_id@meta.data$expansion_2[CD8_expansion_id@meta.data$expansion_2 == "5"] <- "2"
  CD8_expansion_id@meta.data$expansion_2[CD8_expansion_id@meta.data$expansion_2 == "10"] <- "2"
  Total_list <- list()
  for (j in 0:7) {
    Expanded <- sum(CD8_expansion_id@meta.data$expansion_2 == "2" & CD8_expansion_id@meta.data$seurat_cluster == j, na.rm = TRUE)
    Total <- sum(CD8_expansion_id@meta.data$expansion_2 == "2" & CD8_expansion_id@meta.data$seurat_cluster == j, na.rm = TRUE) + 
      sum(CD8_expansion_id@meta.data$expansion_2 == "0" & CD8_expansion_id@meta.data$seurat_cluster == j, na.rm = TRUE) 
    Total_list[[j+1]] <- Total
  }
  
  # Create a dataframe with the results for the current patient
  patient_df <- data.frame(Patient_ID_SPC = i, Cluster = 0:7, Total = unlist(Total_list))
  
  # Append the patient dataframe to the overall results dataframe
  results_df <- rbind(results_df, patient_df)
}

# Use dcast to reshape the dataframe
results_reshaped <- dcast(results_df, Patient_ID_SPC ~ Cluster, value.var = "Total")

# Write the results_reshaped dataframe to a CSV file
write.csv(results_reshaped, file = "D:/DATA/Gastric cancer/CD8 TCR results total counts per cluster.csv", row.names = FALSE)


################################################################################ 각 클러스터별 평균 expanded cells ####
#CD8_expansion <- subset(x = CD8_rank_total, Type == 'Tumor' & Outcome != 'Meta' & Patient_ID_SPC != F)

patient_id_nometa = c(4, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 52, 53, 54, 56, 59, 60, 65, 66, 68, 71, 73, 75)

# Create an empty dataframe to store the results
results_df <- data.frame(Patient_ID_SPC = numeric(),
                         Cluster = numeric(),
                         Fraction = numeric())

for (i in patient_id_nometa){
  CD8_expansion_id <- subset(x = CD8_expansion, Patient_ID_SPC == i)
  CD8_expansion_id@meta.data$expansion_2 <- CD8_expansion_id@meta.data$expansion
  CD8_expansion_id@meta.data$expansion_2[CD8_expansion_id@meta.data$expansion_2 == "2"] <- "2"
  CD8_expansion_id@meta.data$expansion_2[CD8_expansion_id@meta.data$expansion_2 == "5"] <- "2"
  CD8_expansion_id@meta.data$expansion_2[CD8_expansion_id@meta.data$expansion_2 == "10"] <- "2"
  Expanded_list <- list()
  for (j in 0:7) {
    Expanded <- sum(CD8_expansion_id@meta.data$expansion_2 == "2" & CD8_expansion_id@meta.data$seurat_cluster == j, na.rm = TRUE)
    Total <- sum(CD8_expansion_id@meta.data$expansion_2 == "2" & CD8_expansion_id@meta.data$seurat_cluster == j, na.rm = TRUE) + 
      sum(CD8_expansion_id@meta.data$expansion_2 == "0" & CD8_expansion_id@meta.data$seurat_cluster == j, na.rm = TRUE) 
    Expanded_list[[j+1]] <- Expanded
  }
  
  # Create a dataframe with the results for the current patient
  patient_df <- data.frame(Patient_ID_SPC = i, Cluster = 0:7, Expanded = unlist(Expanded_list))
  
  # Append the patient dataframe to the overall results dataframe
  results_df <- rbind(results_df, patient_df)
}

# Use dcast to reshape the dataframe
results_reshaped <- dcast(results_df, Patient_ID_SPC ~ Cluster, value.var = "Expanded")

# Write the results_reshaped dataframe to a CSV file
write.csv(results_reshaped, file = "D:/DATA/Gastric cancer/CD8 TCR results counts per cluster.csv", row.names = FALSE)





# Fraction of expaned clones for each clusters ####
#CD8_expansion <- subset(x = CD8_rank_total, Type == 'Tumor' & Outcome != 'Meta' & Patient_ID_SPC != F)

patient_id_nometa = c(4, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 52, 53, 54, 56, 59, 60, 65, 66, 68, 71, 73, 75)

# Create an empty dataframe to store the results
results_df <- data.frame(Patient_ID_SPC = numeric(),
                         Cluster = numeric(),
                         Fraction = numeric())

for (i in patient_id_nometa){
  CD8_expansion_id <- subset(x = CD8_expansion, Patient_ID_SPC == i)
  CD8_expansion_id@meta.data$expansion_2 <- CD8_expansion_id@meta.data$expansion
  CD8_expansion_id@meta.data$expansion_2[CD8_expansion_id@meta.data$expansion_2 == "2"] <- "2"
  CD8_expansion_id@meta.data$expansion_2[CD8_expansion_id@meta.data$expansion_2 == "5"] <- "2"
  CD8_expansion_id@meta.data$expansion_2[CD8_expansion_id@meta.data$expansion_2 == "10"] <- "2"
  fraction_list <- list()
  for (j in 0:7) {
    Expanded <- sum(CD8_expansion_id@meta.data$expansion_2 == "2" & CD8_expansion_id@meta.data$seurat_cluster == j, na.rm = TRUE)
    Total <- sum(CD8_expansion_id@meta.data$expansion_2 == "2" & CD8_expansion_id@meta.data$seurat_cluster == j, na.rm = TRUE) + 
      sum(CD8_expansion_id@meta.data$expansion_2 == "0" & CD8_expansion_id@meta.data$seurat_cluster == j, na.rm = TRUE) 
    Fraction <- Expanded / Total * 100
    fraction_list[[j+1]] <- Fraction
  }
  
  # Create a dataframe with the results for the current patient
  patient_df <- data.frame(Patient_ID_SPC = i, Cluster = 0:7, Fraction = unlist(fraction_list))
  
  # Append the patient dataframe to the overall results dataframe
  results_df <- rbind(results_df, patient_df)
}


# Use dcast to reshape the dataframe
results_reshaped <- dcast(results_df, Patient_ID_SPC ~ Cluster, value.var = "Fraction")

# Write the results_reshaped dataframe to a CSV file
write.csv(results_reshaped, file = "D:/DATA/Gastric cancer/CD8 TCR results-2.csv", row.names = FALSE)



# Fraction of expaned clones for each clusters ####
#CD8_expansion <- subset(x = CD8_rank_total, Type == 'Tumor' & Outcome != 'Meta' & Patient_ID_SPC != F)

patient_id_nometa = c(4, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 52, 53, 54, 56, 59, 60, 65, 66, 68, 71, 73, 75)

# Create an empty dataframe to store the results
results_df <- data.frame(Patient_ID_SPC = numeric(),
                         Cluster = numeric(),
                         Fraction = numeric())

for (i in patient_id_nometa){
  CD8_expansion_id <- subset(x = CD8_expansion, Patient_ID_SPC == i)
  CD8_expansion_id@meta.data$expansion_2 <- CD8_expansion_id@meta.data$expansion
  CD8_expansion_id@meta.data$expansion_2[CD8_expansion_id@meta.data$expansion_2 == "2"] <- "0"
  CD8_expansion_id@meta.data$expansion_2[CD8_expansion_id@meta.data$expansion_2 == "5"] <- "2"
  CD8_expansion_id@meta.data$expansion_2[CD8_expansion_id@meta.data$expansion_2 == "10"] <- "2"
  fraction_list <- list()
  for (j in 0:7) {
    Expanded <- sum(CD8_expansion_id@meta.data$expansion_2 == "2" & CD8_expansion_id@meta.data$seurat_cluster == j, na.rm = TRUE)
    Total <- sum(CD8_expansion_id@meta.data$expansion_2 == "2" & CD8_expansion_id@meta.data$seurat_cluster == j, na.rm = TRUE) + 
      sum(CD8_expansion_id@meta.data$expansion_2 == "0" & CD8_expansion_id@meta.data$seurat_cluster == j, na.rm = TRUE) 
    Fraction <- Expanded / Total * 100
    fraction_list[[j+1]] <- Fraction
  }
  
  # Create a dataframe with the results for the current patient
  patient_df <- data.frame(Patient_ID_SPC = i, Cluster = 0:7, Fraction = unlist(fraction_list))
  
  # Append the patient dataframe to the overall results dataframe
  results_df <- rbind(results_df, patient_df)
}


# Use dcast to reshape the dataframe
results_reshaped <- dcast(results_df, Patient_ID_SPC ~ Cluster, value.var = "Fraction")

# Write the results_reshaped dataframe to a CSV file
write.csv(results_reshaped, file = "D:/DATA/Gastric cancer/CD8 TCR results-5.csv", row.names = FALSE)


# Fraction of expaned clones for each clusters ####
#CD8_expansion <- subset(x = CD8_rank_total, Type == 'Tumor' & Outcome != 'Meta' & Patient_ID_SPC != F)

patient_id_nometa = c(4, 7, 8, 15, 16, 18, 20, 22, 23, 27, 28, 29, 35, 36, 40, 41, 44, 50, 52, 53, 54, 56, 59, 60, 65, 66, 68, 71, 73, 75)

# Create an empty dataframe to store the results
results_df <- data.frame(Patient_ID_SPC = numeric(),
                         Cluster = numeric(),
                         Fraction = numeric())

for (i in patient_id_nometa){
  CD8_expansion_id <- subset(x = CD8_expansion, Patient_ID_SPC == i)
  CD8_expansion_id@meta.data$expansion_2 <- CD8_expansion_id@meta.data$expansion
  CD8_expansion_id@meta.data$expansion_2[CD8_expansion_id@meta.data$expansion_2 == "2"] <- "0"
  CD8_expansion_id@meta.data$expansion_2[CD8_expansion_id@meta.data$expansion_2 == "5"] <- "0"
  CD8_expansion_id@meta.data$expansion_2[CD8_expansion_id@meta.data$expansion_2 == "10"] <- "2"
  fraction_list <- list()
  for (j in 0:7) {
    Expanded <- sum(CD8_expansion_id@meta.data$expansion_2 == "2" & CD8_expansion_id@meta.data$seurat_cluster == j, na.rm = TRUE)
    Total <- sum(CD8_expansion_id@meta.data$expansion_2 == "2" & CD8_expansion_id@meta.data$seurat_cluster == j, na.rm = TRUE) + 
      sum(CD8_expansion_id@meta.data$expansion_2 == "0" & CD8_expansion_id@meta.data$seurat_cluster == j, na.rm = TRUE) 
    Fraction <- Expanded / Total * 100
    fraction_list[[j+1]] <- Fraction
  }
  
  # Create a dataframe with the results for the current patient
  patient_df <- data.frame(Patient_ID_SPC = i, Cluster = 0:7, Fraction = unlist(fraction_list))
  
  # Append the patient dataframe to the overall results dataframe
  results_df <- rbind(results_df, patient_df)
}


# Use dcast to reshape the dataframe
results_reshaped <- dcast(results_df, Patient_ID_SPC ~ Cluster, value.var = "Fraction")

# Write the results_reshaped dataframe to a CSV file
write.csv(results_reshaped, file = "D:/DATA/Gastric cancer/CD8 TCR results-10.csv", row.names = FALSE)





######################################## Shannon Diversity Index ####
CD8_shannon <- CD8_rank_total

view(CD8_rank_total@meta.data)

CD8_shannon <- subset(x = CD8_rank_total, Type == 'Tumor' & Outcome != 'Meta' & Patient_ID_SPC != F)


output_list <- list()

for (k in 0:7) {
  df <- paste0("CD8_shannon_",k)
  df <- subset(CD8_shannon@meta.data, seurat_clusters == k)
  
  count_df <- df %>% 
    group_by(concatenated) %>% 
    summarise(count = n()) %>% 
    ungroup()
  
  count_df <- count_df[complete.cases(count_df[, "concatenated"]), ]
  count_df <- count_df[count_df$concatenated != "", ]
  
  count_values <- count_df[,"count"]
  
  shannon <- diversity(count_values,"shannon")
  output_list[[k+1]] <- list(k = k, shannon = shannon)
}


output_df <- do.call(rbind, output_list)

ggplot(output_df, aes(x = factor(k), y = output)) + 
  geom_bar(stat = "identity") +
  labs(x = "k value", y = "Output value", title = "Output vs. k value")

















