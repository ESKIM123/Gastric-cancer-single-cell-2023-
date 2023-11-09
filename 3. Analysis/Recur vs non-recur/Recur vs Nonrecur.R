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



# Get data ----
cancer.subset <- readRDS(file="D:/DATA/Gastric cancer/immune_cell_clinical_data_modify.rds")
DimPlot(cancer.subset, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Initial")

view(cancer.subset@meta.data)

# Subset only Stage III with tumor and recurrence data

cancer_III_subset <- subset(x = cancer.subset, Type == 'Tumor' & Patient_ID_SPC != FALSE & Stage == 'Stage III' & Outcome_new != 'NA')
cancer_III_subset_Nonrecur <- subset(x = cancer.subset, Type == 'Tumor' & Patient_ID_SPC != FALSE & Stage == 'Stage III' & Outcome_new == 'Nonrecur')
cancer_III_subset_Recur <- subset(x = cancer.subset, Type == 'Tumor' & Patient_ID_SPC != FALSE & Stage == 'Stage III' & Outcome_new == 'Recur')


# plot difference
DimPlot(cancer_III_subset, reduction = 'umap', group.by = "Outcome_new", raster=FALSE, label = F, cols = c("#6495ED", "lightgrey"))+
  NoAxes()+NoLegend()+ ggtitle("")

DimPlot(cancer_III_subset_Nonrecur, reduction = 'umap', group.by = "Outcome_new", raster=FALSE, label = F, cols = "#6495ED")+
  NoAxes()+NoLegend()+ ggtitle("")
DimPlot(cancer_III_subset_Recur, reduction = 'umap', group.by = "Outcome_new", raster=FALSE, label = F, cols = "Red")+
  NoAxes()+NoLegend()+ ggtitle("")

cancer_III_subset_Nonrecur
cancer_III_subset_Recur

# Cell frequency difference ------------------------------------------------------------------------------------


# SCpubr::do_BarPlot(sample = cancer_III_subset, 
#                    group.by = "clusters", 
#                    split.by = "Outcome_new",
#                    position = "stack",legend.position = "bottom",
#                    plot.title = "", xlab = "",  ylab = "Number of cells per cluster", font.size = 11)
# 
# SCpubr::do_BarPlot(sample = cancer_III_subset, 
#                    group.by = "clusters", 
#                    split.by = "Outcome_new",
#                    position = "Fill",legend.position = "bottom",
#                    plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)



dittoBarPlot(cancer_III_subset, "clusters", group.by = "Outcome_new", scale = "count", main = "")+
  NoLegend()+ xlab(NULL)+ ylab(NULL)+
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))

dittoBarPlot(cancer_III_subset, "clusters", group.by = "Outcome_new", scale = "percent", main = "")+ 
  NoLegend()+ xlab(NULL)+ ylab(NULL)+
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))


# Cell frequency difference ------------------------------------------------------------------------------------

# 
# SCpubr::do_BarPlot(sample = cancer_III_subset_Nonrecur, 
#                    group.by = "clusters", 
#                    split.by = "Patient_ID_SPC",
#                    position = "Fill",legend.position = "bottom",
#                    plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)
# 
# 
# SCpubr::do_BarPlot(sample = cancer_III_subset_Recur, 
#                    group.by = "clusters", 
#                    split.by = "Patient_ID_SPC",
#                    position = "Fill",legend.position = "bottom",
#                    plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)



dittoBarPlot(cancer_III_subset_Nonrecur, "clusters", group.by = "Patient_ID_SPC", scale = "percent", main = "")+ 
  NoLegend()+ xlab(NULL)+ ylab(NULL)+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))

dittoBarPlot(cancer_III_subset_Recur, "clusters", group.by = "Patient_ID_SPC", scale = "percent", main = "")+ 
  NoLegend()+ xlab(NULL)+ ylab(NULL)+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))


# Cell frequency difference ------------------------------------------------------------------------------------

# SCpubr::do_BarPlot(sample = cancer_III_subset_Nonrecur, 
#                    group.by = "clusters", 
#                    split.by = "Patient_ID_SPC",
#                    position = "Fill",legend.position = "bottom",
#                    plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)
# 
# SCpubr::do_BarPlot(sample = cancer_III_subset_Recur, 
#                    group.by = "clusters", 
#                    split.by = "Patient_ID_SPC",
#                    position = "Fill",legend.position = "bottom",
#                    plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)

dittoBarPlot(cancer_III_subset_Nonrecur, "clusters", group.by = "Patient_ID_SPC", scale = "percent", main = "")+ 
  NoLegend()+ xlab(NULL)+ ylab(NULL)+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))

dittoBarPlot(cancer_III_subset_Recur, "clusters", group.by = "Patient_ID_SPC", scale = "percent", main = "")+ 
  NoLegend()+ xlab(NULL)+ ylab(NULL)+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))




dittoBarPlot(cancer_III_subset_Nonrecur, "clusters", group.by = "Patient_ID_SPC", scale = "percent", main = "")

##################################### Stats #########################################
##################################### Dataframe for Outcome
md.cancer_III_subset<- cancer_III_subset@meta.data %>% as.data.table

length(which(md.cancer_III_subset$Patient_ID==4 & md.cancer_III_subset$Type=='Tumor'& md.cancer_III_subset$seurat_clusters==1))

# Vector of sample numbers
sample_numbers <- c(4, 7, 8, 16, 20, 22, 23, 27, 28, 29, 36, 41, 50, 53, 54, 62, 65, 66)

# Vector of cluster names
cluster_vec <- unique(md.cancer_III_subset$clusters)
type_vec <- c("Nonrecur", "Recur")

unique(cancer_III_subset@meta.data$Patient_ID_SPC)

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
      subset_length <- length(which(md.cancer_III_subset$Patient_ID_SPC == number & md.cancer_III_subset$Outcome_new == type & md.cancer_III_subset$clusters == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.cancer_III_subset$Patient_ID_SPC == number & md.cancer_III_subset$Outcome_new == type))
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

# ################################################################################## Graph
# 
# # Create an empty list to store the plots
# plot_list <- list()
# 
# for (cluster in cluster_vec) {
#   # Subset the data for each type_vec
#   subset_df <- output_df %>% subset(cluster_vec == cluster)
#   Nonrecur_df <- subset_df %>% subset(type_vec == "Nonrecur")
#   Recur_df <- subset_df %>% subset(type_vec == "Recur")
#   
#   # Perform the unpaired t-test
#   ttest <- t.test(Nonrecur_df$subset_fraction, Recur_df$subset_fraction)
#   wtest <- wilcox.test(Nonrecur_df$subset_fraction, Recur_df$subset_fraction)
#   
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
# grid.arrange(grobs = plot_list, nrow = 3)
# 
# ############## Stat add
# 
# for (cluster in cluster_vec) {
#   # Subset the data for each type_vec
#   subset_df <- output_df %>% subset(cluster_vec == cluster)
#   Nonrecur_df <- subset_df %>% subset(type_vec == "Nonrecur")
#   Recur_df <- subset_df %>% subset(type_vec == "Recur")
#   
#   # Perform the unpaired t-test
#   ttest <- t.test(Nonrecur_df$subset_fraction, Recur_df$subset_fraction)
#   wtest <- wilcox.test(Nonrecur_df$subset_fraction, Recur_df$subset_fraction)
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
# ggsave("output_graph_unpaired.png", plot = grid.arrange(grobs = plot_list, nrow = 3), width = 13, height = 5, units = "in")
# 



################################################################################## Graph

# Create an empty list to store the plots
plot_list <- list()

for (cluster in cluster_vec) {
  # Subset the data for each type_vec
  subset_df <- output_df %>% subset(cluster_vec == cluster)
  Nonrecur_df <- subset_df %>% subset(type_vec == "Nonrecur")
  Recur_df <- subset_df %>% subset(type_vec == "Recur")
  
  # Perform the unpaired t-test
  ttest <- t.test(Nonrecur_df$subset_fraction, Recur_df$subset_fraction)
  wtest <- wilcox.test(Nonrecur_df$subset_fraction, Recur_df$subset_fraction)
  
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction)) +
    geom_boxplot(color = "black", outlier.shape = NA, fill = c("#E1D4BB", "#537188"), alpha=0.8)+
    geom_point(shape = 21, color = "black", alpha=0.8, aes(fill=type_vec)) + 
    scale_fill_manual(values = c("grey", "brown")) +
    scale_y_continuous(labels = function(x) paste0(x * 100))+
    geom_line(aes(group = sample_numbers), alpha=0.5)+
    ggtitle(cluster) + 
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
  Nonrecur_df <- subset_df %>% subset(type_vec == "Nonrecur")
  Recur_df <- subset_df %>% subset(type_vec == "Recur")
  
  # Perform the unpaired t-test
  ttest <- t.test(Nonrecur_df$subset_fraction, Recur_df$subset_fraction)
  wtest <- wilcox.test(Nonrecur_df$subset_fraction, Recur_df$subset_fraction)
  
  g1<-ggplot(subset_df, aes(x = type_vec, y = subset_fraction)) +
    geom_boxplot(color = "black", outlier.shape = NA, fill = c("grey", "brown"), alpha=0.5)+
    geom_point(shape = 21, color = "black", alpha=0.8, aes(fill=type_vec)) + 
    scale_fill_manual(values = c("grey", "brown")) +
    geom_line(aes(group = sample_numbers), alpha=0.5)+
    ggtitle(cluster) + 
    xlab("") + ylab("") +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          axis.line = element_line(color = "black", size = 0.2),
          axis.ticks = element_line(size = 0.2),
          axis.ticks.length = unit(0.1, "cm"),
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




################################################################################## DEGS #####################################
DEG_cancer_III_subset <- cancer_III_subset
DEG_cancer_III_subset[["Total"]]  <- "All"


# pseudo-bulk workflow -----------------
# Acquiring necessary metrics for aggregation across cells in a sample
# 1. counts matrix - sample level
# counts aggregate to sample level

#View(DEG_cancer_III_subset@meta.data)
DEG_cancer_III_subset$samples <- paste0(DEG_cancer_III_subset$Outcome_new, DEG_cancer_III_subset$Patient_ID_SPC)

DEG_cancer_III_subset
DefaultAssay(object = DEG_cancer_III_subset) <- "RNA"

cts <- AggregateExpression(DEG_cancer_III_subset, 
                           group.by = "samples",
                           assays = 'RNA',
                           slot = "counts",
                           return.seurat = FALSE)

cts <- cts$RNA

# transpose
cts.t <- t(cts)

# convert to data.frame
cts.convert <- as.data.frame(cts)

# Let's run DE analysis with B cells
# 1. Get counts matrix
counts_all <- cts.convert

# 2. generate sample level metadata
colData <- data.frame(samples = colnames(counts_all))

colData <- colData %>%
  mutate(condition = ifelse(grepl('Nonrecur', samples), 'Nonrecur', 'Recur')) %>%
  column_to_rownames(var = 'samples')

# perform DESeq2 -------- All
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_all,
                              colData = colData,
                              design = ~ condition)

# filter
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

# run DESeq2
dds <- DESeq(dds)

# Check the coefficients for the comparison
resultsNames(dds)

# Generate results object
res <- results(dds, name = "condition_Recur_vs_Nonrecur")
res

res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  arrange(padj)

# Check results output
res_tbl 


# Graph
condition.diffgenes <- res_tbl
condition.diffgenes$LOGPvalue <- -log10(condition.diffgenes$padj)
condition.diffgenes$Difference <- condition.diffgenes$log2FoldChange


colnames(condition.diffgenes)

# SCpubr::do_VolcanoPlot(sample = DEG_cancer_III_subset, de_genes = condition.diffgenes, n_genes = 15)
# head(condition.diffgenes, n = 15)


ggplot(condition.diffgenes, aes(x = Difference, y = LOGPvalue)) + 
  geom_point(aes(color = ifelse(LOGPvalue > 1.3 & Difference > 0.25, "blue", 
                                ifelse(LOGPvalue > 1.3 & Difference < -0.25, "grey", "red"))), alpha = 0.8) + 
  scale_color_manual(values = c("red", "blue", "grey")) +
  theme_classic() + 
  labs(x = "Fold change (Log2)", y = "-log(P-value)", title = "") + guides(color = FALSE)+
  geom_vline(xintercept = 0.25, linetype = "dashed") +
  geom_vline(xintercept = -0.25, linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  coord_cartesian(clip = 'off')+
  geom_text_repel(
    data = subset(condition.diffgenes, LOGPvalue > 1.3 & abs(Difference) > 1),
    aes(label = gene),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15))

VlnPlot(DEG_cancer_III_subset, features = "IGKV6-21", split.by = "Outcome_new", pt.size = 0)
# DefaultAssay(object = DEG_cancer_III_subset) <- "integrated"



# 오른쪽이 Non-recur

# perform DESeq2 -------- CD8 T cell
# pseudo-bulk workflow -----------------
# Acquiring necessary metrics for aggregation across cells in a sample
# 1. counts matrix - sample level
# counts aggregate to sample level

View(DEG_cancer_III_subset@meta.data)
DEG_cancer_III_subset$samples <- paste0(DEG_cancer_III_subset$Outcome_new, DEG_cancer_III_subset$Patient_ID_SPC)

DEG_cancer_III_subset
DefaultAssay(object = DEG_cancer_III_subset) <- "RNA"

cts <- AggregateExpression(DEG_cancer_III_subset, 
                           group.by = c("clusters", "samples"),
                           assays = 'RNA',
                           slot = "counts",
                           return.seurat = FALSE)
cts <- cts$RNA

# transpose
cts.t <- t(cts)

# convert to data.frame
cts.t <- as.data.frame(cts.t)

# get values where to split
splitRows <- gsub('_.*', '', rownames(cts.t))

# split data.frame
cts.split <- split.data.frame(cts.t,
                              f = factor(splitRows))

# fix colnames and transpose

cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
  
})

# Let's run DE analysis with CD8 T cells
# 1. Get counts matrix
counts_CD8tcell <- cts.split.modified$`CD8 T cells`


# 2. generate sample level metadata
colData <- data.frame(samples = colnames(counts_CD8tcell))

colData <- colData %>%
  mutate(condition = ifelse(grepl('Nonrecur', samples), 'Nonrecur', 'Recur')) %>%
  column_to_rownames(var = 'samples')


# perform DESeq2 --------
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_CD8tcell,
                              colData = colData,
                              design = ~ condition)

# filter
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

# run DESeq2
dds <- DESeq(dds)
plotDispEsts(dds)

# Check the coefficients for the comparison
resultsNames(dds)

# Generate results object
res <- results(dds, name = "condition_Recur_vs_Nonrecur")
res


res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  arrange(padj)

# Check results output
res_tbl 


# Graph
condition.diffgenes <- res_tbl
condition.diffgenes$LOGPvalue <- -log10(condition.diffgenes$padj)
condition.diffgenes$Difference <- condition.diffgenes$log2FoldChange

colnames(condition.diffgenes)

ggplot(condition.diffgenes, aes(x = Difference, y = LOGPvalue)) + 
  geom_point(aes(color = ifelse(LOGPvalue > 1.3 & Difference > 0.25, "blue", 
                                ifelse(LOGPvalue > 1.3 & Difference < -0.25, "grey", "red"))), alpha = 0.8) + 
  scale_color_manual(values = c("red", "blue", "grey")) +
  theme_classic() + 
  labs(x = "Fold change (Log2)", y = "-log(P-value)", title = "") + guides(color = FALSE)+
  geom_vline(xintercept = 0.25, linetype = "dashed") +
  geom_vline(xintercept = -0.25, linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  coord_cartesian(clip = 'off')+
  geom_text_repel(
    data = subset(condition.diffgenes, LOGPvalue > 1.3 & abs(Difference) > 1),
    aes(label = gene),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15))



# Let's run DE analysis with CD4 T cells-------------------------------------------------------------------------
# 1. Get counts matrix
counts_CD4tcell <- cts.split.modified$`CD4 T cells`


# 2. generate sample level metadata
colData <- data.frame(samples = colnames(counts_CD4tcell))

colData <- colData %>%
  mutate(condition = ifelse(grepl('Nonrecur', samples), 'Nonrecur', 'Recur')) %>%
  column_to_rownames(var = 'samples')


# perform DESeq2 --------
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_CD4tcell,
                              colData = colData,
                              design = ~ condition)

# filter
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

# run DESeq2
dds <- DESeq(dds)
plotDispEsts(dds)

# Check the coefficients for the comparison
resultsNames(dds)

# Generate results object
res <- results(dds, name = "condition_Recur_vs_Nonrecur")
res


res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  arrange(padj)

# Check results output
res_tbl 


# Graph
condition.diffgenes <- res_tbl
condition.diffgenes$LOGPvalue <- -log10(condition.diffgenes$padj)
condition.diffgenes$Difference <- condition.diffgenes$log2FoldChange

colnames(condition.diffgenes)

ggplot(condition.diffgenes, aes(x = Difference, y = LOGPvalue)) + 
  geom_point(aes(color = ifelse(LOGPvalue > 1.3 & Difference > 0.25, "blue", 
                                ifelse(LOGPvalue > 1.3 & Difference < -0.25, "grey", "red"))), alpha = 0.8) + 
  scale_color_manual(values = c("red", "blue", "grey")) +
  theme_classic() + 
  labs(x = "Fold change (Log2)", y = "-log(P-value)", title = "") + guides(color = FALSE)+
  geom_vline(xintercept = 0.25, linetype = "dashed") +
  geom_vline(xintercept = -0.25, linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  coord_cartesian(clip = 'off')+
  geom_text_repel(
    data = subset(condition.diffgenes, LOGPvalue > 1.3 & abs(Difference) > 1),
    aes(label = gene),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15))

# Let's run DE analysis with B cells-------------------------------------------------------------------------
# 1. Get counts matrix
counts_Bcell <- cts.split.modified$`B cells`


# 2. generate sample level metadata
colData <- data.frame(samples = colnames(counts_Bcell))

colData <- colData %>%
  mutate(condition = ifelse(grepl('Nonrecur', samples), 'Nonrecur', 'Recur')) %>%
  column_to_rownames(var = 'samples')


# perform DESeq2 --------
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_Bcell,
                              colData = colData,
                              design = ~ condition)

# filter
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

# run DESeq2
dds <- DESeq(dds)
plotDispEsts(dds)

# Check the coefficients for the comparison
resultsNames(dds)

# Generate results object
res <- results(dds, name = "condition_Recur_vs_Nonrecur")
res


res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  arrange(padj)

# Check results output
res_tbl 


# Graph
condition.diffgenes <- res_tbl
condition.diffgenes$LOGPvalue <- -log10(condition.diffgenes$padj)
condition.diffgenes$Difference <- condition.diffgenes$log2FoldChange

colnames(condition.diffgenes)

ggplot(condition.diffgenes, aes(x = Difference, y = LOGPvalue)) + 
  geom_point(aes(color = ifelse(LOGPvalue > 1.3 & Difference > 0.25, "blue", 
                                ifelse(LOGPvalue > 1.3 & Difference < -0.25, "grey", "red"))), alpha = 0.8) + 
  scale_color_manual(values = c("red", "blue", "grey")) +
  theme_classic() + 
  labs(x = "Fold change (Log2)", y = "-log(P-value)", title = "") + guides(color = FALSE)+
  geom_vline(xintercept = 0.25, linetype = "dashed") +
  geom_vline(xintercept = -0.25, linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  coord_cartesian(clip = 'off')+
  geom_text_repel(
    data = subset(condition.diffgenes, LOGPvalue > 1.3 & abs(Difference) > 1),
    aes(label = gene),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15))


# Let's run DE analysis with Plasma cells-------------------------------------------------------------------------
# 1. Get counts matrix
counts_Plasmacell <- cts.split.modified$`Plasma cells`


# 2. generate sample level metadata
colData <- data.frame(samples = colnames(counts_Plasmacell))

colData <- colData %>%
  mutate(condition = ifelse(grepl('Nonrecur', samples), 'Nonrecur', 'Recur')) %>%
  column_to_rownames(var = 'samples')


# perform DESeq2 --------
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_Plasmacell,
                              colData = colData,
                              design = ~ condition)

# filter
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

# run DESeq2
dds <- DESeq(dds)
plotDispEsts(dds)

# Check the coefficients for the comparison
resultsNames(dds)

# Generate results object
res <- results(dds, name = "condition_Recur_vs_Nonrecur")
res


res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  arrange(padj)

# Check results output
res_tbl 


# Graph
condition.diffgenes <- res_tbl
condition.diffgenes$LOGPvalue <- -log10(condition.diffgenes$padj)
condition.diffgenes$Difference <- condition.diffgenes$log2FoldChange

colnames(condition.diffgenes)

ggplot(condition.diffgenes, aes(x = Difference, y = LOGPvalue)) + 
  geom_point(aes(color = ifelse(LOGPvalue > 1.3 & Difference > 0.25, "blue", 
                                ifelse(LOGPvalue > 1.3 & Difference < -0.25, "grey", "red"))), alpha = 0.8) + 
  scale_color_manual(values = c("red", "blue", "grey")) +
  theme_classic() + 
  labs(x = "Fold change (Log2)", y = "-log(P-value)", title = "") + guides(color = FALSE)+
  geom_vline(xintercept = 0.25, linetype = "dashed") +
  geom_vline(xintercept = -0.25, linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  coord_cartesian(clip = 'off')+
  geom_text_repel(
    data = subset(condition.diffgenes, LOGPvalue > 1.3 & abs(Difference) > 1),
    aes(label = gene),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15))




# Let's run DE analysis with NK cells-------------------------------------------------------------------------
# 1. Get counts matrix
counts_NKcell <- cts.split.modified$`NK cells`


# 2. generate sample level metadata
colData <- data.frame(samples = colnames(counts_NKcell))

colData <- colData %>%
  mutate(condition = ifelse(grepl('Nonrecur', samples), 'Nonrecur', 'Recur')) %>%
  column_to_rownames(var = 'samples')


# perform DESeq2 --------
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_NKcell,
                              colData = colData,
                              design = ~ condition)

# filter
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

# run DESeq2
dds <- DESeq(dds)
plotDispEsts(dds)

# Check the coefficients for the comparison
resultsNames(dds)

# Generate results object
res <- results(dds, name = "condition_Recur_vs_Nonrecur")
res


res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  arrange(padj)

# Check results output
res_tbl 


# Graph
condition.diffgenes <- res_tbl
condition.diffgenes$LOGPvalue <- -log10(condition.diffgenes$padj)
condition.diffgenes$Difference <- condition.diffgenes$log2FoldChange

colnames(condition.diffgenes)

ggplot(condition.diffgenes, aes(x = Difference, y = LOGPvalue)) + 
  geom_point(aes(color = ifelse(LOGPvalue > 1.3 & Difference > 0.25, "blue", 
                                ifelse(LOGPvalue > 1.3 & Difference < -0.25, "grey", "red"))), alpha = 0.8) + 
  scale_color_manual(values = c("red", "blue", "grey")) +
  theme_classic() + 
  labs(x = "Fold change (Log2)", y = "-log(P-value)", title = "") + guides(color = FALSE)+
  geom_vline(xintercept = 0.25, linetype = "dashed") +
  geom_vline(xintercept = -0.25, linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  coord_cartesian(clip = 'off')+
  geom_text_repel(
    data = subset(condition.diffgenes, LOGPvalue > 1.3 & abs(Difference) > 1),
    aes(label = gene),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15))



# Let's run DE analysis with Treg cells-------------------------------------------------------------------------
# 1. Get counts matrix
counts_Tregcell <- cts.split.modified$`Regulatory CD4 T cells`


# 2. generate sample level metadata
colData <- data.frame(samples = colnames(counts_Tregcell))

colData <- colData %>%
  mutate(condition = ifelse(grepl('Nonrecur', samples), 'Nonrecur', 'Recur')) %>%
  column_to_rownames(var = 'samples')


# perform DESeq2 --------
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_Tregcell,
                              colData = colData,
                              design = ~ condition)

# filter
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

# run DESeq2
dds <- DESeq(dds)
plotDispEsts(dds)

# Check the coefficients for the comparison
resultsNames(dds)

# Generate results object
res <- results(dds, name = "condition_Recur_vs_Nonrecur")
res


res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  arrange(padj)

# Check results output
res_tbl 


# Graph
condition.diffgenes <- res_tbl
condition.diffgenes$LOGPvalue <- -log10(condition.diffgenes$padj)
condition.diffgenes$Difference <- condition.diffgenes$log2FoldChange

colnames(condition.diffgenes)

ggplot(condition.diffgenes, aes(x = Difference, y = LOGPvalue)) + 
  geom_point(aes(color = ifelse(LOGPvalue > 1.3 & Difference > 0.25, "blue", 
                                ifelse(LOGPvalue > 1.3 & Difference < -0.25, "grey", "red"))), alpha = 0.8) + 
  scale_color_manual(values = c("red", "blue", "grey")) +
  theme_classic() + 
  labs(x = "Fold change (Log2)", y = "-log(P-value)", title = "") + guides(color = FALSE)+
  geom_vline(xintercept = 0.25, linetype = "dashed") +
  geom_vline(xintercept = -0.25, linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  coord_cartesian(clip = 'off')+
  geom_text_repel(
    data = subset(condition.diffgenes, LOGPvalue > 1.3 & abs(Difference) > 1),
    aes(label = gene),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15))



# Let's run DE analysis with Myeloid cells-------------------------------------------------------------------------
# 1. Get counts matrix
counts_Myeloidcell <- cts.split.modified$`Myeloid cells`


# 2. generate sample level metadata
colData <- data.frame(samples = colnames(counts_Myeloidcell))

colData <- colData %>%
  mutate(condition = ifelse(grepl('Nonrecur', samples), 'Nonrecur', 'Recur')) %>%
  column_to_rownames(var = 'samples')


# perform DESeq2 --------
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_Myeloidcell,
                              colData = colData,
                              design = ~ condition)

# filter
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

# run DESeq2
dds <- DESeq(dds)
plotDispEsts(dds)

# Check the coefficients for the comparison
resultsNames(dds)

# Generate results object
res <- results(dds, name = "condition_Recur_vs_Nonrecur")
res


res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  arrange(padj)

# Check results output
res_tbl 


# Graph
condition.diffgenes <- res_tbl
condition.diffgenes$LOGPvalue <- -log10(condition.diffgenes$padj)
condition.diffgenes$Difference <- condition.diffgenes$log2FoldChange

colnames(condition.diffgenes)

ggplot(condition.diffgenes, aes(x = Difference, y = LOGPvalue)) + 
  geom_point(aes(color = ifelse(LOGPvalue > 1.3 & Difference > 0.25, "blue", 
                                ifelse(LOGPvalue > 1.3 & Difference < -0.25, "grey", "red"))), alpha = 0.8) + 
  scale_color_manual(values = c("red", "blue", "grey")) +
  theme_classic() + 
  labs(x = "Fold change (Log2)", y = "-log(P-value)", title = "") + guides(color = FALSE)+
  geom_vline(xintercept = 0.25, linetype = "dashed") +
  geom_vline(xintercept = -0.25, linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  coord_cartesian(clip = 'off')+
  geom_text_repel(
    data = subset(condition.diffgenes, LOGPvalue > 1.3 & abs(Difference) > 1),
    aes(label = gene),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15))


# Let's run DE analysis with GD T cells-------------------------------------------------------------------------
# 1. Get counts matrix
counts_GDTcell <- cts.split.modified$`Gamma-delta T cells`


# 2. generate sample level metadata
colData <- data.frame(samples = colnames(counts_GDTcell))

colData <- colData %>%
  mutate(condition = ifelse(grepl('Nonrecur', samples), 'Nonrecur', 'Recur')) %>%
  column_to_rownames(var = 'samples')


# perform DESeq2 --------
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_GDTcell,
                              colData = colData,
                              design = ~ condition)

# filter
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

# run DESeq2
dds <- DESeq(dds)
plotDispEsts(dds)

# Check the coefficients for the comparison
resultsNames(dds)

# Generate results object
res <- results(dds, name = "condition_Recur_vs_Nonrecur")
res


res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  arrange(padj)

# Check results output
res_tbl 


# Graph
condition.diffgenes <- res_tbl
condition.diffgenes$LOGPvalue <- -log10(condition.diffgenes$padj)
condition.diffgenes$Difference <- condition.diffgenes$log2FoldChange

colnames(condition.diffgenes)

ggplot(condition.diffgenes, aes(x = Difference, y = LOGPvalue)) + 
  geom_point(aes(color = ifelse(LOGPvalue > 1.3 & Difference > 0.25, "blue", 
                                ifelse(LOGPvalue > 1.3 & Difference < -0.25, "grey", "red"))), alpha = 0.8) + 
  scale_color_manual(values = c("red", "blue", "grey")) +
  theme_classic() + 
  labs(x = "Fold change (Log2)", y = "-log(P-value)", title = "") + guides(color = FALSE)+
  geom_vline(xintercept = 0.25, linetype = "dashed") +
  geom_vline(xintercept = -0.25, linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  coord_cartesian(clip = 'off')+
  geom_text_repel(
    data = subset(condition.diffgenes, LOGPvalue > 1.3 & abs(Difference) > 1),
    aes(label = gene),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15))




























































































################################################################################## DEGS #####################################

DefaultAssay(object = DEG_cancer_III_subset) <- "integrated"

DEG_cancer_III_subset <- cancer_III_subset
Idents(DEG_cancer_III_subset) <- "Outcome_new"
condition.diffgenes <- FindMarkers(DEG_cancer_III_subset, ident.1 = "Nonrecur", ident.2="Recur", min.pct=0, logfc.threshold=-Inf)
condition.diffgenes

###### 
condition.diffgenes <-rownames_to_column(condition.diffgenes, var = "Gene")
condition.diffgenes$LOGPvalue <- -log10(condition.diffgenes$p_val_adj)
condition.diffgenes$Difference <- condition.diffgenes$avg_log2FC

colnames(condition.diffgenes)

SCpubr::do_VolcanoPlot(sample = DEG_cancer_III_subset, de_genes = condition.diffgenes, n_genes = 15)
head(condition.diffgenes, n = 15)


ggplot(condition.diffgenes, aes(x = Difference, y = LOGPvalue)) + 
  geom_point(aes(color = ifelse(LOGPvalue > 1.3 & Difference > 0.25, "blue", 
                                ifelse(LOGPvalue > 1.3 & Difference < -0.25, "grey", "red"))), alpha = 0.8) + 
  scale_color_manual(values = c("red", "blue", "grey")) +
  theme_classic() + 
  labs(x = "Fold change (Log2)", y = "-log(P-value)", title = "") + guides(color = FALSE)+
  geom_vline(xintercept = 0.25, linetype = "dashed") +
  geom_vline(xintercept = -0.25, linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  coord_cartesian(clip = 'off')+
  geom_text_repel(
    data = subset(condition.diffgenes, LOGPvalue > 1.3 & abs(Difference) > 1),
    aes(label = Gene),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15))

VlnPlot(DEG_cancer_III_subset, features = "FBLN1", split.by = "Outcome_new", pt.size = 0)
# 오른쪽이 Non-recur

# GO analysis--------------------

# Set threshold for p-value and fold change
pvalue_threshold <- 0.05
fold_change_threshold <- 0.25

# Filter dataframe for significant genes
condition.diffgenes_significant_genes <- condition.diffgenes[condition.diffgenes$LOGPvalue > -log10(pvalue_threshold) & abs(condition.diffgenes$Difference) > fold_change_threshold, ]
condition.diffgenes_significant_genes_nonrecur <- condition.diffgenes[condition.diffgenes$LOGPvalue > -log10(pvalue_threshold) & condition.diffgenes$Difference > fold_change_threshold, ]
condition.diffgenes_significant_genes_recur <- condition.diffgenes[condition.diffgenes$LOGPvalue > -log10(pvalue_threshold) & condition.diffgenes$Difference < -fold_change_threshold, ]

sum(condition.diffgenes_significant_genes$Difference > 0.25)
sum(condition.diffgenes_significant_genes$Difference < -0.25)

write.csv(condition.diffgenes_significant_genes_nonrecur, file = "D:/DATA/condition.diffgenes_significant_genes_nonrecur.csv", row.names = FALSE)
write.csv(condition.diffgenes_significant_genes_recur, file = "D:/DATA/condition.diffgenes_significant_genes_recur.csv", row.names = FALSE)



################################################################################### GO ####
# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db" 

keytypes(org.Hs.eg.db)

BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

# reading in data 
df <- condition.diffgenes_significant_genes
df$log2FoldChange <- df$Difference

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$Gene

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")


gse_BP <- gseGO(geneList=gene_list, 
                ont ="BP", 
                keyType = "SYMBOL", 
                nPerm = 10000, 
                minGSSize = 3, 
                maxGSSize = 800, 
                pvalueCutoff = 0.05, 
                verbose = TRUE, 
                OrgDb = organism, 
                pAdjustMethod = "none")


write.csv(gse@result, file = "D:/DATA/condition.diffgenes_significant_genes_total_go_result.csv", row.names = FALSE)

require(DOSE)


gse_BP_GO <- gse_BP
view(gse_BP_GO@result)

######################## Draw figure (dotplot)

dotplot(gse_BP_GO, showCategory=10, split=".sign", x="NES", decreasing = TRUE, font.size = 6, label_format = 30)+
  scale_x_continuous(limit = c(-3, 3), breaks=seq(-3,3,1)) +
  geom_vline(xintercept = 0, linetype="dotted", color = "#444444", size=1) +
  theme( 
    legend.text=element_text(size=5),
    legend.position = "bottom",
    legend.box = "vertical",
    text = element_text(size=6))

dotplot(gse_BP, showCategory=9, split=".sign", x="NES", decreasing = TRUE, font.size = 6, label_format = 30)






























