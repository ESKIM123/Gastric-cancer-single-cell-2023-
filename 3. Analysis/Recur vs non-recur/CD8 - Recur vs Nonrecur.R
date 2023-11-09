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
library(DESeq2)


# Get data ----
CD8_subset_anno <- readRDS(file="D:/DATA/Gastric cancer/GS_subsets/CD8_subset.rds")
CD8_cancer.subset <- CD8_subset_anno
DimPlot(CD8_cancer.subset, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Initial")

view(CD8_cancer.subset@meta.data)

# Subset only Stage III with tumor and recurrence data

CD8_cancer_III_subset <- subset(x = CD8_cancer.subset, Type == 'Tumor' & Patient_ID_SPC != FALSE & Stage == 'Stage III' & Outcome_new != 'NA')
CD8_cancer_III_subset_Nonrecur <- subset(x = CD8_cancer.subset, Type == 'Tumor' & Patient_ID_SPC != FALSE & Stage == 'Stage III' & Outcome_new == 'Nonrecur')
CD8_cancer_III_subset_Recur <- subset(x = CD8_cancer.subset, Type == 'Tumor' & Patient_ID_SPC != FALSE & Stage == 'Stage III' & Outcome_new == 'Recur')


# plot difference
DimPlot(CD8_cancer_III_subset, reduction = 'umap', group.by = "Outcome_new", raster=FALSE, label = F, cols = c("#6495ED", "lightgrey"))+
  NoAxes()+NoLegend()+ ggtitle("")

DimPlot(CD8_cancer_III_subset_Nonrecur, reduction = 'umap', group.by = "Outcome_new", raster=FALSE, label = F, cols = "#6495ED")+
  NoAxes()+NoLegend()+ ggtitle("")
DimPlot(CD8_cancer_III_subset_Recur, reduction = 'umap', group.by = "Outcome_new", raster=FALSE, label = F, cols = "Red")+
  NoAxes()+NoLegend()+ ggtitle("")

CD8_cancer_III_subset_Nonrecur
CD8_cancer_III_subset_Recur

# Cell frequency difference ------------------------------------------------------------------------------------


# SCpubr::do_BarPlot(sample = CD8_cancer_III_subset, 
#                    group.by = "clusters", 
#                    split.by = "Outcome_new",
#                    position = "stack",legend.position = "bottom",
#                    plot.title = "", xlab = "",  ylab = "Number of cells per cluster", font.size = 11)
# 
# SCpubr::do_BarPlot(sample = CD8_cancer_III_subset, 
#                    group.by = "clusters", 
#                    split.by = "Outcome_new",
#                    position = "Fill",legend.position = "bottom",
#                    plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)



dittoBarPlot(CD8_cancer_III_subset, "seurat_clusters", group.by = "Outcome_new", scale = "count", main = "")+
  NoLegend()+ xlab(NULL)+ ylab(NULL)+
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))

dittoBarPlot(CD8_cancer_III_subset, "seurat_clusters", group.by = "Outcome_new", scale = "percent", main = "")+ 
  NoLegend()+ xlab(NULL)+ ylab(NULL)+
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))



# Cell frequency difference ------------------------------------------------------------------------------------

# SCpubr::do_BarPlot(sample = CD8_cancer_III_subset_Nonrecur, 
#                    group.by = "clusters", 
#                    split.by = "Patient_ID_SPC",
#                    position = "Fill",legend.position = "bottom",
#                    plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)
# 
# SCpubr::do_BarPlot(sample = CD8_cancer_III_subset_Recur, 
#                    group.by = "clusters", 
#                    split.by = "Patient_ID_SPC",
#                    position = "Fill",legend.position = "bottom",
#                    plot.title = "", xlab = "",  ylab = "Fraction of cells per cluster", font.size = 11)

dittoBarPlot(CD8_cancer_III_subset_Nonrecur, "seurat_clusters", group.by = "Patient_ID_SPC", scale = "percent", main = "")+ 
  NoLegend()+ xlab(NULL)+ ylab(NULL)+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))

dittoBarPlot(CD8_cancer_III_subset_Recur, "seurat_clusters", group.by = "Patient_ID_SPC", scale = "percent", main = "")+ 
  NoLegend()+ xlab(NULL)+ ylab(NULL)+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(face = "italic"), panel.border = element_rect(color = "black", fill = NA, size = 0.2))




dittoBarPlot(CD8_cancer_III_subset_Nonrecur, "seurat_clusters", group.by = "Patient_ID_SPC", scale = "percent", main = "")

##################################### Stats #########################################
##################################### Dataframe for Outcome
md.CD8_cancer_III_subset<- CD8_cancer_III_subset@meta.data %>% as.data.table

length(which(md.CD8_cancer_III_subset$Patient_ID==4 & md.CD8_cancer_III_subset$Type=='Tumor'& md.CD8_cancer_III_subset$seurat_clusters==1))

# Vector of sample numbers
sample_numbers <- c(4, 7, 8, 16, 20, 22, 23, 27, 28, 29, 36, 41, 50, 53, 54, 62, 65, 66)

# Vector of cluster names
cluster_vec <- unique(md.CD8_cancer_III_subset$seurat_clusters)
type_vec <- c("Nonrecur", "Recur")

unique(CD8_cancer_III_subset@meta.data$Patient_ID_SPC)

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
      subset_length <- length(which(md.CD8_cancer_III_subset$Patient_ID_SPC == number & md.CD8_cancer_III_subset$Outcome_new == type & md.CD8_cancer_III_subset$seurat_clusters == cluster))
      # Calculate the subset fraction
      total_count <- length(which(md.CD8_cancer_III_subset$Patient_ID_SPC == number & md.CD8_cancer_III_subset$Outcome_new == type))
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
DEG_CD8_cancer_III_subset <- CD8_cancer_III_subset
DEG_CD8_cancer_III_subset[["Total"]]  <- "All"


# pseudo-bulk workflow -----------------
# Acquiring necessary metrics for aggregation across cells in a sample
# 1. counts matrix - sample level
# counts aggregate to sample level

View(DEG_CD8_cancer_III_subset@meta.data)
DEG_CD8_cancer_III_subset$samples <- paste0(DEG_CD8_cancer_III_subset$Outcome_new, DEG_CD8_cancer_III_subset$Patient_ID_SPC)

DEG_CD8_cancer_III_subset
DefaultAssay(object = DEG_CD8_cancer_III_subset) <- "RNA"

cts <- AggregateExpression(DEG_CD8_cancer_III_subset, 
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

# SCpubr::do_VolcanoPlot(sample = DEG_CD8_cancer_III_subset, de_genes = condition.diffgenes, n_genes = 15)
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

VlnPlot(DEG_CD8_cancer_III_subset, features = "", split.by = "Outcome_new", pt.size = 0)
DefaultAssay(object = DEG_CD8_cancer_III_subset) <- "integrated"



# 오른쪽이 Non-recur

# perform DESeq2 -------- C0
# pseudo-bulk workflow -----------------
# Acquiring necessary metrics for aggregation across cells in a sample
# 1. counts matrix - sample level
# counts aggregate to sample level

View(DEG_CD8_cancer_III_subset@meta.data)
DEG_CD8_cancer_III_subset$samples <- paste0(DEG_CD8_cancer_III_subset$Outcome_new, DEG_CD8_cancer_III_subset$Patient_ID_SPC)

DEG_CD8_cancer_III_subset
DefaultAssay(object = DEG_CD8_cancer_III_subset) <- "RNA"

cts <- AggregateExpression(DEG_CD8_cancer_III_subset, 
                           group.by = c("seurat_clusters", "samples"),
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

# Let's run DE analysis with C0
# 1. Get counts matrix
counts_C0 <- cts.split.modified$`C0`
counts_C0

# 2. generate sample level metadata
colData <- data.frame(samples = colnames(counts_C0))

colData <- colData %>%
  mutate(condition = ifelse(grepl('Nonrecur', samples), 'Nonrecur', 'Recur')) %>%
  column_to_rownames(var = 'samples')


# perform DESeq2 --------
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_C0,
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


# perform DESeq2 -------- C1---------------------------------------------------------------------------------
# 1. Get counts matrix
counts_C1 <- cts.split.modified$`C1`
counts_C1

# 2. generate sample level metadata
colData <- data.frame(samples = colnames(counts_C1))

colData <- colData %>%
  mutate(condition = ifelse(grepl('Nonrecur', samples), 'Nonrecur', 'Recur')) %>%
  column_to_rownames(var = 'samples')


# perform DESeq2 --------
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_C1,
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


# perform DESeq2 -------- C2---------------------------------------------------------------------------------
# 1. Get counts matrix
counts_C2 <- cts.split.modified$`C2`
counts_C2

# 2. generate sample level metadata
colData <- data.frame(samples = colnames(counts_C2))

colData <- colData %>%
  mutate(condition = ifelse(grepl('Nonrecur', samples), 'Nonrecur', 'Recur')) %>%
  column_to_rownames(var = 'samples')


# perform DESeq2 --------
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_C2,
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



# perform DESeq2 -------- C3---------------------------------------------------------------------------------
# 1. Get counts matrix
counts_C3 <- cts.split.modified$`C3`
counts_C3

# 2. generate sample level metadata
colData <- data.frame(samples = colnames(counts_C3))

colData <- colData %>%
  mutate(condition = ifelse(grepl('Nonrecur', samples), 'Nonrecur', 'Recur')) %>%
  column_to_rownames(var = 'samples')


# perform DESeq2 --------
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_C3,
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


# perform DESeq2 -------- C4---------------------------------------------------------------------------------
# 1. Get counts matrix
counts_C4 <- cts.split.modified$`C4`
counts_C4

# 2. generate sample level metadata
colData <- data.frame(samples = colnames(counts_C4))

colData <- colData %>%
  mutate(condition = ifelse(grepl('Nonrecur', samples), 'Nonrecur', 'Recur')) %>%
  column_to_rownames(var = 'samples')


# perform DESeq2 --------
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_C4,
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
