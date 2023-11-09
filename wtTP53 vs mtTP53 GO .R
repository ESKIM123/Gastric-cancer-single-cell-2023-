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
library(ggrepel)
library(enrichplot)


# Get data ----
immune_annotation <- readRDS(file="D:/DATA/Gastric cancer/immune_cell_annotation_final_OLD.rds")
DimPlot(immune_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Initial")

##################################### wtp53 vs mtp53 #########################################

view(cancer_MSS_subset@meta.data)

cancer_subset_MSIMSS <- subset(x = immune_annotation, Type == 'Tumor' & Patient_ID_SPC != FALSE)
cancer_subset_MSIMSS <- subset(x = cancer_subset_MSIMSS, Stage != 'Stage IV')
cancer_MSI_subset <- subset(x = cancer_subset_MSIMSS, MSI_type == 'MSI')
cancer_MSS_subset <- subset(x = cancer_subset_MSIMSS, MSI_type == 'MSS')

cancer_MSS_subset[["p53_type"]]  <- cancer_MSS_subset[["Molecular_type"]]

cancer_MSS_subset@meta.data$p53_type[cancer_MSS_subset@meta.data$p53_type == "GS"] <- "wtTP53"
cancer_MSS_subset@meta.data$p53_type[cancer_MSS_subset@meta.data$p53_type == "CIN"] <- "mtTP53"


##################################### wtp53 vs mtp53 #########################################

cancer_MSS_subset_DEG <- cancer_MSS_subset
cancer_MSS_subset_DEG <- SetIdent(cancer_MSS_subset_DEG, value = cancer_MSS_subset_DEG@meta.data$p53_type)


TP53_DEG <- FindMarkers(cancer_MSS_subset_DEG, ident.1 = "wtTP53", ident.2 ="mtTP53", method = "DESeq2")


TP53_DEG <-rownames_to_column(TP53_DEG, var = "Gene")
TP53_DEG$LOGPvalue <- -log10(TP53_DEG$p_val_adj)
TP53_DEG$Difference <- TP53_DEG$avg_log2FC

colnames(TP53_DEG)

SCpubr::do_VolcanoPlot(sample = cancer_MSS_subset_DEG, de_genes = TP53_DEG, n_genes = 15)


#VlnPlot(cancer_MSS_subset_DEG, features = "CTSB", split.by = "p53_type", pt.size = 0)
# 왼쪽이 mtTP53, 오른쪽이 wtTP53



ggplot(TP53_DEG, aes(x = Difference, y = LOGPvalue)) + 
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
    data = subset(TP53_DEG, LOGPvalue > 1.3 & abs(Difference) > 1),
    aes(label = Gene),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15))



# Selected genes


ggplot(TP53_DEG, aes(x = Difference, y = LOGPvalue)) + 
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
    data = subset(TP53_DEG, Gene %in% c("CTSB", "S100A6", "LGALS1")),
    aes(label = Gene),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15))
  

# Set threshold for p-value and fold change
pvalue_threshold <- 0.05
fold_change_threshold <- 0.25

# Filter dataframe for significant genes
TP53_DEG_significant_genes <- TP53_DEG[TP53_DEG$LOGPvalue > -log10(pvalue_threshold) & abs(TP53_DEG$Difference) > fold_change_threshold, ]
TP53_DEG_significant_genes_up <- TP53_DEG[TP53_DEG$LOGPvalue > -log10(pvalue_threshold) & TP53_DEG$Difference > fold_change_threshold, ]
TP53_DEG_significant_genes_down <- TP53_DEG[TP53_DEG$LOGPvalue > -log10(pvalue_threshold) & TP53_DEG$Difference < -fold_change_threshold, ]

sum(TP53_DEG_significant_genes$Difference > 0.25)
sum(TP53_DEG_significant_genes$Difference < -0.25)




################################################################################### GO ####
# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db" 

keytypes(org.Hs.eg.db)

BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

# reading in data 
df <- TP53_DEG_significant_genes
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

view(gse_BP_GO@result)
write.csv(gse@result, file = "D:/DATA/TP53_total_go_result.csv", row.names = FALSE)

require(DOSE)


gse_BP_GO <- gse_BP
gse_BP_GO@result$Check <- NA
ids_selected <- c("GO:0000165", "GO:0048041", "GO:1901889")


if()

gse_BP_GO@result$ID <- subset(gse_BP_GO@result$ID, gse_BP_GO@result$ID == ids_selected)


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



















# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)

ridgeplot(gse) + labs(x = "enrichment distribution")

# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(gse, by = "all", title = gse$Description[100], geneSetID = 1)

view(gse@result)
gene_list



