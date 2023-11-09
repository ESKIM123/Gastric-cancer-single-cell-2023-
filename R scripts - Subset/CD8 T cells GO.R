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
devtools::install_github("immunogenomics/presto")
library(presto)


# Get data ----
CD8_subset_anno <- readRDS(file="D:/DATA/Gastric cancer/GS_subsets/CD8_subset.rds")
DimPlot(CD8_subset_anno, reduction = 'umap', raster=FALSE, label = TRUE) + ggtitle("Initial")

TCR.subset <- CD8_subset_anno

TCR.subset[["TCR_type"]]  <- TCR.subset[["Patient_ID_SPC"]]


TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 22] <- "C0_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 59] <- "C0_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 8] <- "C0_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 66] <- "C0_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 28] <- "C0_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 18] <- "C0_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 65] <- "C0_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 16] <- "C0_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 44] <- "C0_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 53] <- "C0_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 56] <- "C0_EXP"

TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 73] <- "C1_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 15] <- "C1_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 71] <- "C1_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 54] <- "C1_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 23] <- "C1_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 60] <- "C1_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 20] <- "C1_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 36] <- "C1_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 40] <- "C1_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 29] <- "C1_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 68] <- "C1_EXP"

TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 52] <- "C2_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 35] <- "C2_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 50] <- "C2_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 7] <- "C2_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 27] <- "C2_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 75] <- "C2_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 4] <- "C2_EXP"
TCR.subset@meta.data["TCR_type"][TCR.subset@meta.data["TCR_type"] == 41] <- "C2_EXP"

CD8_TCRtpye <- TCR.subset

cancer_subset <- subset(x = CD8_TCRtpye, Type == 'Tumor' & Patient_ID_SPC != FALSE & Stage != 'Stage IV')
##################################### TCR type comparison #########################################



##################################### C0_EXP #########################################

cancer_subset_GO <- cancer_subset
cancer_subset_GO <- SetIdent(cancer_subset_GO, value = cancer_subset_GO@meta.data$TCR_type)

###### Findmarkers
C0_EXP_DEG <- FindMarkers(cancer_subset_GO, ident.1 = "C0_EXP", ident.2 = c("C1_EXP", "C2_EXP"), method = "DESeq2")

###### 
C0_EXP_DEG <-rownames_to_column(C0_EXP_DEG, var = "Gene")
C0_EXP_DEG$LOGPvalue <- -log10(C0_EXP_DEG$p_val_adj)
C0_EXP_DEG$Difference <- C0_EXP_DEG$avg_log2FC

colnames(C0_EXP_DEG)

SCpubr::do_VolcanoPlot(sample = cancer_subset_GO, de_genes = C0_EXP_DEG, n_genes = 15)


VlnPlot(cancer_subset_GO, features = "ITGA9", split.by = "TCR_type", pt.size = 0)
# 오른쪽이 C0_EXP_DEG



ggplot(C0_EXP_DEG, aes(x = Difference, y = LOGPvalue)) + 
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
    data = subset(C0_EXP_DEG, LOGPvalue > 1.3 & abs(Difference) > 1),
    aes(label = Gene),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15))



# Selected genes


ggplot(C0_EXP_DEG, aes(x = Difference, y = LOGPvalue)) + 
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
    data = subset(C0_EXP_DEG, Gene %in% c("CTSB", "S100A6", "LGALS1")),
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
C0_EXP_DEG_significant_genes <- C0_EXP_DEG[C0_EXP_DEG$LOGPvalue > -log10(pvalue_threshold) & abs(C0_EXP_DEG$Difference) > fold_change_threshold, ]
C0_EXP_DEG_significant_genes_up <- C0_EXP_DEG[C0_EXP_DEG$LOGPvalue > -log10(pvalue_threshold) & C0_EXP_DEG$Difference > fold_change_threshold, ]
C0_EXP_DEG_significant_genes_down <- C0_EXP_DEG[C0_EXP_DEG$LOGPvalue > -log10(pvalue_threshold) & C0_EXP_DEG$Difference < -fold_change_threshold, ]

sum(C0_EXP_DEG_significant_genes$Difference > 0.25)
sum(C0_EXP_DEG_significant_genes$Difference < -0.25)

write.csv(C0_EXP_DEG_significant_genes_up, file = "D:/DATA/C0_EXP_DEG_significant_genes_up.csv", row.names = FALSE)




##################################### C1_EXP #########################################

cancer_subset_GO <- cancer_subset
cancer_subset_GO <- SetIdent(cancer_subset_GO, value = cancer_subset_GO@meta.data$TCR_type)

###### Findmarkers
C1_EXP_DEG <- FindMarkers(cancer_subset_GO, ident.1 = "C1_EXP", ident.2 = c("C0_EXP", "C2_EXP"), method = "DESeq2")

###### 
C1_EXP_DEG <-rownames_to_column(C1_EXP_DEG, var = "Gene")
C1_EXP_DEG$LOGPvalue <- -log10(C1_EXP_DEG$p_val_adj)
C1_EXP_DEG$Difference <- C1_EXP_DEG$avg_log2FC

colnames(C1_EXP_DEG)

SCpubr::do_VolcanoPlot(sample = cancer_subset_GO, de_genes = C1_EXP_DEG, n_genes = 15)


VlnPlot(cancer_subset_GO, features = "ITGA9", split.by = "TCR_type", pt.size = 0)
# 오른쪽이 C1_EXP_DEG



ggplot(C1_EXP_DEG, aes(x = Difference, y = LOGPvalue)) + 
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
    data = subset(C1_EXP_DEG, LOGPvalue > 1.3 & abs(Difference) > 1),
    aes(label = Gene),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15))



# Selected genes


ggplot(C1_EXP_DEG, aes(x = Difference, y = LOGPvalue)) + 
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
    data = subset(C1_EXP_DEG, Gene %in% c("CTSB", "S100A6", "LGALS1")),
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
C1_EXP_DEG_significant_genes <- C1_EXP_DEG[C1_EXP_DEG$LOGPvalue > -log10(pvalue_threshold) & abs(C1_EXP_DEG$Difference) > fold_change_threshold, ]
C1_EXP_DEG_significant_genes_up <- C1_EXP_DEG[C1_EXP_DEG$LOGPvalue > -log10(pvalue_threshold) & C1_EXP_DEG$Difference > fold_change_threshold, ]
C1_EXP_DEG_significant_genes_down <- C1_EXP_DEG[C1_EXP_DEG$LOGPvalue > -log10(pvalue_threshold) & C1_EXP_DEG$Difference < -fold_change_threshold, ]

sum(C1_EXP_DEG_significant_genes$Difference > 0.25)
sum(C1_EXP_DEG_significant_genes$Difference < -0.25)

write.csv(C1_EXP_DEG_significant_genes_up, file = "D:/DATA/C1_EXP_DEG_significant_genes_up.csv", row.names = FALSE)



##################################### C2_EXP #########################################

cancer_subset_GO <- cancer_subset
cancer_subset_GO <- SetIdent(cancer_subset_GO, value = cancer_subset_GO@meta.data$TCR_type)

###### Findmarkers
C2_EXP_DEG <- FindMarkers(cancer_subset_GO, ident.1 = "C2_EXP", ident.2 = c("C0_EXP", "C1_EXP"), method = "DESeq2")

###### 
C2_EXP_DEG <-rownames_to_column(C2_EXP_DEG, var = "Gene")
C2_EXP_DEG$LOGPvalue <- -log10(C2_EXP_DEG$p_val_adj)
C2_EXP_DEG$Difference <- C2_EXP_DEG$avg_log2FC

colnames(C2_EXP_DEG)

SCpubr::do_VolcanoPlot(sample = cancer_subset_GO, de_genes = C2_EXP_DEG, n_genes = 15)


VlnPlot(cancer_subset_GO, features = "ITGA9", split.by = "TCR_type", pt.size = 0)
# 오른쪽이 C2_EXP_DEG



ggplot(C2_EXP_DEG, aes(x = Difference, y = LOGPvalue)) + 
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
    data = subset(C2_EXP_DEG, LOGPvalue > 1.3 & abs(Difference) > 1),
    aes(label = Gene),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15))



# Selected genes


ggplot(C2_EXP_DEG, aes(x = Difference, y = LOGPvalue)) + 
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
    data = subset(C2_EXP_DEG, Gene %in% c("CTSB", "S100A6", "LGALS1")),
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
C2_EXP_DEG_significant_genes <- C2_EXP_DEG[C2_EXP_DEG$LOGPvalue > -log10(pvalue_threshold) & abs(C2_EXP_DEG$Difference) > fold_change_threshold, ]
C2_EXP_DEG_significant_genes_up <- C2_EXP_DEG[C2_EXP_DEG$LOGPvalue > -log10(pvalue_threshold) & C2_EXP_DEG$Difference > fold_change_threshold, ]
C2_EXP_DEG_significant_genes_down <- C2_EXP_DEG[C2_EXP_DEG$LOGPvalue > -log10(pvalue_threshold) & C2_EXP_DEG$Difference < -fold_change_threshold, ]

sum(C2_EXP_DEG_significant_genes$Difference > 0.25)
sum(C2_EXP_DEG_significant_genes$Difference < -0.25)

write.csv(C2_EXP_DEG_significant_genes_up, file = "D:/DATA/C2_EXP_DEG_significant_genes_up.csv", row.names = FALSE)



################################################################################### GO ####
# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db" 

keytypes(org.Hs.eg.db)

BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

# reading in data 
df <- C0_EXP_DEG_significant_genes
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


write.csv(gse@result, file = "D:/DATA/CO_EXP_total_go_result.csv", row.names = FALSE)

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



######################## Draw figure (dotplot)
library(msigdbr)
msigdbr_show_species()
print(msigdbr_collections(), n=30)
m_df<- msigdbr(species = "Homo sapiens", category = "H") 
m_df<- msigdbr(species = "Homo sapiens", category = "H") 
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)



#불러온 Gene Set을 살펴봅니다.
head(m_df)
dplyr::count(myData.genes, group)
#그런 다음에 제가 원하는 그룹을 뽑아서 logFC와 AUC score를 빼낸 뒤 정리해봅니다. 저는 CD14+ Mono_cALD 그룹에서 어떠한 pathway들이 더 혹은 덜 발현되는지 알아보고 싶습니다.

myData.genes %>%
  dplyr::filter(group == "CD14+ Mono_Control") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)

#그리고 뽑은 정보들을 auc에 따라 1등부터 꼴등까지 줄을 세워봅니다.
CD14Mono.genes<- myData.genes %>%
  dplyr::filter(group == "CD14+ Mono_Control") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

ranks<- deframe(CD14Mono.genes)
head(ranks)
이제 실제로 GSEA를 실행합니다. 저희는 이번에 fgsea라는 패키지를 사용할 겁니다. 바이오컨덕터로 간편하게 설치 가능합니다.

C2_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, gene_symbol)
head(C2_t2g)

em2 <- GSEA(geneList=gene_list, TERM2GENE = C2_t2g)
head(em2)

view(em2@result)

# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(em2, by = "all", title = em2$Description[1], geneSetID = 1)


################################################################################### Load data

convac = read.csv("D:/Myocarditis/GSEA/convac_significant_genes.csv", header=TRUE)
conmyo = read.csv("D:/Myocarditis/GSEA/conmyo_significant_genes.csv", header=TRUE)
myovac = read.csv("D:/Myocarditis/GSEA/myovac_significant_genes.csv", header=TRUE)

################################################################################### Extract gene lists

convac_original_gene_list <- convac$log2FoldChange
names(convac_original_gene_list) <- convac$X
convac_gene_list <- na.omit(convac_original_gene_list)
convac_gene_list = sort(convac_gene_list, decreasing = TRUE)
convac_gene_list

conmyo_original_gene_list <- conmyo$log2FoldChange
names(conmyo_original_gene_list) <- conmyo$X
conmyo_gene_list <- na.omit(conmyo_original_gene_list)
conmyo_gene_list = sort(conmyo_gene_list, decreasing = TRUE)
conmyo_gene_list

myovac_original_gene_list <- myovac$log2FoldChange
names(myovac_original_gene_list) <- myovac$X
myovac_gene_list <- na.omit(myovac_original_gene_list)
myovac_gene_list = sort(myovac_gene_list, decreasing = TRUE)
myovac_gene_list

################################################################################### Set up GSEA 

library(msigdbr)
msigdbr_show_species()
print(msigdbr_collections(), n=25)

C5_msig <- msigdbr(species = "Homo sapiens", category = "C5") %>% dplyr::select(gs_name, gene_symbol)
All_msig <- msigdbr(species = "Homo sapiens") %>% dplyr::select(gs_name, gene_symbol)

Selected_gs <- All_msig[grep("MYOCARDIAL_FIBROSIS", All_msig$gs_name), ]


################################################################################### Run GSEA

convac_GSEA <- GSEA(geneList=convac_gene_list, TERM2GENE = Selected_gs, pvalueCutoff = 1)
conmyo_GSEA <- GSEA(geneList=conmyo_gene_list, TERM2GENE = Selected_gs)
myovac_GSEA <- GSEA(geneList=myovac_gene_list, TERM2GENE = All_msig)

view(convac_GSEA@result)
view(conmyo_GSEA@result)
view(myovac_GSEA@result)

myovac_GSEA <- gseGO(geneList=myovac_gene_list, 
                     ont ="ALL", 
                     keyType = "SYMBOL", 
                     nPerm = 10000, 
                     minGSSize = 3, 
                     maxGSSize = 800, 
                     pvalueCutoff = 0.05, 
                     verbose = TRUE, 
                     OrgDb = All_msig, 
                     pAdjustMethod = "none")



################################################################################### Plot GSEA

gseaplot(myovac_GSEA, by = "all", title = myovac_GSEA$Description[1], geneSetID = 1)
































gse_BP_GO@result$Check <- NA
ids_selected <- c("GO:0000165", "GO:0048041", "GO:1901889")


if()
  
  gse_BP_GO@result$ID <- subset(gse_BP_GO@result$ID, gse_BP_GO@result$ID == ids_selected)













# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)

ridgeplot(gse) + labs(x = "enrichment distribution")

# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(gse, by = "all", title = gse$Description[100], geneSetID = 1)

view(gse@result)
gene_list
































