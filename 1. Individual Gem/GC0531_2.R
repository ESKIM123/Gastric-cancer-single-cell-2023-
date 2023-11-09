# Setup the Seurat Object
library(dplyr)
library(plyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(DoubletFinder)
library(Matrix)
library(ggeasy)
library(tidyr)
library(glmGamPoi)

# ggplot title centered
theme_update(plot.title = element_text(hjust = 0.5))
theme_update(plot.title = element_text(face = "bold"))


## Load GC0531_3 ----
################################################################################# Load the dataset 
GC0531_3.data <- Read10X(data.dir = "~/DATA_GC/TBD220852_15884_20220910 - 5' GEX/02_cellranger_file/0531-3/outs/filtered_feature_bc_matrix/")
GC0531_3.data.FB <- Read10X(data.dir = "~/DATA_GC/TBD220852_15883_20220906 - 5' FB/02_cellranger_file/0531-3/outs/raw_feature_bc_matrix/")
GC0531_3.data.TCR <- read.csv("~/DATA_GC/TBD220852_15812_20220829 - TCR/02_cellranger_file/0531-3/outs/all_contig_annotations.csv", header = T)

# Dividing FB data into Hashtag and ADT(CITE) 
GC0531_3.data.FB_Hashtag <- GC0531_3.data.FB[grep("Hashtag",rownames(GC0531_3.data.FB)),]; 
GC0531_3.data.FB_CITE <- GC0531_3.data.FB[-grep("Hashtag",rownames(GC0531_3.data.FB)),] 

## HTO
GC0531_3.joint.bcs.HTO <- intersect(colnames(GC0531_3.data), colnames(GC0531_3.data.FB_Hashtag)) 
GC0531_3.data.HTO <- GC0531_3.data[, GC0531_3.joint.bcs.HTO] 
GC0531_3.data.FB_Hashtag <- as.matrix(GC0531_3.data.FB_Hashtag[, GC0531_3.joint.bcs.HTO]) 

# Setup Seurat object
GC0531_3.hashtag <- CreateSeuratObject(counts = GC0531_3.data.HTO) 

GC0531_3.hashtag <- NormalizeData(GC0531_3.hashtag) 
GC0531_3.hashtag <- FindVariableFeatures(GC0531_3.hashtag, selection.method = "mean.var.plot") 
GC0531_3.hashtag <- ScaleData(GC0531_3.hashtag, features = VariableFeatures(GC0531_3.hashtag)) 

# Add HTO data as a new assay independent from RNA
GC0531_3.hashtag[["HTO"]] <- CreateAssayObject(counts = GC0531_3.data.FB_Hashtag) 

GC0531_3.hashtag <- NormalizeData(GC0531_3.hashtag, assay = "HTO", normalization.method = "CLR")
GC0531_3.hashtag <- HTODemux(GC0531_3.hashtag, assay = "HTO", positive.quantile = 0.99)
Idents(GC0531_3.hashtag) <- "HTO_maxID"
Idents(GC0531_3.hashtag) <- "HTO_classification.global"

GC0531_3.hashtag.subset <- subset(GC0531_3.hashtag, idents = "Negative", invert = TRUE)

DefaultAssay(GC0531_3.hashtag.subset) <- "HTO"
GC0531_3.hashtag.subset <- ScaleData(GC0531_3.hashtag.subset, features = rownames(GC0531_3.hashtag.subset), verbose = FALSE)
GC0531_3.hashtag.subset <- RunPCA(GC0531_3.hashtag.subset, features = rownames(GC0531_3.hashtag.subset), approx = FALSE)
GC0531_3.hashtag.subset <- RunTSNE(GC0531_3.hashtag.subset, dims = 1:8, perplexity = 100)
GC0531_3.singlet <- subset(GC0531_3.hashtag, idents = "Singlet")
GC0531_3.singlet <- FindVariableFeatures(GC0531_3.singlet, selection.method = "mean.var.plot")
GC0531_3.singlet <- ScaleData(GC0531_3.singlet, features = VariableFeatures(GC0531_3.singlet))
GC0531_3.singlet <- RunPCA(GC0531_3.singlet, features = VariableFeatures(GC0531_3.singlet))
GC0531_3.singlet <- FindNeighbors(GC0531_3.singlet, reduction = "pca", dims = 1:10)
GC0531_3.singlet <- FindClusters(GC0531_3.singlet, resolution = 0.6, verbose = FALSE)
GC0531_3.singlet <- RunTSNE(GC0531_3.singlet, reduction = "pca", dims = 1:10)


# Global classification results
table(GC0531_3.hashtag$HTO_classification.global)
VlnPlot(GC0531_3.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE) + ggtitle("nCount_RNA (Gem:GC0531_3)")
HTOHeatmap(GC0531_3.hashtag, assay = "HTO", ncells = 10000) + ggtitle("HTO_Heatmap (Gem:GC0531_3)") 
# Projecting singlet identities on TSNE visualization
DimPlot(GC0531_3.singlet, group.by = "HTO_classification") + ggtitle("HTO_classification (Gem:GC0531_3)") 

## HTO 完----

## ADT
GC0531_3.joint.bcs.ADT <- intersect(colnames(GC0531_3.data), colnames(GC0531_3.data.FB_CITE))
GC0531_3.data.ADT <- GC0531_3.data[, GC0531_3.joint.bcs.ADT]
GC0531_3.data.FB_CITE <- as.matrix(GC0531_3.data.FB_CITE[, GC0531_3.joint.bcs.ADT])

# Setup Seurat object
GC0531_3.CITE <- CreateSeuratObject(counts = GC0531_3.data.ADT)
GC0531_3.adt_assay <- CreateAssayObject(counts = GC0531_3.data.FB_CITE)

# Add this assay to the previously created Seurat object
GC0531_3.CITE[["ADT"]] <- GC0531_3.adt_assay

# Perform visualization and clustering steps
GC0531_3.CITE <- NormalizeData(GC0531_3.CITE)
GC0531_3.CITE <- FindVariableFeatures(GC0531_3.CITE)
GC0531_3.CITE <- ScaleData(GC0531_3.CITE)
GC0531_3.CITE <- RunPCA(GC0531_3.CITE, verbose = FALSE)
GC0531_3.CITE <- FindNeighbors(GC0531_3.CITE, dims = 1:30)
GC0531_3.CITE <- FindClusters(GC0531_3.CITE, resolution = 0.8, verbose = FALSE)
GC0531_3.CITE <- RunUMAP(GC0531_3.CITE, dims = 1:30)

# Normalize ADT data
GC0531_3.CITE <- NormalizeData(GC0531_3.CITE, normalization.method = "CLR", margin = 2, assay = "ADT")
DimPlot(GC0531_3.CITE, label = TRUE) + ggtitle("DimPlot_CITE (Gem:GC0531_3)") + ggeasy::easy_center_title()

## ADT 完----

## HTO ADT Integration 

GC0531_3.CITE@meta.data$hash.id<-
  mapvalues(rownames(GC0531_3.CITE@meta.data), 
            from=rownames(GC0531_3.singlet@meta.data), 
            to=as.character(GC0531_3.singlet$hash.ID))

rownames(GC0531_3.CITE@meta.data%>%filter(!hash.id %in% unique(as.character(GC0531_3.singlet$hash.ID))))

hash_na_cell<-rownames(GC0531_3.CITE@meta.data%>%filter(!hash.id %in% unique(as.character(GC0531_3.singlet$hash.ID))))
GC0531_3.CITE@meta.data[hash_na_cell,"hash.id"]<-NA

################################################################################# Hashtag name
GC0531_3.CITE$Patient_ID<-mapvalues(GC0531_3.CITE$hash.id, 
                                    from=c("Hashtag1", "Hashtag2", "Hashtag3", "Hashtag4", "Hashtag5", "Hashtag6" ),
                                    to=c("15", "20", "22", "40", "41", "50" ))

GC0531_3.CITE$seurat_clusters <- NULL # Metadata column 지우기
GC0531_3.CITE$RNA_snn_res.0.8 <- NULL

## HTO ADT Integration 完----

## TCR 

GC0531_3.data.TCR_true <- subset(GC0531_3.data.TCR, productive == "true")
GC0531_3_BAR <- unique(GC0531_3.data.TCR_true$barcode)
GC0531_3.data.TCR_true_BAR <- data.frame(matrix(NA, nrow = length(GC0531_3_BAR), ncol = 3))
colnames(GC0531_3.data.TCR_true_BAR) <- c("barcode", "alpha", "beta")
GC0531_3.data.TCR_true_BAR$barcode <- GC0531_3_BAR

for(i in 1:length(GC0531_3_BAR)){
  GC0531_3.data.TCR_true.temp.alpha <- GC0531_3.data.TCR_true %>% subset(barcode == GC0531_3_BAR[i] & chain == "TRA")
  if(nrow(GC0531_3.data.TCR_true.temp.alpha) == 1){
    GC0531_3.data.TCR_true_BAR$alpha[i] <- GC0531_3.data.TCR_true.temp.alpha$cdr3_nt
  }
  else if(nrow(GC0531_3.data.TCR_true.temp.alpha) == 0){
    GC0531_3.data.TCR_true_BAR$alpha[i] <- NA
  }
  else {
    GC0531_3.data.TCR_true.temp.alpha <- GC0531_3.data.TCR_true.temp.alpha %>% subset(umis == max(GC0531_3.data.TCR_true.temp.alpha$umis))
    GC0531_3.data.TCR_true_BAR$alpha[i] <- GC0531_3.data.TCR_true.temp.alpha$cdr3_nt
  }
  
  GC0531_3.data.TCR_true.temp.beta <- GC0531_3.data.TCR_true %>% subset(barcode == GC0531_3_BAR[i] & chain == "TRB")
  if(nrow(GC0531_3.data.TCR_true.temp.beta) == 1){
    GC0531_3.data.TCR_true_BAR$beta[i] <- GC0531_3.data.TCR_true.temp.beta$cdr3_nt
  }
  else if(nrow(GC0531_3.data.TCR_true.temp.beta) == 0){
    GC0531_3.data.TCR_true_BAR$beta[i] <- NA
  }
  else {
    GC0531_3.data.TCR_true.temp.beta <- GC0531_3.data.TCR_true.temp.beta %>% subset(umis == max(GC0531_3.data.TCR_true.temp.beta$umis))
    GC0531_3.data.TCR_true_BAR$beta[i] <- GC0531_3.data.TCR_true.temp.beta$cdr3_nt
  }
}

GC0531_3.data.TCR_true_BAR <- GC0531_3.data.TCR_true_BAR %>% subset(is.na(alpha) == F & is.na(beta) == F)
GC0531_3.data.TCR_true_BAR <- data.frame(GC0531_3.data.TCR_true_BAR, row.names = 1)

nrow(GC0531_3.data.TCR_true)
nrow(GC0531_3.data.TCR_true_BAR)
nrow(GC0531_3.CITE@meta.data)


## TCR 完----

## TCR Integration 

GC0531_3.CITE@meta.data$alpha<-
  mapvalues(rownames(GC0531_3.CITE@meta.data), 
            from=rownames(GC0531_3.data.TCR_true_BAR), 
            to=as.character(GC0531_3.data.TCR_true_BAR$alpha))

alpha_cell<-rownames(GC0531_3.CITE@meta.data%>%filter(!alpha %in% unique(as.character(GC0531_3.data.TCR_true_BAR$alpha))))
GC0531_3.CITE@meta.data[alpha_cell,"alpha"]<-NA


GC0531_3.CITE@meta.data$beta<-
  mapvalues(rownames(GC0531_3.CITE@meta.data), 
            from=rownames(GC0531_3.data.TCR_true_BAR), 
            to=as.character(GC0531_3.data.TCR_true_BAR$beta))

beta_cell<-rownames(GC0531_3.CITE@meta.data%>%filter(!beta %in% unique(as.character(GC0531_3.data.TCR_true_BAR$beta))))
GC0531_3.CITE@meta.data[alpha_cell,"beta"]<-NA


## TCR Integration 完----

## DoubletFinder
################################################################################# Create Seurat object

GC0531_3.CITE[["Gem"]] <- 'GC0531_3'
GC0531_3.CITE[["Type"]] <- 'Tumor'
GC0531_3.CITE[["Outcome"]] <- 'Non_recur'

# Quality control 
# % MT reads

GC0531_3.CITE [["percent.mt"]] <- PercentageFeatureSet(GC0531_3.CITE, pattern = "^MT-")
Idents(object = GC0531_3.CITE)=0

((VlnPlot(GC0531_3.CITE, features = c("nFeature_RNA"), pt.size = 0, raster=FALSE, y.max = 6000) + theme(legend.position = 'none')) |
  (VlnPlot(GC0531_3.CITE, features = c("nCount_RNA"), pt.size = 0, raster=FALSE, y.max = 15000) + theme(legend.position = 'none')) |
  (VlnPlot(GC0531_3.CITE, features = c("percent.mt"), pt.size = 0, raster=FALSE, y.max = 15) + theme(legend.position = 'none')))

((FeatureScatter(GC0531_3.CITE, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
  + geom_smooth(method = 'lm'))+ theme(legend.position = 'none')) | 
  (FeatureScatter(GC0531_3.CITE, feature1 = "nCount_RNA", feature2 = "percent.mt")+ theme(legend.position = 'none'))

# Filtering 
GC0531_3.CITE <- subset(GC0531_3.CITE, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
GC0531_3.CITE <- NormalizeData(GC0531_3.CITE, normalization.method = "LogNormalize", scale.factor = 10000)
GC0531_3.CITE <- FindVariableFeatures(GC0531_3.CITE, selection.method = "vst", nfeatures = 2000)
LabelPoints(plot = VariableFeaturePlot(GC0531_3.CITE), points = head(VariableFeatures(GC0531_3.CITE), 10), repel = TRUE)

# Scaling
GC0531_3.CITE <- ScaleData(GC0531_3.CITE, features = rownames(GC0531_3.CITE))

# Linear dimensionality reduction (PCA)
GC0531_3.CITE <- RunPCA(GC0531_3.CITE, features = VariableFeatures(object = GC0531_3.CITE))
DimHeatmap(GC0531_3.CITE, dims = 1, cells = 500, balanced = TRUE)

pct_GC0531_3.CITE <- GC0531_3.CITE[["pca"]]@stdev / sum(GC0531_3.CITE[["pca"]]@stdev) * 100
cumu_GC0531_3.CITE <- cumsum(pct_GC0531_3.CITE)
co1_GC0531_3.CITE <- which(cumu_GC0531_3.CITE > 90 & pct_GC0531_3.CITE < 5)[1]
co2_GC0531_3.CITE <- sort(which((pct_GC0531_3.CITE[1:length(pct_GC0531_3.CITE) - 1] - pct_GC0531_3.CITE[2:length(pct_GC0531_3.CITE)]) > 0.1), decreasing = T)[1] + 1
pcs_GC0531_3.CITE <- min(co1_GC0531_3.CITE, co2_GC0531_3.CITE)

# Create a dataframe with values
plot_df_GC0531_3.CITE <- data.frame(pct_GC0531_3.CITE = pct_GC0531_3.CITE, 
                               cumu_GC0531_3.CITE = cumu_GC0531_3.CITE, 
                               rank = 1:length(pct_GC0531_3.CITE))
# Elbow plot to visualize 
ElbowPlot(GC0531_3.CITE, ndims = 50) |
  ggplot(plot_df_GC0531_3.CITE, aes(cumu_GC0531_3.CITE, pct_GC0531_3.CITE, label = rank, color = rank > pcs_GC0531_3.CITE)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct_GC0531_3.CITE[pct_GC0531_3.CITE > 5]), color = "grey") +
  theme_bw()

# Optimal dimensionality
co3_GC0531_3.CITE <- co2_GC0531_3.CITE +1

# Clustering & UMAP
GC0531_3.CITE <- FindNeighbors(object = GC0531_3.CITE, dims = 1:co3_GC0531_3.CITE)
GC0531_3.CITE <- FindClusters(object = GC0531_3.CITE)
GC0531_3.CITE <- RunUMAP(object = GC0531_3.CITE, dims = 1:co3_GC0531_3.CITE)
DimPlot(GC0531_3.CITE, reduction = "umap")

# DoubletFinder 
## pK Identification (no ground-truth)
sweep.res.list_GC0531_3.CITE <- paramSweep_v3(GC0531_3.CITE, PCs = 1:co3_GC0531_3.CITE, sct = FALSE)
sweep.stats_GC0531_3.CITE <- summarizeSweep(sweep.res.list_GC0531_3.CITE, GT = FALSE)
bcmvn_GC0531_3.CITE <- find.pK(sweep.stats_GC0531_3.CITE)

ggplot(bcmvn_GC0531_3.CITE, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line() + ggtitle("pK value (Gem:GC0531_3)")

pK <- bcmvn_GC0531_3.CITE %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate
annotations_GC0531_3.CITE <- GC0531_3.CITE@meta.data$seurat_clusters
homotypic.prop_GC0531_3.CITE <- modelHomotypic(annotations_GC0531_3.CITE)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi_GC0531_3.CITE <- round(0.054*nrow(GC0531_3.CITE@meta.data))  ## Assuming 5.4 % doublet formation rate 
nExp_poi.adj_GC0531_3.CITE <- round(nExp_poi_GC0531_3.CITE*(1-homotypic.prop_GC0531_3.CITE))

# run doubletFinder 
GC0531_3.CITE <- doubletFinder_v3(GC0531_3.CITE, 
                             PCs = 1:co3_GC0531_3.CITE, 
                             pN = 0.25, 
                             pK = pK, 
                             nExp = nExp_poi.adj_GC0531_3.CITE,
                             reuse.pANN = FALSE, sct = FALSE)
names(GC0531_3.CITE@meta.data)
pK
DimPlot(GC0531_3.CITE, reduction = 'umap', group.by = colnames(GC0531_3.CITE@meta.data)[grep("DF.classifications", colnames(GC0531_3.CITE@meta.data))], raster=FALSE)

################################################################################# Visualize doublets 
GC0531_3.CITE$DF.classifications <- GC0531_3.CITE$DF.classifications_0.25_0.27_192
GC0531_3.CITE$DF.classifications_0.25_0.27_192 <- NULL
GC0531_3.CITE$pANN_0.25_0.27_192 <- NULL

# number of singlets and doublets
table(GC0531_3.CITE@meta.data$DF.classifications)
GC0531_3.CITE <- subset(GC0531_3.CITE, subset = DF.classifications == 'Singlet')
table(GC0531_3.CITE@meta.data$DF.classifications)

## DoubletFinder 完----


## Using sctransform in Seurat

# run sctransform

GC0531_3.CITE <- SCTransform(GC0531_3.CITE, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)

# These are now standard steps in the Seurat workflow for visualization and clustering
GC0531_3.CITE <- RunPCA(GC0531_3.CITE, verbose = FALSE)
GC0531_3.CITE <- RunUMAP(GC0531_3.CITE, dims = 1:30, verbose = FALSE)

GC0531_3.CITE <- FindNeighbors(GC0531_3.CITE, dims = 1:30, verbose = FALSE)
GC0531_3.CITE <- FindClusters(GC0531_3.CITE, verbose = FALSE)
DimPlot(GC0531_3.CITE, label = TRUE) + NoLegend() + ggtitle("SCTransform (Gem:GC0531_3)") 

## SCTransform 完----

saveRDS(GC0531_3.CITE, file="~/Single_cell_sample_filter_singlet/GC0531_3.CITE_filter_singlet.rds")

## Save 完----
