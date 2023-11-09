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


## Load GC0711_1_1 ----
################################################################################# Load the dataset 
GC0711_1_1.data <- Read10X(data.dir = "~/DATA_GC/TBD220852_15884_20220910 - 5' GEX/02_cellranger_file/0711_1-1/outs/filtered_feature_bc_matrix/")
GC0711_1_1.data.FB <- Read10X(data.dir = "~/DATA_GC/TBD220852_15883_20220906 - 5' FB/02_cellranger_file/0711_1-1/outs/raw_feature_bc_matrix/")
GC0711_1_1.data.TCR <- read.csv("~/DATA_GC/TBD220852_15812_20220829 - TCR/02_cellranger_file/0711_1-1/outs/all_contig_annotations.csv", header = T)

# Dividing FB data into Hashtag and ADT(CITE) 
GC0711_1_1.data.FB_Hashtag <- GC0711_1_1.data.FB[grep("Hashtag",rownames(GC0711_1_1.data.FB)),]; 
GC0711_1_1.data.FB_CITE <- GC0711_1_1.data.FB[-grep("Hashtag",rownames(GC0711_1_1.data.FB)),] 

## HTO
GC0711_1_1.joint.bcs.HTO <- intersect(colnames(GC0711_1_1.data), colnames(GC0711_1_1.data.FB_Hashtag)) 
GC0711_1_1.data.HTO <- GC0711_1_1.data[, GC0711_1_1.joint.bcs.HTO] 
GC0711_1_1.data.FB_Hashtag <- as.matrix(GC0711_1_1.data.FB_Hashtag[, GC0711_1_1.joint.bcs.HTO]) 

# Setup Seurat object
GC0711_1_1.hashtag <- CreateSeuratObject(counts = GC0711_1_1.data.HTO) 

GC0711_1_1.hashtag <- NormalizeData(GC0711_1_1.hashtag) 
GC0711_1_1.hashtag <- FindVariableFeatures(GC0711_1_1.hashtag, selection.method = "mean.var.plot") 
GC0711_1_1.hashtag <- ScaleData(GC0711_1_1.hashtag, features = VariableFeatures(GC0711_1_1.hashtag)) 

# Add HTO data as a new assay independent from RNA
GC0711_1_1.hashtag[["HTO"]] <- CreateAssayObject(counts = GC0711_1_1.data.FB_Hashtag) 

GC0711_1_1.hashtag <- NormalizeData(GC0711_1_1.hashtag, assay = "HTO", normalization.method = "CLR")
GC0711_1_1.hashtag <- HTODemux(GC0711_1_1.hashtag, assay = "HTO", positive.quantile = 0.99)
Idents(GC0711_1_1.hashtag) <- "HTO_maxID"
Idents(GC0711_1_1.hashtag) <- "HTO_classification.global"

GC0711_1_1.hashtag.subset <- subset(GC0711_1_1.hashtag, idents = "Negative", invert = TRUE)

DefaultAssay(GC0711_1_1.hashtag.subset) <- "HTO"
GC0711_1_1.hashtag.subset <- ScaleData(GC0711_1_1.hashtag.subset, features = rownames(GC0711_1_1.hashtag.subset), verbose = FALSE)
GC0711_1_1.hashtag.subset <- RunPCA(GC0711_1_1.hashtag.subset, features = rownames(GC0711_1_1.hashtag.subset), approx = FALSE)
GC0711_1_1.hashtag.subset <- RunTSNE(GC0711_1_1.hashtag.subset, dims = 1:8, perplexity = 100)
GC0711_1_1.singlet <- subset(GC0711_1_1.hashtag, idents = "Singlet")
GC0711_1_1.singlet <- FindVariableFeatures(GC0711_1_1.singlet, selection.method = "mean.var.plot")
GC0711_1_1.singlet <- ScaleData(GC0711_1_1.singlet, features = VariableFeatures(GC0711_1_1.singlet))
GC0711_1_1.singlet <- RunPCA(GC0711_1_1.singlet, features = VariableFeatures(GC0711_1_1.singlet))
GC0711_1_1.singlet <- FindNeighbors(GC0711_1_1.singlet, reduction = "pca", dims = 1:10)
GC0711_1_1.singlet <- FindClusters(GC0711_1_1.singlet, resolution = 0.6, verbose = FALSE)
GC0711_1_1.singlet <- RunTSNE(GC0711_1_1.singlet, reduction = "pca", dims = 1:10)


# Global classification results
table(GC0711_1_1.hashtag$HTO_classification.global)
VlnPlot(GC0711_1_1.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE) + ggtitle("nCount_RNA (Gem:GC0711_1_1)")
HTOHeatmap(GC0711_1_1.hashtag, assay = "HTO", ncells = 10000) + ggtitle("HTO_Heatmap (Gem:GC0711_1_1)") 
# Projecting singlet identities on TSNE visualization
DimPlot(GC0711_1_1.singlet, group.by = "HTO_classification") + ggtitle("HTO_classification (Gem:GC0711_1_1)") 

## HTO 完----

## ADT
GC0711_1_1.joint.bcs.ADT <- intersect(colnames(GC0711_1_1.data), colnames(GC0711_1_1.data.FB_CITE))
GC0711_1_1.data.ADT <- GC0711_1_1.data[, GC0711_1_1.joint.bcs.ADT]
GC0711_1_1.data.FB_CITE <- as.matrix(GC0711_1_1.data.FB_CITE[, GC0711_1_1.joint.bcs.ADT])

# Setup Seurat object
GC0711_1_1.CITE <- CreateSeuratObject(counts = GC0711_1_1.data.ADT)
GC0711_1_1.adt_assay <- CreateAssayObject(counts = GC0711_1_1.data.FB_CITE)

# Add this assay to the previously created Seurat object
GC0711_1_1.CITE[["ADT"]] <- GC0711_1_1.adt_assay

# Perform visualization and clustering steps
GC0711_1_1.CITE <- NormalizeData(GC0711_1_1.CITE)
GC0711_1_1.CITE <- FindVariableFeatures(GC0711_1_1.CITE)
GC0711_1_1.CITE <- ScaleData(GC0711_1_1.CITE)
GC0711_1_1.CITE <- RunPCA(GC0711_1_1.CITE, verbose = FALSE)
GC0711_1_1.CITE <- FindNeighbors(GC0711_1_1.CITE, dims = 1:30)
GC0711_1_1.CITE <- FindClusters(GC0711_1_1.CITE, resolution = 0.8, verbose = FALSE)
GC0711_1_1.CITE <- RunUMAP(GC0711_1_1.CITE, dims = 1:30)

# Normalize ADT data
GC0711_1_1.CITE <- NormalizeData(GC0711_1_1.CITE, normalization.method = "CLR", margin = 2, assay = "ADT")
DimPlot(GC0711_1_1.CITE, label = TRUE) + ggtitle("DimPlot_CITE (Gem:GC0711_1_1)") + ggeasy::easy_center_title()

## ADT 完----

## HTO ADT Integration 

GC0711_1_1.CITE@meta.data$hash.id<-
  mapvalues(rownames(GC0711_1_1.CITE@meta.data), 
            from=rownames(GC0711_1_1.singlet@meta.data), 
            to=as.character(GC0711_1_1.singlet$hash.ID))

rownames(GC0711_1_1.CITE@meta.data%>%filter(!hash.id %in% unique(as.character(GC0711_1_1.singlet$hash.ID))))

hash_na_cell<-rownames(GC0711_1_1.CITE@meta.data%>%filter(!hash.id %in% unique(as.character(GC0711_1_1.singlet$hash.ID))))
GC0711_1_1.CITE@meta.data[hash_na_cell,"hash.id"]<-NA

################################################################################# Hashtag name
GC0711_1_1.CITE$Patient_ID<-mapvalues(GC0711_1_1.CITE$hash.id, 
                                    from=c("Hashtag1", "Hashtag2", "Hashtag3", "Hashtag4", "Hashtag5", "Hashtag6" ),
                                    to=c("52", "59", "60", "65", "66", "75" ))

GC0711_1_1.CITE$seurat_clusters <- NULL # Metadata column 지우기
GC0711_1_1.CITE$RNA_snn_res.0.8 <- NULL

## HTO ADT Integration 完----

## TCR 

GC0711_1_1.data.TCR_true <- subset(GC0711_1_1.data.TCR, productive == "true")
GC0711_1_1_BAR <- unique(GC0711_1_1.data.TCR_true$barcode)
GC0711_1_1.data.TCR_true_BAR <- data.frame(matrix(NA, nrow = length(GC0711_1_1_BAR), ncol = 3))
colnames(GC0711_1_1.data.TCR_true_BAR) <- c("barcode", "alpha", "beta")
GC0711_1_1.data.TCR_true_BAR$barcode <- GC0711_1_1_BAR

for(i in 1:length(GC0711_1_1_BAR)){
  GC0711_1_1.data.TCR_true.temp.alpha <- GC0711_1_1.data.TCR_true %>% subset(barcode == GC0711_1_1_BAR[i] & chain == "TRA")
  if(nrow(GC0711_1_1.data.TCR_true.temp.alpha) == 1){
    GC0711_1_1.data.TCR_true_BAR$alpha[i] <- GC0711_1_1.data.TCR_true.temp.alpha$cdr3_nt
  }
  else if(nrow(GC0711_1_1.data.TCR_true.temp.alpha) == 0){
    GC0711_1_1.data.TCR_true_BAR$alpha[i] <- NA
  }
  else {
    GC0711_1_1.data.TCR_true.temp.alpha <- GC0711_1_1.data.TCR_true.temp.alpha %>% subset(umis == max(GC0711_1_1.data.TCR_true.temp.alpha$umis))
    GC0711_1_1.data.TCR_true_BAR$alpha[i] <- GC0711_1_1.data.TCR_true.temp.alpha$cdr3_nt
  }
  
  GC0711_1_1.data.TCR_true.temp.beta <- GC0711_1_1.data.TCR_true %>% subset(barcode == GC0711_1_1_BAR[i] & chain == "TRB")
  if(nrow(GC0711_1_1.data.TCR_true.temp.beta) == 1){
    GC0711_1_1.data.TCR_true_BAR$beta[i] <- GC0711_1_1.data.TCR_true.temp.beta$cdr3_nt
  }
  else if(nrow(GC0711_1_1.data.TCR_true.temp.beta) == 0){
    GC0711_1_1.data.TCR_true_BAR$beta[i] <- NA
  }
  else {
    GC0711_1_1.data.TCR_true.temp.beta <- GC0711_1_1.data.TCR_true.temp.beta %>% subset(umis == max(GC0711_1_1.data.TCR_true.temp.beta$umis))
    GC0711_1_1.data.TCR_true_BAR$beta[i] <- GC0711_1_1.data.TCR_true.temp.beta$cdr3_nt
  }
}

GC0711_1_1.data.TCR_true_BAR <- GC0711_1_1.data.TCR_true_BAR %>% subset(is.na(alpha) == F & is.na(beta) == F)
GC0711_1_1.data.TCR_true_BAR <- data.frame(GC0711_1_1.data.TCR_true_BAR, row.names = 1)

nrow(GC0711_1_1.data.TCR_true)
nrow(GC0711_1_1.data.TCR_true_BAR)
nrow(GC0711_1_1.CITE@meta.data)


## TCR 完----

## TCR Integration 

GC0711_1_1.CITE@meta.data$alpha<-
  mapvalues(rownames(GC0711_1_1.CITE@meta.data), 
            from=rownames(GC0711_1_1.data.TCR_true_BAR), 
            to=as.character(GC0711_1_1.data.TCR_true_BAR$alpha))

alpha_cell<-rownames(GC0711_1_1.CITE@meta.data%>%filter(!alpha %in% unique(as.character(GC0711_1_1.data.TCR_true_BAR$alpha))))
GC0711_1_1.CITE@meta.data[alpha_cell,"alpha"]<-NA


GC0711_1_1.CITE@meta.data$beta<-
  mapvalues(rownames(GC0711_1_1.CITE@meta.data), 
            from=rownames(GC0711_1_1.data.TCR_true_BAR), 
            to=as.character(GC0711_1_1.data.TCR_true_BAR$beta))

beta_cell<-rownames(GC0711_1_1.CITE@meta.data%>%filter(!beta %in% unique(as.character(GC0711_1_1.data.TCR_true_BAR$beta))))
GC0711_1_1.CITE@meta.data[alpha_cell,"beta"]<-NA


## TCR Integration 完----

## DoubletFinder
################################################################################# Create Seurat object

GC0711_1_1.CITE[["Gem"]] <- 'GC0711_1_1'
GC0711_1_1.CITE[["Type"]] <- 'Tumor'
GC0711_1_1.CITE[["Outcome"]] <- 'Non_recur'

# Quality control 
# % MT reads

GC0711_1_1.CITE [["percent.mt"]] <- PercentageFeatureSet(GC0711_1_1.CITE, pattern = "^MT-")
Idents(object = GC0711_1_1.CITE)=0

((VlnPlot(GC0711_1_1.CITE, features = c("nFeature_RNA"), pt.size = 0, raster=FALSE, y.max = 6000) + theme(legend.position = 'none')) |
  (VlnPlot(GC0711_1_1.CITE, features = c("nCount_RNA"), pt.size = 0, raster=FALSE, y.max = 15000) + theme(legend.position = 'none')) |
  (VlnPlot(GC0711_1_1.CITE, features = c("percent.mt"), pt.size = 0, raster=FALSE, y.max = 15) + theme(legend.position = 'none')))

((FeatureScatter(GC0711_1_1.CITE, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
  + geom_smooth(method = 'lm'))+ theme(legend.position = 'none')) | 
  (FeatureScatter(GC0711_1_1.CITE, feature1 = "nCount_RNA", feature2 = "percent.mt")+ theme(legend.position = 'none'))

# Filtering 
GC0711_1_1.CITE <- subset(GC0711_1_1.CITE, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
GC0711_1_1.CITE <- NormalizeData(GC0711_1_1.CITE, normalization.method = "LogNormalize", scale.factor = 10000)
GC0711_1_1.CITE <- FindVariableFeatures(GC0711_1_1.CITE, selection.method = "vst", nfeatures = 2000)
LabelPoints(plot = VariableFeaturePlot(GC0711_1_1.CITE), points = head(VariableFeatures(GC0711_1_1.CITE), 10), repel = TRUE)

# Scaling
GC0711_1_1.CITE <- ScaleData(GC0711_1_1.CITE, features = rownames(GC0711_1_1.CITE))

# Linear dimensionality reduction (PCA)
GC0711_1_1.CITE <- RunPCA(GC0711_1_1.CITE, features = VariableFeatures(object = GC0711_1_1.CITE))
DimHeatmap(GC0711_1_1.CITE, dims = 1, cells = 500, balanced = TRUE)

pct_GC0711_1_1.CITE <- GC0711_1_1.CITE[["pca"]]@stdev / sum(GC0711_1_1.CITE[["pca"]]@stdev) * 100
cumu_GC0711_1_1.CITE <- cumsum(pct_GC0711_1_1.CITE)
co1_GC0711_1_1.CITE <- which(cumu_GC0711_1_1.CITE > 90 & pct_GC0711_1_1.CITE < 5)[1]
co2_GC0711_1_1.CITE <- sort(which((pct_GC0711_1_1.CITE[1:length(pct_GC0711_1_1.CITE) - 1] - pct_GC0711_1_1.CITE[2:length(pct_GC0711_1_1.CITE)]) > 0.1), decreasing = T)[1] + 1
pcs_GC0711_1_1.CITE <- min(co1_GC0711_1_1.CITE, co2_GC0711_1_1.CITE)

# Create a dataframe with values
plot_df_GC0711_1_1.CITE <- data.frame(pct_GC0711_1_1.CITE = pct_GC0711_1_1.CITE, 
                               cumu_GC0711_1_1.CITE = cumu_GC0711_1_1.CITE, 
                               rank = 1:length(pct_GC0711_1_1.CITE))
# Elbow plot to visualize 
ElbowPlot(GC0711_1_1.CITE, ndims = 50) |
  ggplot(plot_df_GC0711_1_1.CITE, aes(cumu_GC0711_1_1.CITE, pct_GC0711_1_1.CITE, label = rank, color = rank > pcs_GC0711_1_1.CITE)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct_GC0711_1_1.CITE[pct_GC0711_1_1.CITE > 5]), color = "grey") +
  theme_bw()

# Optimal dimensionality
co3_GC0711_1_1.CITE <- co2_GC0711_1_1.CITE +1

# Clustering & UMAP
GC0711_1_1.CITE <- FindNeighbors(object = GC0711_1_1.CITE, dims = 1:co3_GC0711_1_1.CITE)
GC0711_1_1.CITE <- FindClusters(object = GC0711_1_1.CITE)
GC0711_1_1.CITE <- RunUMAP(object = GC0711_1_1.CITE, dims = 1:co3_GC0711_1_1.CITE)
DimPlot(GC0711_1_1.CITE, reduction = "umap")

# DoubletFinder 
## pK Identification (no ground-truth)
sweep.res.list_GC0711_1_1.CITE <- paramSweep_v3(GC0711_1_1.CITE, PCs = 1:co3_GC0711_1_1.CITE, sct = FALSE)
sweep.stats_GC0711_1_1.CITE <- summarizeSweep(sweep.res.list_GC0711_1_1.CITE, GT = FALSE)
bcmvn_GC0711_1_1.CITE <- find.pK(sweep.stats_GC0711_1_1.CITE)

ggplot(bcmvn_GC0711_1_1.CITE, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line() + ggtitle("pK value (Gem:GC0711_1_1)")

pK <- bcmvn_GC0711_1_1.CITE %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate
annotations_GC0711_1_1.CITE <- GC0711_1_1.CITE@meta.data$seurat_clusters
homotypic.prop_GC0711_1_1.CITE <- modelHomotypic(annotations_GC0711_1_1.CITE)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi_GC0711_1_1.CITE <- round(0.054*nrow(GC0711_1_1.CITE@meta.data))  ## Assuming 5.4 % doublet formation rate 
nExp_poi.adj_GC0711_1_1.CITE <- round(nExp_poi_GC0711_1_1.CITE*(1-homotypic.prop_GC0711_1_1.CITE))

# run doubletFinder 
GC0711_1_1.CITE <- doubletFinder_v3(GC0711_1_1.CITE, 
                             PCs = 1:co3_GC0711_1_1.CITE, 
                             pN = 0.25, 
                             pK = pK, 
                             nExp = nExp_poi.adj_GC0711_1_1.CITE,
                             reuse.pANN = FALSE, sct = FALSE)
names(GC0711_1_1.CITE@meta.data)
pK
DimPlot(GC0711_1_1.CITE, reduction = 'umap', group.by = colnames(GC0711_1_1.CITE@meta.data)[grep("DF.classifications", colnames(GC0711_1_1.CITE@meta.data))], raster=FALSE)

################################################################################# Visualize doublets 
GC0711_1_1.CITE$DF.classifications <- GC0711_1_1.CITE$DF.classifications_0.25_0.24_265
GC0711_1_1.CITE$DF.classifications_0.25_0.24_265 <- NULL
GC0711_1_1.CITE$pANN_0.25_0.24_265 <- NULL

# number of singlets and doublets
table(GC0711_1_1.CITE@meta.data$DF.classifications)
GC0711_1_1.CITE <- subset(GC0711_1_1.CITE, subset = DF.classifications == 'Singlet')
table(GC0711_1_1.CITE@meta.data$DF.classifications)

## DoubletFinder 完----


## Using sctransform in Seurat

# run sctransform

GC0711_1_1.CITE <- SCTransform(GC0711_1_1.CITE, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)

# These are now standard steps in the Seurat workflow for visualization and clustering
GC0711_1_1.CITE <- RunPCA(GC0711_1_1.CITE, verbose = FALSE)
GC0711_1_1.CITE <- RunUMAP(GC0711_1_1.CITE, dims = 1:30, verbose = FALSE)

GC0711_1_1.CITE <- FindNeighbors(GC0711_1_1.CITE, dims = 1:30, verbose = FALSE)
GC0711_1_1.CITE <- FindClusters(GC0711_1_1.CITE, verbose = FALSE)
DimPlot(GC0711_1_1.CITE, label = TRUE) + NoLegend() + ggtitle("SCTransform (Gem:GC0711_1_1)") 

## SCTransform 完----

saveRDS(GC0711_1_1.CITE, file="~/Single_cell_sample_filter_singlet/GC0711_1_1.CITE_filter_singlet.rds")

## Save 完----
