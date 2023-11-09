# load libraries
library(dplyr)
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
library(AUCell)
library(plotly)
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
library(SCpubr)
library(gt)
library(data.table)
library(magrittr)

install.packages("gt")

# Get data
immune_clusters <- readRDS(file="~/Single_cell_sample_merged/integrated_seurat_filtered_new.rds")

# DefaultAssay(immune_clusters) <- "integrated"
immune_clusters <- RunUMAP(immune_clusters, reduction = "pca", dims = 1:50) 
immune_clusters <- FindNeighbors(immune_clusters, reduction = "pca", dims = 1:50) 
immune_clusters <- FindClusters(immune_clusters, resolution = 5)

DimPlot(immune_clusters, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes()
DimPlot(immune_clusters, reduction = 'umap', raster=FALSE, label = FALSE) + NoLegend() +NoAxes()

immune_clusters_annotation <- immune_clusters

################################################################################# Annotation
#head(FindMarkers(immune_clusters, ident.1 = "27", ident.2 = "4", min.pct = 0.5, logfc.threshold = log(2)), n=20)

# B cells
#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "CD19", legend.position="none", min.cutoff =1, max.cutoff =3, pt.size = 2)
#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "MS4A1", legend.position="none", min.cutoff =1, max.cutoff =4, pt.size = 2)

immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '66' = "B cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '54' = "B cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '38' = "B cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `17` = "B cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `75` = "B cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `82` = "B cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `25` = "B cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `0` = "B cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `12` = "B cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `6` = "B cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `1` = "B cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `77` = "B cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `73` = "Proliferating B cells")

#DimPlot(immune_clusters_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes()


# Plasma cell
#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "SDC1", legend.position="none", min.cutoff =1, max.cutoff =3, pt.size = 2)
#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "JCHAIN", legend.position="none", min.cutoff =1, max.cutoff =4, pt.size = 2)

immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '70' = "Plasma cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '15' = "Plasma cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '36' = "Plasma cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `48` = "Plasma cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `83` = "Plasma cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `65` = "Plasma cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `68` = "Plasma cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `62` = "Plasma cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `79` = "Plasma cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `57` = "Plasma cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `78` = "Plasma cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, `24` = "Plasma cells")

#DimPlot(immune_clusters_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes()


# Mast cell
#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "KIT", legend.position="none", min.cutoff =1, max.cutoff =3, pt.size = 2)
#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "IL1RL1", legend.position="none", min.cutoff =1, max.cutoff =4, pt.size = 2)

immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '32' = "Mast cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '50' = "Mast cells")

#DimPlot(immune_clusters_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes()


# Myeloid cell
#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "ITGB2", legend.position="none", min.cutoff =0, max.cutoff =4, pt.size = 2)
#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "CD14", legend.position="none", min.cutoff =1, max.cutoff =4, pt.size = 2)
#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "FCN1", legend.position="none", min.cutoff =1, max.cutoff =4, pt.size = 2)
#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "CSF1R", legend.position="none", min.cutoff =1, max.cutoff =4, pt.size = 2)
#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "CD86", legend.position="none", min.cutoff =1, max.cutoff =10, pt.size = 2)
#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "CD1C", legend.position="none", min.cutoff =1, max.cutoff =15, pt.size = 2)

immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '9' = "Monocytes")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '26' = "Monocytes")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '37' = "Monocytes")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '22' = "Monocytes")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '55' = "Monocytes")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '81' = "Monocytes")

immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '45' = "Dendritic cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '43' = "Dendritic cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '80' = "Dendritic cells")

DimPlot(immune_clusters_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes()


# T/NK cell

#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "CD3G", legend.position="none", min.cutoff =0, max.cutoff =3, pt.size = 2)
#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "CD3E", legend.position="none", min.cutoff =0, max.cutoff =3, pt.size = 2)

#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "GNLY", legend.position="none", min.cutoff =0, max.cutoff =15, pt.size = 2)
#SCpubr::do_FeaturePlot(sample = immune_clusters, features = "NCAM1", legend.position="none", min.cutoff =6, max.cutoff =8, pt.size = 2)

immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '31' = "NK cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '42' = "NK cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '23' = "NK cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '29' = "NK cells")

#DimPlot(immune_clusters_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes()


#SCpubr::do_FeaturePlot(sample = immune_clusters_annotation, 
#                       idents.highlight = levels(immune_clusters_annotation)[!(levels(immune_clusters_annotation) 
#                                                                               %in% c('Dendritic cells','Monocytes','Plasma cells','B cells','Mast cells','Proliferating B cells'))], 
#                       features = c("MKI67"),min.cutoff =1, max.cutoff =10, pt.size = 2, legend.position="none")

#SCpubr::do_FeaturePlot(sample = immune_clusters_annotation, 
#                       idents.highlight = levels(immune_clusters_annotation)[!(levels(immune_clusters_annotation) 
#                                                                               %in% c('Dendritic cells','Monocytes','Plasma cells','B cells','Mast cells','Proliferating B cells'))], 
#                       features = c("HSPB1"),min.cutoff =1, max.cutoff =4, pt.size = 2, legend.position="none")

immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '64' = "Proliferating T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '59' = "Stressed / dying T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '41' = "Stressed / dying T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '71' = "Stressed / dying T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '76' = "Stressed / dying T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '40' = "Stressed / dying T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '30' = "Stressed / dying T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '69' = "Stressed / dying T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '52' = "Stressed / dying T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '53' = "Stressed / dying T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '58' = "Stressed / dying T cells")

#DimPlot(immune_clusters_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes()


# T cells specific

#SCpubr::do_FeaturePlot(sample = immune_clusters_annotation, 
#                       idents.highlight = levels(immune_clusters_annotation)[!(levels(immune_clusters_annotation) 
#                                                                               %in% c('Dendritic cells','Monocytes','Plasma cells','B cells','Mast cells','Proliferating B cells'))], 
#                       features = c("CD4"),min.cutoff =1, max.cutoff =3, pt.size = 2, legend.position="none")

#SCpubr::do_FeaturePlot(sample = immune_clusters_annotation, 
#                       idents.highlight = levels(immune_clusters_annotation)[!(levels(immune_clusters_annotation) 
#                                                                               %in% c('Dendritic cells','Monocytes','Plasma cells','B cells','Mast cells','Proliferating B cells'))], 
#                       features = c("CD8A"),min.cutoff =1, max.cutoff =4, pt.size = 2, legend.position="none")

#SCpubr::do_FeaturePlot(sample = immune_clusters_annotation, 
#                       idents.highlight = levels(immune_clusters_annotation)[!(levels(immune_clusters_annotation) 
#                                                                               %in% c('Dendritic cells','Monocytes','Plasma cells','B cells','Mast cells','Proliferating B cells'))], 
#                       features = c("FOXP3"),min.cutoff =1, max.cutoff =4, pt.size = 2, legend.position="none")


immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '4' = "Regulatory CD4 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '33' = "Regulatory CD4 T cells")
#DimPlot(immune_clusters_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes()

immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '49' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '47' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '72' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '46' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '3' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '7' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '8' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '61' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '16' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '21' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '19' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '14' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '28' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '51' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '67' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '10' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '63' = "CD8 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '13' = "CD8 T cells")

#DimPlot(immune_clusters_annotation, reduction = 'umap', raster=FALSE, label = TRUE) + NoLegend()+NoAxes()


immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '2' = "CD4 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '56' = "CD4 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '39' = "CD4 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '60' = "CD4 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '27' = "CD4 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '20' = "CD4 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '35' = "CD4 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '5' = "CD4 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '18' = "CD4 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '74' = "CD4 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '44' = "CD4 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '84' = "CD4 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '34' = "CD4 T cells")
immune_clusters_annotation <- RenameIdents(immune_clusters_annotation, '11' = "CD4 T cells")

DimPlot(immune_clusters_annotation, reduction = 'umap', raster=FALSE, label = TRUE, label.size = 6) + NoLegend()+NoAxes() +RotatedAxis()

immune_clusters_annotation[["immune_cluster"]]  <- immune_clusters_annotation@active.ident


################################################################################# Dot plot

genes <- c("CD3G","CD3E","CD4","CD8A","CCR7","FOXP3","GNLY","NCAM1","NKG7","HSPB1","HSPH1", "GADD45B","MKI67",
           "CD86","CD1C","ITGB2","CD14","FCN1","CSF1R","KIT","IL1RL1","SDC1","JCHAIN","CD19","MS4A1")

SCpubr::do_ExpressionHeatmap(sample = immune_clusters_annotation,
                             features = genes,
                             flip = F,
                             cluster_cols = FALSE,
                             cluster_rows = TRUE,
                             enforce_symmetry = TRUE,
                             use_viridis = FALSE,
                             max.cutoff = 6,
                             min.cutoff = 0
)


SCpubr::do_DotPlot(sample = immune_clusters_annotation, 
                   features = genes)

DefaultAssay(immune_clusters_annotation) <-"ADT"
Proteins <- c('CD45','CD8','CD4','CD25','CD56','CD11c','CD14','CD16','CD19')

SCpubr::do_DotPlot(sample = immune_clusters_annotation, 
                   features = Proteins)

VlnPlot(immune_clusters_annotation, "CD4", raster=FALSE, pt.size = 0) + NoLegend()
VlnPlot(immune_clusters_annotation, "CD8", raster=FALSE, pt.size = 0) + NoLegend()
VlnPlot(immune_clusters_annotation, "CD25", raster=FALSE, pt.size = 0, y.max = 2) + NoLegend()
VlnPlot(immune_clusters_annotation, "CD56", raster=FALSE, pt.size = 0, y.max = 2) + NoLegend()

VlnPlot(immune_clusters_annotation, "CD14", raster=FALSE, pt.size = 0, y.max = 3) + NoLegend()
VlnPlot(immune_clusters_annotation, "CD19", raster=FALSE, pt.size = 0, y.max = 5) + NoLegend()

SCpubr::do_ViolinPlot(sample = immune_clusters_annotation,features = "CD4", plot_boxplot = FALSE)
SCpubr::do_ViolinPlot(sample = immune_clusters_annotation,features = "CD8", plot_boxplot = FALSE)
SCpubr::do_ViolinPlot(sample = immune_clusters_annotation,features = "CD25", plot_boxplot = FALSE, y_cut = rep(2, length(immune_clusters_annotation)))


DefaultAssay(immune_clusters_annotation) <-"integrated"



saveRDS(immune_clusters_annotation, file="~/Single_cell_sample_merged/immune_clusters_annotation.rds")
saveRDS(immune_clusters_annotation@reductions$umap@cell.embeddings, "~/Single_cell_sample_merged/immune_clusters_annotation_umap_data_save_path.rds")
saveRDS(immune_clusters_annotation@meta.data, "~/Single_cell_sample_merged/immune_clusters_annotation_meta_data.rds")

################################################################################# Tumor vs Normal 

SCpubr::do_BarPlot(sample = immune_clusters_annotation, 
                   group.by = "immune_cluster", 
                   legend.position = "none", 
                   plot.title = "Number of cells per cluster")


SCpubr::do_BarPlot(immune_clusters_annotation,
                   group.by = 'immune_cluster',
                   split.by = "Type",
                   plot.title = "Frequency of immune clusters",
                   position = "fill",
                   legend.position = "none",
                   flip = FALSE)



de_genes <- tibble::tibble(Seurat::FindAllMarkers(object = immune_clusters_annotation))



SCpubr::do_GroupwiseDEPlot(sample = immune_clusters_annotation,
                           de_genes = de_genes,
                           min.cutoff = 1)


SCpubr::do_GroupwiseDEPlot(sample = immune_clusters_annotation,
                           de_genes = de_genes,
                           top_genes = 10)

SCpubr::do_GroupwiseDEPlot(sample = sample,
                           de_genes = de_genes,
                           column_title = "Title A",
                           row_title_p_values = "Title B",
                           row_title_logfc = "Title C",
                           row_title_expression = "Title D",
                           row_title_rot = 0)
