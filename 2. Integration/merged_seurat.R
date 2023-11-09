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

# Merge data ----
# Load filtered singlet Seurat object 
GC0516_1.CITE_filter_singlet <- readRDS(file="~/Single_cell_sample_filter_singlet/GC0516_1.CITE_filter_singlet.rds")
GC0516_2.CITE_filter_singlet <- readRDS(file="~/Single_cell_sample_filter_singlet/GC0516_2.CITE_filter_singlet.rds")
GC0516_3.CITE_filter_singlet <- readRDS(file="~/Single_cell_sample_filter_singlet/GC0516_3.CITE_filter_singlet.rds")

GC0517_1.CITE_filter_singlet <- readRDS(file="~/Single_cell_sample_filter_singlet/GC0517_1.CITE_filter_singlet.rds")
GC0517_2.CITE_filter_singlet <- readRDS(file="~/Single_cell_sample_filter_singlet/GC0517_2.CITE_filter_singlet.rds")
GC0517_3.CITE_filter_singlet <- readRDS(file="~/Single_cell_sample_filter_singlet/GC0517_3.CITE_filter_singlet.rds")

GC0530_1.CITE_filter_singlet <- readRDS(file="~/Single_cell_sample_filter_singlet/GC0530_1.CITE_filter_singlet.rds")
GC0530_2.CITE_filter_singlet <- readRDS(file="~/Single_cell_sample_filter_singlet/GC0530_2.CITE_filter_singlet.rds")
GC0530_3.CITE_filter_singlet <- readRDS(file="~/Single_cell_sample_filter_singlet/GC0530_3.CITE_filter_singlet.rds")

GC0531_1.CITE_filter_singlet <- readRDS(file="~/Single_cell_sample_filter_singlet/GC0531_1.CITE_filter_singlet.rds")
GC0531_2.CITE_filter_singlet <- readRDS(file="~/Single_cell_sample_filter_singlet/GC0531_2.CITE_filter_singlet.rds")
GC0531_3.CITE_filter_singlet <- readRDS(file="~/Single_cell_sample_filter_singlet/GC0531_3.CITE_filter_singlet.rds")

GC0711_1_1.CITE_filter_singlet <- readRDS(file="~/Single_cell_sample_filter_singlet/GC0711_1_1.CITE_filter_singlet.rds")
GC0711_1_2.CITE_filter_singlet <- readRDS(file="~/Single_cell_sample_filter_singlet/GC0711_1_2.CITE_filter_singlet.rds")
GC0711_1_3.CITE_filter_singlet <- readRDS(file="~/Single_cell_sample_filter_singlet/GC0711_1_3.CITE_filter_singlet.rds")

GC0711_2_1.CITE_filter_singlet <- readRDS(file="~/Single_cell_sample_filter_singlet/GC0711_2_1.CITE_filter_singlet.rds")
GC0711_2_2.CITE_filter_singlet <- readRDS(file="~/Single_cell_sample_filter_singlet/GC0711_2_2.CITE_filter_singlet.rds")

GC0712_1_1.CITE_filter_singlet <- readRDS(file="~/Single_cell_sample_filter_singlet/GC0712_1_1.CITE_filter_singlet.rds")
GC0712_1_2.CITE_filter_singlet <- readRDS(file="~/Single_cell_sample_filter_singlet/GC0712_1_2.CITE_filter_singlet.rds")

GC0712_2.CITE_filter_singlet <- readRDS(file="~/Single_cell_sample_filter_singlet/GC0712_2.CITE_filter_singlet.rds")

GC0718_1_1.CITE_filter_singlet <- readRDS(file="~/Single_cell_sample_filter_singlet/GC0718_1_1.CITE_filter_singlet.rds")
GC0718_1_2.CITE_filter_singlet <- readRDS(file="~/Single_cell_sample_filter_singlet/GC0718_1_2.CITE_filter_singlet.rds")

GC0718_2.CITE_filter_singlet <- readRDS(file="~/Single_cell_sample_filter_singlet/GC0718_2.CITE_filter_singlet.rds")

# Merge Seurat object
merged_seurat_singlet <- merge(GC0516_1.CITE_filter_singlet, y = c(GC0516_2.CITE_filter_singlet, GC0516_3.CITE_filter_singlet, GC0517_1.CITE_filter_singlet, 
                                                              GC0517_2.CITE_filter_singlet, GC0517_3.CITE_filter_singlet, GC0530_1.CITE_filter_singlet, 
                                                              GC0530_2.CITE_filter_singlet, GC0530_3.CITE_filter_singlet, GC0531_1.CITE_filter_singlet, 
                                                              GC0531_2.CITE_filter_singlet, GC0531_3.CITE_filter_singlet, GC0711_1_1.CITE_filter_singlet, 
                                                              GC0711_1_2.CITE_filter_singlet, GC0711_1_3.CITE_filter_singlet, GC0711_2_1.CITE_filter_singlet, 
                                                              GC0711_2_2.CITE_filter_singlet, GC0712_1_1.CITE_filter_singlet, GC0712_1_2.CITE_filter_singlet, GC0712_2.CITE_filter_singlet,
                                                              GC0718_1_1.CITE_filter_singlet, GC0718_1_2.CITE_filter_singlet, GC0718_2.CITE_filter_singlet
                                                              ),
                               project = "GC")

saveRDS(merged_seurat_singlet, file="~/Single_cell_sample_merged/merged_seurat_singlet.rds")