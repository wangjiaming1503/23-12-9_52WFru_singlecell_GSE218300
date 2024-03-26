library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(maartenutils)
library(SeuratDisk)
gen_file_overview("SeuratObjects")
View(maartenutils::gen_file_overview("SeuratObjects"))
load("SeuratObjects/2024-03-26-15-21combined_seurat_filter_combinedN5_scvi_tsne.RData")
library(SingleCellExperiment)

combined_seurat_filter
Layers(combined_seurat_filter)
combined_seurat_filter[["RNA"]]  <- JoinLayers(combined_seurat_filter[["RNA"]])

sce<- as.SingleCellExperiment(combined_seurat_filter)
