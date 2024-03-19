library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
load("SeuratObjects/2024-03-11-00-49combined_seurat-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden_all_joined-tsne-anotated-alltsne.Rdata")
library(SingleCellExperiment)
load("./sceobject/2024-03-13-15-44miloobj-all_joined-tsne-anotated-alltsne.Rdata")
load("~/work/23-12-9_52WFru_singlecell_GSE218300/sceobject/2024-03-13-14-52-SCE-combined_seurat-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden_all_joined-tsne-anotated-alltsne.Rdata")

combined_seurat
table(Idents(combined_seurat))
