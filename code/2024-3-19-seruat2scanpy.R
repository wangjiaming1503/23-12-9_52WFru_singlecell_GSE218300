# failed seurat disk to anndata---------------------------
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
# library(SeuratDisk)
setwd("/home/rstudio/work/23-12-9_52WFru_singlecell_GSE218300")
View(maartenutils::gen_file_overview("SeuratObjects"))
load("./SeuratObjects/2024-03-11-00-49combined_seurat-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden_all_joined-tsne-anotated-alltsne.Rdata")
get_time()
DefaultAssay(combined_seurat) <- "RNA"
combined_seurat
DefaultAssay(combined_seurat)
Layers(combined_seurat)


save.dir <- "SeuratObjects"
SaveH5Seurat(combined_seurat, filename = paste0(
    save.dir,
    "/", get_time(), "anotated-alltsne-combined_seurat.h5Seurat"
)) # only able to save in h5Seurat format,not h5ad


# Convert("SeuratObjects/2024-03-19-09-07anotated-alltsne-combined_seurat.h5Seurat", dest = "h5ad") #faied
# convert seurat to anndata---------------------------------------------------
install.packages("scCustomize")
install.packages("rliger")
# Load Packages
library(Seurat)
library(rliger)
library(scCustomize)
library(qs)
get_time()
# as.anndata(x = pbmc, file_path = "~/Desktop", file_name = "pbmc_anndata.h5ad")
as.anndata(combined_seurat, file_path = "SeuratObjects", file_name = "2024-03-19-09-35anotated-alltsne-combined_seurat.h5ad")

DefaultAssay(combined_seurat) <- "SCT"
get_time()
as.anndata(combined_seurat, file_path = "SeuratObjects", file_name = "2024-03-19-09-39anotated-alltsne-combined_seurat-SCT.h5ad") # same result
