setwd("/home/rstudio/work/23-12-9_52WFru_singlecell_GSE218300")
.libPaths("/home/rstudio/R/x86_64-pc-linux-gnu-library/4.3")
.libPaths()
library(Seurat)
library(SeuratWrappers)
options(future.globals.maxSize = 1e9)
load("~/work/23-12-9_52WFru_singlecell_GSE218300/SeuratObjects/combined_seurat_combined_umap_N4-N3-anatoedd-cluster-sct-harmony-ladden-diff-resolution-all_0_1-2-after-leiden-v2-2024-03-03.RData")
library(sceasy)
library(reticulate)
library(anndata)

use_condaenv("/home/rstudio/work/anaconda3/envs/scvi-env-2")
Sys.which("python")
combined_seurat <- IntegrateLayers(
  object = combined_seurat,
  method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "../anaconda3/envs/scvi-env-2",
  verbose = FALSE
)

