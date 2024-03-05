setwd("/home/rstudio/work/23-12-9_52WFru_singlecell_GSE218300")
.libPaths("/home/rstudio/R/x86_64-pc-linux-gnu-library/4.3")
.libPaths()
library(Seurat)
library(SeuratWrappers)
options(future.globals.maxSize = 1e9)
library(sceasy)
library(reticulate)
library(anndata)

#导出SCT的metadata整合到未sct的数据中---------------------------------------
load("~/work/23-12-9_52WFru_singlecell_GSE218300/SeuratObjects/combined_seurat_combined_umap_N4-N3-anatoedd-cluster-sct-harmony-ladden-diff-resolution-all_0_1-2-after-leiden-v2-2024-03-03.RData")

str(combined_seurat@meta.data)
metadata.sct <- combined_seurat@meta.data
head(metadata.sct)
str(metadata.sct)
save(metadata.sct, file = "SeuratObjects/metadata.sct-2024-03-05.RData")

load(metadata.sct, file = "SeuratObjects/metadata.sct-2024-03-05.RData")

load("~/work/23-12-9_52WFru_singlecell_GSE218300/SeuratObjects/combined_seurat_combined_umap_N4_harmony_cca_rcpa-mnn2024-03-01.RData")

combined_seurat@meta.data <- metadata.sct

save(combined_seurat, file = "SeuratObjects/combined_seurat_combined_umap_N4_harmony_cca_rcpa-mnn-withoutsct-addmeta-2024-03-05.RData")


#scVIIntegration-V1,用未SCT的数据----------------------------------------
load("SeuratObjects/combined_seurat_combined_umap_N4_harmony_cca_rcpa-mnn-withoutsct-addmeta-2024-03-05.RData")
#用SCT的数据
combined_seurat@active.assay

DefaultAssay(combined_seurat)
#也许可以用SCT的方法来分析，待定。

combined_seurat@active.assay <- "RNA"
use_condaenv("/home/rstudio/work/miniforge3/envs/scvi-env-4")
Sys.which("python")
combined_seurat <- IntegrateLayers(
  object = combined_seurat,
  method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "/home/rstudio/work/miniforge3/envs/scvi-env-4",
  verbose = FALSE
)

colnames(combined_seurat@meta.data)

save.dir <-"SeuratObjects"
save(
  combined_seurat,
  file = paste0(
    save.dir,
    "/combined_seurat_combined_umap_N4_harmony_cca_rcpa-mnn-withoutsct-addmeta-SCVI-intergraed",
    Sys.Date(),
    ".RData"
  )
)



#进入root登陆，解除限制。对umap不行
#use_condaenv("/home/rstudio/work/root_conda/anaconda3/envs/leiden")
use_condaenv("/home/rstudio/work/miniforge3/envs/scvi-env-4")
save.dir <-"SeuratObjects"
load(file = paste0(
    save.dir,
    "/combined_seurat_combined_umap_N4_harmony_cca_rcpa-mnn-withoutsct-addmeta-SCVI-intergraed",
    Sys.Date(),
    ".RData"
  ))


combined_seurat <-
  FindNeighbors(combined_seurat, reduction = "integrated.scvi", dims = 1:30)

combined_seurat <-
  RunUMAP(
    combined_seurat,
    reduction = "integrated.scvi",
    dims = 1:30,
    reduction.name = "umap.scvi"
  )

save(combined_seurat,file = paste0(
    save.dir,
    "/combined_seurat_combined_umap_N4_harmony_cca_rcpa-mnn-withoutsct-addmeta-SCVI-intergraed-umap",
    Sys.Date(),
    ".RData"
  ))

#umap-learn有bug
combined_seurat <-
  RunUMAP(
    combined_seurat,
    reduction = "integrated.scvi",
    dims = 1:30,
    umap.method="umap-learn",
    metric = "cosine",
    reduction.name = "umap-learn.scvi"
  )


load(file = paste0(
    save.dir,
    "/combined_seurat_combined_umap_N4_harmony_cca_rcpa-mnn-withoutsct-addmeta-SCVI-intergraed-umap",
    Sys.Date(),
    ".RData"
  ))
Reselute_clusters_leiden <- function(combined_seurat, resolution) {
  combined_seurat <-
    FindClusters(combined_seurat,
                 resolution = resolution,
                 # cluster.name = "sct_harmony_2_clusters",
                 algorithm = 4)
  return(combined_seurat)
}

Reselute_clusters_slm <- function(combined_seurat, resolution) {
  combined_seurat <-
    FindClusters(combined_seurat,
                 resolution = resolution,
                 # cluster.name = "sct_harmony_2_clusters",
                 algorithm = 3)
  return(combined_seurat)
}


Reselute_clusters_louvain_muti <-
  function(combined_seurat, resolution) {
    combined_seurat <-
      FindClusters(combined_seurat,
                   resolution = resolution,
                   # cluster.name = "sct_harmony_2_clusters",
                   algorithm = 2)
    return(combined_seurat)
  }

Reselute_clusters_louvain <- function(combined_seurat, resolution) {
  combined_seurat <-
    FindClusters(combined_seurat,
                 resolution = resolution,
                 # cluster.name = "sct_harmony_2_clusters",
                 algorithm = 1)
  return(combined_seurat)
}

resolution <-seq(0.1, 2, by = 0.1)

combined_seurat <-
  Reselute_clusters_slm(combined_seurat, resolution)

colnames(combined_seurat@meta.data)

# 通过gusb函数替换含有"SCT_snn_res"的列名
colnames(combined_seurat@meta.data) <-
  gsub(
    "RNA_snn_res.",
    "RNA_snn_scvi_slm_res.",
    colnames(combined_seurat@meta.data)
  )
colnames(combined_seurat@meta.data)
save(combined_seurat,file = paste0(
    save.dir,
    "/combined_seurat_combined_umap_N4_harmony_cca_rcpa-mnn-withoutsct-addmeta-SCVI-intergraed-umap-N1-slm",
    Sys.Date(),
    ".RData"
  ))


combined_seurat <-
  Reselute_clusters_louvain(combined_seurat, resolution)

colnames(combined_seurat@meta.data)

# 通过gusb函数替换含有"SCT_snn_res"的列名
colnames(combined_seurat@meta.data) <-
  gsub(
    "RNA_snn_res.",
    "RNA_snn_scvi_louvain_res.",
    colnames(combined_seurat@meta.data)
  )
colnames(combined_seurat@meta.data)
save(combined_seurat,file = paste0(
    save.dir,
    "/combined_seurat_combined_umap_N4_harmony_cca_rcpa-mnn-withoutsct-addmeta-SCVI-intergraed-umap-N2-slm-louvain",
    Sys.Date(),
    ".RData"
  ))



combined_seurat <-
  Reselute_clusters_louvain_muti(combined_seurat, resolution)

# 通过gusb函数替换含有"SCT_snn_res"的列名
colnames(combined_seurat@meta.data) <-
  gsub(
    "RNA_snn_res.",
    "RNA_snn_scvi_miti_louvain_res.",
    colnames(combined_seurat@meta.data)
  )
colnames(combined_seurat@meta.data)
save(combined_seurat,file = paste0(
    save.dir,
    "/combined_seurat_combined_umap_N4_harmony_cca_rcpa-mnn-withoutsct-addmeta-SCVI-intergraed-umap-N2-slm-louvain",
    Sys.Date(),
    ".RData"
  ))

tree_cluster <-
  function(combined_seurat,
           prefix,
           width = 25,
           height = 16) {
    library(clustree)
    # install.packages("clustree")
    p_cluster <-
      clustree(combined_seurat, prefix = prefix)

    # # install.packages("ploty")
    # install.packages("plotly")
    #
    # library(plotly)
    #
    # ggplotly(p_cluster)

    # width = 25
    # height = 16
    ggsave(
      width = width,
      height = height,
      paste0(
        prefix,
        "width_",
        width,
        "_",
        "height",
        height,
        "_",
        Sys.Date(),
        ".pdf"
      )
    )
    return(p_cluster)
  }

colnames(combined_seurat@meta.data)

