setwd("/home/rstudio/work/23-12-9_52WFru_singlecell_GSE218300")

.libPaths()
library(Seurat)
library(SeuratWrappers)
options(future.globals.maxSize = 1e9)
library(sceasy)
library(reticulate)
library(anndata)
library(ggplot2)
get_time <- function(my.tz = "Asia/Shanghai") {
  # 获取当前时间，并按照指定时区和格式化字符串
  datetime <- format(Sys.time(), "%Y-%m-%d-%H-%M", tz = my.tz)

  # 返回格式化的日期时间字符串
  return(datetime)
}
get_time()

save.dir <- "SeuratObjects"


load("SeuratObjects/combined_seurat_combined_umap_N4_harmony_cca_rcpa-mnn-keepsct-SCVI-intergraed-umap-add-learn-N3-slm-louvain-muti_joined2024-03-06-00-27.RData")

combined_seurat

# 使用leiden
use_condaenv("/home/rstudio/work/root_conda/anaconda3/envs/leiden")

Reselute_clusters_leiden <- function(combined_seurat, resolution) {
  combined_seurat <-
    FindClusters(combined_seurat,
      resolution = resolution,
      algorithm = 4
    )
  return(combined_seurat)
}

resolution <- seq(0.1, 2, by = 0.1)
resolution[-10]

combined_seurat <-
  Reselute_clusters_leiden(combined_seurat, 1)

colnames(combined_seurat@meta.data)
colnames(combined_seurat@meta.data) <-
  gsub(
    "RNA_snn_res.",
    "RNA_snn_scvi_leiden_res.",
    colnames(combined_seurat@meta.data)
  )

save(combined_seurat, file = paste0(
  save.dir,
  "/combined_seurat_combined_umap_N4_harmony_cca_rcpa-mnn-keepsct-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden-1_joined",
  get_time(),
  ".RData"
))
library("loupeR")

create_loupe_from_seurat(
  combined_seurat,
  output_name = paste0(
    "SeuratObjects/combined_seurat_combined_umap_N4_harmony_cca_rcpa-mnn-keepsct-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden-1_joined",
    format(Sys.Date(), "%Y-%m-%d")
  )
)

combined_seurat <-
  Reselute_clusters_leiden(combined_seurat, resolution[-10])

colnames(combined_seurat@meta.data)
colnames(combined_seurat@meta.data) <-
  gsub(
    "RNA_snn_res.",
    "RNA_snn_scvi_leiden_res.",
    colnames(combined_seurat@meta.data)
  )

get_time()
save(combined_seurat, file = paste0(
  save.dir,
  "/combined_seurat_combined_umap_N4_harmony_cca_rcpa-mnn-keepsct-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden_all_joined",
  get_time(),
  ".RData"
))

create_loupe_from_seurat(
  combined_seurat,
  output_name = paste0(
    "SeuratObjects/combined_seurat_combined_umap_N4_harmony_cca_rcpa-mnn-keepsct-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden_all_joined",
    format(Sys.Date(), "%Y-%m-%d")
  )
)

DefaultAssay(combined_seurat)

DefaultAssay(combined_seurat) <- "SCT"

create_loupe_from_seurat(
  combined_seurat,
  output_name = paste0(
    "SeuratObjects/combined_seurat_combined_umap_N4_harmony_cca_rcpa-mnn-keepsct-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden_all_joined_SCTassey",
    get_time()
  )
)
library(clustree)
tree_cluster <-
  function(combined_seurat,
           prefix,
           width = 25,
           height = 20) {
    library(clustree)

    p_cluster <-
      clustree(combined_seurat, prefix = prefix)


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
        get_time(),
        ".pdf"
      )
    )
    return(p_cluster)
  }
colnames(combined_seurat@meta.data)
str(combined_seurat@meta.data$RNA_snn_scvi_leiden_res.2)

p_cluster_muti_louvain <-
  tree_cluster(combined_seurat,
    prefix = "RNA_snn_scvi_leiden_res.",
    width = 25,
    height = 20
  )
#这里命名没改过来，其实是leiden
p_cluster_muti_louvain

save(
  list = ls(pattern = "p_cluster"),
  file = paste0(
    "plotobject/cluster_scvi_leiden",
    get_time(),
    ".RData"
  )
)
