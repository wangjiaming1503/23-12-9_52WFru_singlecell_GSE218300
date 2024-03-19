library(Seurat)
# library(SeuratData)
library(SeuratWrappers)
# library(Azimuth)
library(ggplot2)
library(patchwork)
options(future.globals.maxSize = 1e9)

library(Seurat)
library(SeuratWrappers)
options(future.globals.maxSize = 1e9)
library(sceasy)
library(reticulate)
library(anndata)
library(ggplot2)


load(
  "~/work/23-12-9_52WFru_singlecell_GSE218300/SeuratObjects/2024-03-08combined_seurat-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden_all_joined-tsne2024-03-08-13-38.RData"
)
Layers(combined_seurat)
combined_seurat

combined_seurat[["RNA"]] <-
  split(
    combined_seurat[["RNA"]],
    f = combined_seurat$orig.ident,
    layers = c('counts', 'data', 'scale.data')
  )
combined_seurat[["SCT"]] <-
  split(
    combined_seurat[["SCT"]],
    f = combined_seurat$orig.ident,
    layers = c('counts', 'data', 'scale.data')
  )


combined_seurat <- JoinLayers(combined_seurat)
DefaultAssay(combined_seurat)

DefaultAssay(combined_seurat) <- "RNA"
DefaultAssay(combined_seurat) <- "SCT"
save(combined_seurat, file = "./SeuratObjects/combined_seurat-2024-03-18-23.RData")

combined_seurat
DefaultAssay(combined_seurat)
Layers(combined_seurat)

use_condaenv("/home/rstudio/work/miniforge3/envs/scvi-env-4")

Sys.which("python")

combined_seurat <- IntegrateLayers(
  object = combined_seurat,
  method = scVIIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.sct.scvi",
  conda_env = "/home/rstudio/work/miniforge3/envs/scvi-env-4",
  normalization.method = "SCT",
  verbose = FALSE
)
save.dir <- "SeuratObjects"
combined_seurat <-
  FindNeighbors(combined_seurat, dims = 1:30, reduction = "integrated.sct.scvi") #这里等同于RNA激活assey，

combined_seurat <- RunTSNE(
  combined_seurat,
  reduction = "integrated.sct.scvi",
  cells = NULL,
  dims = 1:30,
  features = NULL,
  seed.use = 1,
  tsne.method = "Rtsne",
  dim.embed = 2,
  distance.matrix = NULL,
  reduction.name = "tsne.sct.scvi",
  reduction.key = "sct.scvi.tSNE_"
)


combined_seurat <-
  RunUMAP(
    combined_seurat,
    dims = 1:30,
    reduction = "integrated.sct.scvi",
    reduction.name = "umap.sct.scvi"
  )


combined_seurat <-
  Reselute_clusters_slm(combined_seurat, seq(0.1, 2, by = 0.3))
combined_seurat <-
  Reselute_clusters_slm(combined_seurat, 2)

colnames(combined_seurat@meta.data)

# 通过gusb函数替换含有"SCT_snn_res"的列名
colnames(combined_seurat@meta.data) <-
  gsub("SCT_snn_res",
       "SCT_snn_scvi_slm_res",
       colnames(combined_seurat@meta.data))

# 确认替换成功
colnames(combined_seurat@meta.data)

combined_seurat[["RNA"]] <- JoinLayers(combined_seurat[["RNA"]])
combined_seurat[["SCT"]] <- JoinLayers(combined_seurat[["SCT"]])

library("loupeR")
create_loupe_from_seurat(
  combined_seurat,
  output_name = paste0(
    "combined_seurat-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden_all_joined-tsne-SCVI-sct_",
    format(Sys.Date(), "%Y-%m-%d")
  )
)

save(
  combined_seurat,
  file = paste0(
    save.dir,
    get_time(),
    "/combined_seurat-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden_all_joined-tsne-SCVI-sct",
    Sys.Date(),
    ".RData"
  )
)

#发现整合出来SCVI的降维还是按照RNA的降维，极为相似，所以只要defaultassey是RNA，那么整合的就是RNA。

# history
load(
  "~/work/23-12-9_52WFru_singlecell_GSE218300/SeuratObjects/combined_seurat_anotated_mouse_blu_immgen_hpa_v5_fine_tsne_joined_2024-02-24.RData"
)
# Load necessary libraries
library(Seurat)
View(combined_seurat)
reticulate::repl_python()
load(
  "~/work/23-12-9_52WFru_singlecell_GSE218300/SeuratObjects/2024-03-08combined_seurat-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden_all_joined-tsne2024-03-08-13-38.RData"
)
ls()
library(Seurat)
library(Seurat)
library(Seurat)
# library(SeuratData)
library(SeuratWrappers)
# library(Azimuth)
library(ggplot2)
library(patchwork)
options(future.globals.maxSize = 1e9)
Layers(combined_seurat)
combined_seurat[["RNA"]] <-
  split(combined_seurat[["RNA"]], f = combined_seurat$orig.ident)
DefaultAssay(combined_seurat)
DefaultAssay(combined_seurat) <- "SCT"
Layers(combined_seurat)
combined_seurat[["RNA"]] <-
  split(combined_seurat[["RNA"]], f = combined_seurat$orig.ident)
Layers(combined_seurat)
combined_seurat
load(
  "~/work/23-12-9_52WFru_singlecell_GSE218300/SeuratObjects/combined_seurat-2024-03-18-23.RData"
)
library(Seurat)
library(SeuratWrappers)
options(future.globals.maxSize = 1e9)
library(reticulate)
library(ggplot2)
use_condaenv("/home/rstudio/work/miniforge3/envs/scvi-env-4")
Layers(combined_seurat)
DefaultAssay(combined_seurat)
View(combined_seurat)
Sys.which("python")
combined_seurat <- IntegrateLayers(
  object = combined_seurat,
  method = scVIIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.sct.scvi",
  conda_env = "/home/rstudio/work/miniforge3/envs/scvi-env-4",
  verbose = FALSE,
  assay = "SCT"
)
combined_seurat[["SCT"]] <-
  split(
    combined_seurat[["SCT"]],
    f = combined_seurat$orig.ident,
    layers = c('counts', 'data', 'scale.data')
  )
Layers(combined_seurat)
DefaultAssay(combined_seurat)
combined_seurat <- IntegrateLayers(
  object = combined_seurat,
  method = scVIIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.sct.scvi",
  conda_env = "/home/rstudio/work/miniforge3/envs/scvi-env-4",
  verbose = FALSE,
  assay = "SCT"
)
rlang::last_trace()
rlang::last_trace(drop = FALSE)
combined_seurat <- IntegrateLayers(
  object = combined_seurat,
  method = scVIIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.sct.scvi",
  conda_env = "/home/rstudio/work/miniforge3/envs/scvi-env-4",
  normalization.method = "SCT",
  verbose = FALSE,
  assay = "SCT"
)
combined_seurat <- IntegrateLayers(
  object = combined_seurat,
  method = scVIIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.sct.scvi",
  conda_env = "/home/rstudio/work/miniforge3/envs/scvi-env-4",
  normalization.method = "SCT",
  verbose = FALSE
)
combined_seurat <- IntegrateLayers(
  object = combined_seurat,
  method = scVIIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.sct.scvi",
  conda_env = "/home/rstudio/work/miniforge3/envs/scvi-env-4",
  #normalization.method = "SCT",
  verbose = FALSE
)
DefaultAssay(combined_seurat) <- "RNA"
combined_seurat <- IntegrateLayers(
  object = combined_seurat,
  method = scVIIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.sct.scvi",
  conda_env = "/home/rstudio/work/miniforge3/envs/scvi-env-4",
  #normalization.method = "SCT",
  verbose = FALSE
)
combined_seurat <- IntegrateLayers(
  object = combined_seurat,
  method = scVIIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.sct.scvi",
  conda_env = "/home/rstudio/work/miniforge3/envs/scvi-env-4",
  normalization.method = "SCT",
  verbose = FALSE
)
combined_seurat
DefaultAssay(combined_seurat)
DefaultAssay(combined_seurat) <- "SCT"
combined_seurat
combined_seurat <-
  FindNeighbors(combined_seurat, dims = 1:30, reduction = "integrated.sct.scvi")
combined_seurat <- RunTSNE(
  combined_seurat,
  reduction = "integrated.sct.scvi",
  cells = NULL,
  dims = 1:30,
  features = NULL,
  seed.use = 1,
  tsne.method = "Rtsne",
  dim.embed = 2,
  distance.matrix = NULL,
  reduction.name = "tsne.sct.scvi",
  reduction.key = "sct.scvi.tSNE_"
)
combined_seurat <-
  RunUMAP(
    combined_seurat,
    dims = 1:30,
    reduction = "tsne.sct.scvi",
    reduction.name = "umap.sct.scvi"
  )
combined_seurat <-
  RunUMAP(
    combined_seurat,
    dims = 1:30,
    reduction = "integrated.sct.scvi",
    reduction.name = "umap.sct.scvi"
  )
# 这是后面加进去的，不同的cluster算法如何避免记住参数，把其包装出来-------------------------------------------
Reselute_clusters_slm <- function(combined_seurat, resolution) {
  combined_seurat <-
    FindClusters(combined_seurat,
                 resolution = resolution,
                 algorithm = 3)
  return(combined_seurat)
}
seq(0.5, 2, by = 0.3)
seq(0.1, 2, by = 0.3)
Reselute_clusters_louvain_muti <-
  function(combined_seurat, resolution) {
    combined_seurat <-
      FindClusters(combined_seurat,
                   resolution = resolution,
                   algorithm = 2)
    return(combined_seurat)
  }
Reselute_clusters_leiden <- function(combined_seurat, resolution) {
  combined_seurat <-
    FindClusters(combined_seurat,
                 resolution = resolution,
                 algorithm = 4)
  return(combined_seurat)
}
combined_seurat <-
  Reselute_clusters_slm(combined_seurat, 2)
colnames(combined_seurat@meta.data)
# 通过gusb函数替换含有"SCT_snn_res"的列名
colnames(combined_seurat@meta.data) <-
  gsub("SCT_snn_res",
       "SCT_snn_scvi_slm_res",
       colnames(combined_seurat@meta.data))
# 确认替换成功
colnames(combined_seurat@meta.data)
library("loupeR")
loupeR::setup()
create_loupe_from_seurat(
  combined_seurat,
  output_name = paste0(
    "combined_seurat-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden_all_joined-tsne-SCVI-sct_",
    format(Sys.Date(), "%Y-%m-%d")
  )
)
combined_seurat[["SCT"]] <- JoinLayers(combined_seurat[["SCT"]])
Layers(combined_seurat)
DefaultAssay(combined_seurat) <- "RNA"
Layers(combined_seurat)
DefaultAssay(combined_seurat) <- "SCT"
create_loupe_from_seurat(
  combined_seurat,
  output_name = paste0(
    "combined_seurat-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden_all_joined-tsne-SCVI-sct_",
    format(Sys.Date(), "%Y-%m-%d")
  )
)
