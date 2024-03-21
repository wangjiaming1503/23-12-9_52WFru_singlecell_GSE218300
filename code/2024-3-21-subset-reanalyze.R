library(Seurat)
library(tidyverse)
library(SeuratWrappers)
options(future.globals.maxSize = 1e9)
library(sceasy)
library(reticulate)
library(anndata)
load("SeuratObjects/2024-03-21-11-23anotated-alltsne-combined_seurat-quality_control_doublet_normalized_addmeta.RData")
# 将 outlier 转换为因子型
combined_seurat@meta.data$outlier <- factor(combined_seurat@meta.data$outlier, levels = c("TRUE", "FALSE"))

# 将 mt_outlier 转换为因子型 
combined_seurat@meta.data$mt_outlier <- factor(combined_seurat@meta.data$mt_outlier, levels = c("TRUE", "FALSE"))
combined_seurat
DefaultAssay(combined_seurat) 
colnames(combined_seurat@meta.data)
DefaultAssay(combined_seurat)  <- "RNA"
colnames(combined_seurat@meta.data)
DefaultAssay(combined_seurat)  <- "SCT"
colnames(combined_seurat@meta.data)
DefaultAssay(combined_seurat)  <- "RNA"

#subset后就不容易切换assay了，也许有同时subset的方法，应该是一堆在内存里搞错了
combined_seurat_filter <-
  subset(combined_seurat, subset = outlier == "FALSE" &
    mt_outlier == "FALSE" & doublet_class == "singlet")
combined_seurat_filter
colnames(combined_seurat_filter@meta.data)
DefaultAssay(combined_seurat_filter)  <- "SCT"
combined_seurat_filter
colnames(combined_seurat_filter@meta.data)
DefaultAssay(combined_seurat_filter)  <- "RNA"
combined_seurat_filter

library(Seurat)
# library(SeuratData)
library(SeuratWrappers)
# library(Azimuth)
library(ggplot2)
library(patchwork)
options(future.globals.maxSize = 1e9)

combined_seurat_filter <- NormalizeData(combined_seurat_filter)
combined_seurat_filter <- FindVariableFeatures(combined_seurat_filter)
combined_seurat_filter <- ScaleData(combined_seurat_filter)
combined_seurat_filter <- RunPCA(combined_seurat_filter)


combined_seurat_filter <-
  FindNeighbors(combined_seurat_filter, dims = 1:30, reduction = "pca")
combined_seurat_filter <-
  FindClusters(combined_seurat_filter,
    resolution = 2,
    cluster.name = "ft_unintegrated_clusters"
  )

combined_seurat_filter <-
  RunUMAP(
    combined_seurat_filter,
    dims = 1:30,
    reduction = "pca",
    reduction.name = "ft_umap.unintegrated"
  )
# visualize by batch and cell type annotation

Layers(combined_seurat_filter)

combined_seurat_filter[["RNA"]] <- split(combined_seurat_filter[["RNA"]], f = combined_seurat_filter$orig.ident)
Layers(combined_seurat_filter)

combined_seurat_filter <- IntegrateLayers(
  object = combined_seurat_filter,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "ft_integrated.cca",
  verbose = FALSE
)


combined_seurat_filter <- IntegrateLayers(
  object = combined_seurat_filter,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "ft_integrated.rpca",
  verbose = FALSE
)

save.dir <- "SeuratObjects"
# 保存临时整合结果，整合layers-CCA-RPCA
save(
  combined_seurat_filter,
  file = paste0(
    save.dir,get_time(),
    "/combined_seurat_filter_combinedN2_",
    Sys.Date(),
    ".RData"
  )
)

# 整合layers-harmony
combined_seurat_filter <- IntegrateLayers(
  object = combined_seurat_filter,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "ft_integrated.harmony",
  verbose = FALSE

)

# 保存临时整合结果，整合layers-CCA-RPCA-harmony
save(
  combined_seurat_filter,
  file = paste0(
    save.dir,
    "/combined_seurat_filter_combinedN3_cca-rcca-harmony",
    Sys.Date(),
    ".RData"
  )
)

# 整合layers-mmn
combined_seurat_filter <- IntegrateLayers(
  object = combined_seurat_filter,
  method = FastMNNIntegration,
  new.reduction = "ft_integrated.mnn",
  verbose = FALSE
)
# 保存临时整合结果，整合layers-CCA-RPCA-harmony-mmn
save(
  combined_seurat_filter,
  file = paste0(
    save.dir,"/",get_time(),
    "combined_seurat_filter_combinedN4_mnn",
    ".RData"
  )
)

use_condaenv("/home/rstudio/work/miniforge3/envs/scvi-env-4")
Sys.which("python")
combined_seurat_filter <- IntegrateLayers(
  object = combined_seurat_filter,
  method = scVIIntegration,
  new.reduction = "ft_integrated.scvi",
  conda_env = "/home/rstudio/work/miniforge3/envs/scvi-env-4",
  verbose = FALSE
)
save(
  combined_seurat_filter,
  file = paste0(
    save.dir,"/",get_time(),
    "combined_seurat_filter_combinedN5_scvi",
    ".RData"
  )
)


ElbowPlot(combined_seurat_filter,reduction = "pca",ndims = 50)

print(combined_seurat_filter[["pca"]], dims = 1:50, nfeatures = 10)
DimHeatmap(combined_seurat_filter, dims = 1:30, cells = 500, balanced = TRUE)
#DimHeatmap(combined_seurat_filter, dims = 1:40, cells = 500, balanced = TRUE)

#cluster
combined_seurat_filter <- FindNeighbors(combined_seurat_filter, dims = 1:30, reduction = "ft_integrated.scvi")
#combined_seurat_filter <- FindClusters(combined_seurat_filter, resolution = 2, cluster.name = "ft_integrated_clusters")

combined_seurat_filter <- RunUMAP(combined_seurat_filter, dims = 1:30, reduction = "ft_integrated.scvi", reduction.name = "ft_umap.scvi")
#tsne
combined_seurat_filter <- RunTSNE(combined_seurat_filter, dims = 1:30, reduction = "ft_integrated.scvi", reduction.name = "ft_tsne.scvi")


Reselute_clusters_slm <- function(combined_seurat, resolution) {
  combined_seurat <-
    FindClusters(combined_seurat,
      resolution = resolution,
      algorithm = 3
    )
  return(combined_seurat)
}
resolution  <- seq(0.1, 2, by = 0.1)
combined_seurat_filter <- Reselute_clusters_slm(combined_seurat_filter, resolution)
#合并
combined_seurat_filter[["RNA"]]  <- JoinLayers(combined_seurat_filter[["RNA"]])
library(loupeR)
create_loupe_from_seurat(
  combined_seurat_filter,
  output_name = paste0(
    save.dir,"/",get_time(),
    "combined_seurat_filter_combinedN6_reselute_clusters_slm"
  )
)
colnames(combined_seurat_filter@meta.data)

# 通过gusb函数替换含有"SCT_snn_res"的列名
colnames(combined_seurat_filter@meta.data) <-
  gsub(
    "RNA_snn_res.",
    "ft_RNA_snn_scvi_slm_res.",
    colnames(combined_seurat_filter@meta.data)
  )
colnames(combined_seurat_filter@meta.data)

save(
  combined_seurat_filter,
  file = paste0(
    save.dir,"/",get_time(),
    "combined_seurat_filter_combinedN6_reselute_clusters_slm",
    ".RData"
  )
)

# List of resolution values to loop over
res_values <-resolution
  
# Loop over each resolution value
for (res in res_values) {
  col_name <- paste0("ft_RNA_snn_scvi_slm_res.", res)
  
  # Check if the column exists in the metadata
  if (col_name %in% names(combined_seurat_filter@meta.data)) {
    # Extract the factor and store as a character vector
    metadata_col <-
      FetchData(combined_seurat_filter, vars = col_name)[, col_name]
    
    # Convert the factor levels to numeric, sort them, and convert back to a factor with the new levels
    combined_seurat_filter <- Seurat::AddMetaData(
      combined_seurat_filter,
      metadata = factor(metadata_col, levels = sort(as.numeric(
        levels(metadata_col)
      ))),
      col.name = col_name
    )
  }
}
str(combined_seurat_filter@meta.data)

tree_cluster <-
  function(combined_seurat,
           prefix,
           width = 25,
           height = 16) {
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
        Sys.Date(),
        ".pdf"
      )
    )
    return(p_cluster)
  }

p_cluster_muti_louvain <- tree_cluster(combined_seurat_filter,
                                       prefix = "ft_RNA_snn_scvi_slm_res.",
                                       width = 25,
                                       height = 16
)
