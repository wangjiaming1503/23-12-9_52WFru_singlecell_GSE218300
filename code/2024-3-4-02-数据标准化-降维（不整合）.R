# 数据读取完毕，进行数据降维标准化，先尝试不整合看分析效果----------------------------
# Normalize data, find variable features, and perform linear dimension reduction
load(
  "~/work/23-12-9_52WFru_singlecell_GSE218300/SeuratObjects/combined_seurat_af_raw_merge_fillter_2024-03-01.RData"
)

library(Seurat)
# library(SeuratData)
library(SeuratWrappers)
# library(Azimuth)
library(ggplot2)
library(patchwork)
options(future.globals.maxSize = 1e9)

# combined_seurat[["RNA"]] <- split(combined_seurat[["RNA"]], f = combined_seurat$orig.ident)
combined_seurat

combined_seurat <- NormalizeData(combined_seurat)
combined_seurat <- FindVariableFeatures(combined_seurat)
combined_seurat <- ScaleData(combined_seurat)
combined_seurat <- RunPCA(combined_seurat)

combined_seurat <-
  FindNeighbors(combined_seurat, dims = 1:30, reduction = "pca")
combined_seurat <-
  FindClusters(combined_seurat,
    resolution = 2,
    cluster.name = "unintegrated_clusters"
  )

combined_seurat <-
  RunUMAP(
    combined_seurat,
    dims = 1:30,
    reduction = "pca",
    reduction.name = "umap.unintegrated"
  )
# visualize by batch and cell type annotation
# cell type annotations were previously added by Azimuth
DimPlot(
  combined_seurat,
  reduction = "integrated.cca",
  group.by = c("orig.ident", "seurat_clusters"),
  label = T
)