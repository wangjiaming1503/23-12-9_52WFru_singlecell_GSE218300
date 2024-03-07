
combined_seurat <- load(".Rdata")
library(ggplot2)
# # 我想跑tSNE
combined_seurat <- RunTSNE(combined_seurat)
# # 画tSNE图
tryCatch(
  {
    save(combined_seurat, file = "SeuratObjects/combined_seurat_tsne_2024-02-24.RData")
  },
  error = function(e) {
    print(e)
  }
)

save(combined_seurat, file = "SeuratObjects/combined_seurat_tsne_2024-02-24.RData")
DimPlot(combined_seurat, reduction = "tsne", group.by = "seurat_clusters")
save(combined_seurat, file = "SeuratObjects/combined_seurat_tsne_2024-02-24.RData")
# # convert the SeuratObject named `seurat_obj` to a Loupe file


combined_seurat <-
  RunUMAP(
    combined_seurat,
    dims = 1:30,
    reduction = "pca",
    reduction.name = "umap.unintegrated"
  )

combined_seurat <-
  RunUMAP(
    combined_seurat,
    reduction = "integrated.cca",
    dims = 1:30,
    reduction.name = "umap.cca")


combined_seurat <-
  RunUMAP(
    combined_seurat,
    reduction = "harmony",
    dims = 1:30,
    reduction.name = "umap.harmony"
  )

combined_seurat <-
  RunUMAP(
    combined_seurat,
    reduction = "integrated.mnn",
    dims = 1:30,
    reduction.name = "umap.mnn"
  )

combined_seurat <-
  RunUMAP(
    combined_seurat,
    reduction = "integrated.sct.rpca",
    dims = 1:30,
    reduction.name = "umap.sct.rpca"
  )
# 保存了harmony的SCTumap，不过可以看看有啥参数，这个可能是设定了dims
combined_seurat <-
  RunUMAP(
    combined_seurat,
    reduction = "integrated.sct.harmony",
    dims = 1:30,
    reduction.name = "umap.2.sct.harmony"
  )
combined_seurat <-
  RunUMAP(
    combined_seurat,
    reduction = "integrated.sct.cca",
    dims = 1:30,
    reduction.name = "umap.sct.cca"
  )

#这个第一次用SCT的默认跑的umap，可以看到用了CCA，分辨率0.6---------------------------------------------------

combined_seurat <-
  IntegrateLayers(
    object = combined_seurat,
    method = CCAIntegration,
    normalization.method = "SCT",
    verbose = F
  )

# Find neighbors and clusters
combined_seurat <-
  FindNeighbors(combined_seurat, reduction = "integrated.dr", dims = 1:30) %>%
  FindClusters(resolution = 0.6)

# Run UMAP again on the integrated data
combined_seurat <-
  RunUMAP(combined_seurat, dims = 1:30, reduction = "integrated.dr")
# Visualize UMAP results post integration
p4 <-
  DimPlot(combined_seurat,
    reduction = "umap",
    group.by = c("orig.ident", "group")
  )
p4
# Compare UMAP plots and save them with the current date in the filename
umap_compare_plot <-  p4
umap_compare_plot

ggsave(
  paste0("umap_merged_intergrated_", Sys.Date(), ".png"),
  umap_compare_plot,
  width = 10,
  height = 10
)
save(combined_seurat,
  file = paste0("combined_seurat_", Sys.Date(), ".RData")
)

#后续探索一下参数