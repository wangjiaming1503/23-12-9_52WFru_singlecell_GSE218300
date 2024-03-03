
# SCTtransform====完成了常规标准化下的所有整合方法--------------------------------------------------------
combined_seurat
options(future.globals.maxSize = 3e+09)
combined_seurat <- SCTransform(combined_seurat)
combined_seurat <- RunPCA(combined_seurat, verbose = F)




combined_seurat <- IntegrateLayers(
  object = combined_seurat,
  method = RPCAIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca",
  new.reduction = "integrated.sct.rpca",
  verbose = F
)
save(
  combined_seurat,
  file = paste0(
    save.dir,
    "/combined_seurat_combined_umap_N4--rpca_harmony_cca_rcpa-mnn-SCT-N1-rpca",
    Sys.Date(),
    ".RData"
  )
)
# hamony
combined_seurat <- IntegrateLayers(
  object = combined_seurat,
  method = HarmonyIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca",
  new.reduction = "integrated.sct.harmony",
  verbose = F
)

save(
  combined_seurat,
  file = paste0(
    save.dir,
    "/combined_seurat_combined_umap_N4--rpca_harmony_cca_rcpa-mnn-SCT-N2-harmony",
    Sys.Date(),
    ".RData"
  )
)



# cca
combined_seurat <- IntegrateLayers(
  object = combined_seurat,
  method = CCAIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca",
  new.reduction = "integrated.sct.cca",
  verbose = F
)

save(
  combined_seurat,
  file = paste0(
    save.dir,
    "/combined_seurat_combined_umap_N4--rpca_harmony_cca_rcpa-mnn-SCT-N3-cca",
    Sys.Date(),
    ".RData"
  )
)



# # FastMNN
# combined_seurat <- IntegrateLayers(
#   object = combined_seurat,
#   method = FastMNNIntegration,
#   normalization.method = "SCT",orig.reduction = "pca", new.reduction = "integrated.sct.mnn",
#   verbose = F
# )

# save(combined_seurat,file = paste0(save.dir,"/combined_seurat_combined_umap_N4--rpca_harmony_cca_rcpa-mnn-SCT-N4-mnn", Sys.Date(), ".RData"))
# # scVI
# combined_seurat <- IntegrateLayers(
#   object = combined_seurat,
#   method = scVIIntegration,
#   normalization.method = "SCT",orig.reduction = "pca", new.reduction = "integrated.sct.scvi",
#   verbose = F
# )

# scVI
# combined_seurat <- IntegrateLayers(
#   object = combined_seurat,
#   method = scVIIntegration,
#   normalization.method = "SCT",orig.reduction = "pca", new.reduction = "integrated.sct.scvi",
#   verbose = F,conda_env = "/home/rstudio/work/miniconda3/envs/scvi-env-2"
# )
#
# save(combined_seurat,file = paste0(save.dir,"/combined_seurat_combined_umap_N4--rpca_harmony_cca_rcpa-mnn-SCT-N5-scvi", Sys.Date(), ".RData"))

# SCT转换后，整合用了harmony，cca两种，接下来是探究了Findcluster的不同算法，默认为1louvein，
# 后续2，3，4，slm,mitilouvein,leiden做了尝试，不同分辨率也尝试0.1-2
combined_seurat <-
  FindNeighbors(combined_seurat, dims = 1:30, reduction = "integrated.sct.rpca")
combined_seurat <-
  FindClusters(combined_seurat,
    resolution = 2,
    cluster.name = "sct_rcca_clusters"
  )
combined_seurat <-
  RunUMAP(
    combined_seurat,
    reduction = "integrated.sct.rpca",
    dims = 1:30,
    reduction.name = "umap.sct.rpca"
  )

p1 <- DimPlot(
  combined_seurat,
  reduction = "umap.sct.rpca",
  group.by = c("orig.ident", "sct_cca_clusters"),
  combine = FALSE,
  label.size = 2
)


#
p2 <- DimPlot(
  combined_seurat,
  reduction = "umap.sct.harmony",
  group.by = c("orig.ident", "sct_harmony_clusters"),
  combine = FALSE, label.size = 2
)

# 一些旧代码，应该是为了比一下整合方法的差异，能不能把来源的细胞都整合起来，
# 不过可以用louprbower，这里P1,2还是不太好，应该命名清楚，这里比较harmony和cca整合对样本的样子。
library(patchwork)
wrap_plots(c(p1, p2), ncol = 2, byrow = F)
library(ggplot2)
# save.dir <-
ggsave(
  paste0(
    "combined_seurat_combined_umap_N4--rpca_harmony_cca_rcpa-mnn-SCT-N2-cca-harmony",
    Sys.Date(),
    ".png"
  ),
  width = 10,
  height = 10
)



# 此处FindNeighBors,把后续的cluster的纬度换为了harmony整合。
combined_seurat <-
  FindNeighbors(combined_seurat, dims = 1:30, reduction = "integrated.sct.harmony")