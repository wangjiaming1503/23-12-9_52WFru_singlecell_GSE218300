
# 整合看效果，前面已经有不整合的，步骤就是多一个整合layers---------------------
# 整合layers-CCA-RPCA------------------------
combined_seurat <- IntegrateLayers(
  object = combined_seurat,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE
)


combined_seurat <- IntegrateLayers(
  object = combined_seurat,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  verbose = FALSE
)
# 保存临时整合结果，整合layers-CCA-RPCA
save(
  combined_seurat,
  file = paste0(
    save.dir,
    "/combined_seurat_combinedN2_",
    Sys.Date(),
    ".RData"
  )
)

# 整合layers-harmony
combined_seurat <- IntegrateLayers(
  object = combined_seurat,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = FALSE
)
# 保存临时整合结果，整合layers-CCA-RPCA-harmony
save(
  combined_seurat,
  file = paste0(
    save.dir,
    "/combined_seurat_combinedN3_cca-rcca-harmony",
    Sys.Date(),
    ".RData"
  )
)

# 整合layers-mmn
combined_seurat <- IntegrateLayers(
  object = combined_seurat,
  method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)
# 保存临时整合结果，整合layers-CCA-RPCA-harmony-mmn
save(
  combined_seurat,
  file = paste0(
    save.dir,
    "/combined_seurat_combinedN4_mnn",
    Sys.Date(),
    ".RData"
  )
)

# SCVI包安装准备，这部分其实是SCVI用在R中的教程，不过SCVI有bug====================================

# install.packages("Seurat")
# install.packages("reticulate")
# install.packages("cowplot")
# install.packages("devtools")
# devtools::install_github("satijalab/seurat-data")
# SeuratData::InstallData("pbmc3k")
# install.packages(
#   "https://seurat.nygenome.org/src/contrib/ifnb.SeuratData_3.0.0.tar.gz",
#   repos = NULL,
#   type = "source"
# )
# SeuratData::InstallData("ifnb")
# devtools::install_github("cellgeni/sceasy")


load(
  "SeuratObjects/combined_seurat_combined_umap_N4--SCT-N3-cca-N3-louvain_slm_muti2024-03-02.RData"
)

# library(reticulate)

# 这可以指定python
# path_to_python <- "/home/rstudio/work/miniconda3/envs/scvi-env-2/bin/python"
# use_python(path_to_python)
library(Seurat)
library(SeuratWrappers)
# 整合前设定最大的对象可以避免报错
options(future.globals.maxSize = 1e9)

# library(reticulate)
# library(anndata)
# install.packages("anndata")
# library(sceasy)

#
# # library(SeuratData)
# remotes::install_github('satijalab/seurat-wrappers')
# BiocManager::install(c("LoomExperiment", "SingleCellExperiment"))
# devtools::install_github("cellgeni/sceasy")

# conda用这个方法会好一点最好先有conda指定再查询python位置，不然容易锁定
# Sys.which("python")
library(reticulate)
use_condaenv("/home/rstudio/work/anaconda3/envs/scvi-env")
# 尝试了机械上的bioinfo的anaconda
combined_seurat <- IntegrateLayers(
  object = combined_seurat,
  method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "../anaconda3/envs/scvi-env",
  verbose = FALSE
)

# 尝试了机械上的bioinfo的miniconda
combined_seurat <- IntegrateLayers(
  object = combined_seurat,
  method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "/home/rstudio/work/miniconda3/envs/scvi-env",
  verbose = FALSE
)
# 如果成功了的话就保存命名为N5 SCVI
save.dir <- "SeuratObjects"
save(combined_seurat, file = paste0(save.dir, "/combined_seurat_combinedN5_scvi", Sys.Date(), ".RData"))

# 用cca整合来找neighbor和cluster-获得了cca-cluster-分辨率为2的聚类，名字为cca_clusters

# 整合确定reduction_
combined_seurat <-
  FindNeighbors(combined_seurat, reduction = "integrated.cca", dims = 1:30)
combined_seurat <-
  FindClusters(combined_seurat,
    resolution = 2,
    cluster.name = "cca_clusters"
  )

# 1 cca整合来画了umap_cca,这里画umap图，然后用来源来看整合效果，此处是在确定要用什么整合的方法-------------------------------
combined_seurat <-
  RunUMAP(
    combined_seurat,
    reduction = "integrated.cca",
    dims = 1:30,
    reduction.name = "umap.cca"
  )
p1 <- DimPlot(
  combined_seurat,
  reduction = "umap.cca",
  group.by = c("orig.ident", "cca_clusters"),
  combine = FALSE,
  label.size = 2
)

# combined_seurat <- FindNeighbors(combined_seurat, reduction = "integrated.scvi", dims = 1:30)
# combined_seurat <- FindClusters(combined_seurat, resolution = 2, cluster.name = "scvi_clusters")
# combined_seurat <- RunUMAP(combined_seurat, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi")
# p2 <- DimPlot(
#   combined_seurat,
#   reduction = "umap.scvi",
#   group.by = c("orig.ident", "scvi_clusters"),
#   combine = FALSE, label.size = 2
# )

# 2 rpca整合---------------------------
combined_seurat <-
  FindNeighbors(combined_seurat, reduction = "integrated.rpca", dims = 1:30)
combined_seurat <-
  FindClusters(combined_seurat,
    resolution = 2,
    cluster.name = "rpca_clusters"
  )
combined_seurat <-
  RunUMAP(
    combined_seurat,
    reduction = "integrated.rpca",
    dims = 1:30,
    reduction.name = "umap.rpca"
  )
p3 <- DimPlot(
  combined_seurat,
  reduction = "umap.rpca",
  group.by = c("orig.ident", "rpca_clusters"),
  combine = FALSE,
  label.size = 2
)

# 3 harmony整合---------------------------------
combined_seurat <-
  FindNeighbors(combined_seurat, reduction = "harmony", dims = 1:30)
combined_seurat <-
  FindClusters(combined_seurat,
    resolution = 2,
    cluster.name = "harmony_clusters"
  )
combined_seurat <-
  RunUMAP(
    combined_seurat,
    reduction = "harmony",
    dims = 1:30,
    reduction.name = "umap.harmony"
  )
p4 <- DimPlot(
  combined_seurat,
  reduction = "umap.harmony",
  group.by = c("orig.ident", "harmony_clusters"),
  combine = FALSE,
  label.size = 2
)

# 4 mnn-------------------------------------
combined_seurat <-
  FindNeighbors(combined_seurat, reduction = "integrated.mnn", dims = 1:30)
combined_seurat <-
  FindClusters(combined_seurat,
    resolution = 2,
    cluster.name = "mnn_clusters"
  )
combined_seurat <-
  RunUMAP(
    combined_seurat,
    reduction = "integrated.mnn",
    dims = 1:30,
    reduction.name = "umap.mnn"
  )
p5 <- DimPlot(
  combined_seurat,
  reduction = "umap.mnn",
  group.by = c("orig.ident", "mnn_clusters"),
  combine = FALSE,
  label.size = 2
)

library(patchwork)
# p1 + p2 / p3 + p4 / p5
# p1 + p3 / p4 + p5
# 这里的图片形式因为采用了combine=FALSE，所以是两张图分来画出来的，结合
# wrap的话比较合适，后期采用loupe的方法查看就不需要画出来比较了

wrap_plots(c(p4, p5), ncol = 2, byrow = F) + NoLegend()


save(
  combined_seurat,
  file = paste0(
    save.dir,
    "/combined_seurat_combined_umap_N4_harmony_cca_rcpa-mnn",
    Sys.Date(),
    ".RData"
  )
)

wrap_plots(c(p4, p5), ncol = 2, byrow = F) + NoLegend()
ggsave(
  paste0(
    save.dir,
    "/combined_seurat_combined_umap_N4_harmony-mnn",
    Sys.Date(),
    ".png"
  ),
  width = 10,
  height = 10
)

wrap_plots(c(p1, p3), ncol = 2, byrow = F) + NoLegend()

p1
ggsave(
  paste0(
    save.dir,
    "/combined_seurat_combined_umap_N4_cca-rcpa",
    Sys.Date(),
    ".png"
  ),
  width = 10,
  height = 10
)

# 前面采用了seruat的整合教程，用除了SCVI的方法做了整合，之前试着用SCTtransform后
# 再整合有bug，应该是启用了SCT之后就不容易回去了，而且如果join之后再spit也不太容易，
# 如果要做转换的话可以注意保存join之前的版本
