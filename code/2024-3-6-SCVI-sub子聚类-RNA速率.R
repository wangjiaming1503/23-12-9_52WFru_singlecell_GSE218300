library(Seurat)
library(SeuratWrappers)
library(tidyverse)
setwd("/home/rstudio/work/23-12-9_52WFru_singlecell_GSE218300")
View(maartenutils::gen_file_overview("SeuratObjects"))
load("SeuratObjects/combined_seurat_combined_umap_N4_harmony_cca_rcpa-mnn-keepsct-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden_all_joined2024-03-07-08-46.RData")

get_time()

# 优化无rstudio的使用环境。

list.files(path = "SeuratObjects", full.names = TRUE)
View(maartenutils::gen_file_overview("SeuratObjects"))
View(maartenutils::gen_file_overview("SeuratObjects"))
detach("package:maartenutils", unload = TRUE)
ls()
str(combined_seurat, max.level = 3)

names(combined_seurat)
slotNames(combined_seurat)

combined_seurat@reductions$integrated.scvi

combined_seurat <- RunTSNE(
  combined_seurat,
  reduction = "integrated.scvi",
  cells = NULL,
  dims = 1:30,
  features = NULL,
  seed.use = 1,
  tsne.method = "Rtsne",
  dim.embed = 2,
  distance.matrix = NULL,
  reduction.name = "tsne.scvi",
  reduction.key = "tSNE_"
)

DimPlot(
  combined_seurat,
  dims = c(1, 2),
  cells = NULL,
  cols = NULL,
  pt.size = NULL,
  reduction = "tsne.scvi",
  group.by = c(
    "RNA_snn_scvi_leiden_res.2",
    "ImmGen_SingleR.labels",
    "MouseRNA_SingleR.labels"
  ),
  split.by = "group",
  shape.by = NULL,
  order = NULL,
  shuffle = FALSE,
  seed = 1,
  label = TRUE,
  label.size = 4,
  label.color = "black",
  label.box = FALSE,
  repel = TRUE,
  alpha = 0.6,
  cells.highlight = NULL,
  cols.highlight = "#DE2D26",
  sizes.highlight = 1,
  na.value = "grey50",
  ncol = 2,
  combine = TRUE,
  raster = NULL,
  raster.dpi = c(512, 512)
)

# umap.learn.scvi
DimPlot(
  combined_seurat,
  dims = c(1, 2),
  cells = NULL,
  cols = NULL,
  pt.size = NULL,
  reduction = "umap.learn.scvi",
  group.by = c(
    "RNA_snn_scvi_leiden_res.2",
    "ImmGen_SingleR.labels",
    "MouseRNA_SingleR.labels"
  ),
  split.by = "group",
  shape.by = NULL,
  order = NULL,
  shuffle = FALSE,
  seed = 1,
  label = TRUE,
  label.size = 4,
  label.color = "black",
  label.box = FALSE,
  repel = TRUE,
  alpha = 0.6,
  cells.highlight = NULL,
  cols.highlight = "#DE2D26",
  sizes.highlight = 1,
  na.value = "grey50",
  ncol = 2,
  combine = TRUE,
  raster = NULL,
  raster.dpi = c(512, 512)
)


str(combined_seurat, max.level = 3)
colnames(combined_seurat@meta.data)

save.dir <- "SeuratObjects"
View(maartenutils::gen_file_overview("SeuratObjects"))

saveRDS(combined_seurat, paste0(
  save.dir,
  "/", Sys.Date(),
  "combined_seurat-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden_all_joined-tsne",
  get_time(),
  ".rds"
),
ascii = FALSE, version = NULL,
compress = TRUE, refhook = NULL
)

save(combined_seurat,
  file = paste0(
    save.dir,
    "/", Sys.Date(),
    "combined_seurat-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden_all_joined-tsne",
    get_time(),
    ".RData"
  )
)

library("loupeR")
loupeR::setup()

create_loupe_from_seurat(
  combined_seurat,
  output_name = paste0(
    save.dir,
    "/", Sys.Date(),
    "combined_seurat-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden_all_joined-tsne",
    get_time()
  )
)
library(httpgd)
hgd_browse()

features <- c("Lcn2")

str(combined_seurat@meta.data)
load("SeuratObjects/2024-03-08combined_seurat-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden_all_joined-tsne2024-03-08-13-38.RData")
FeaturePlot(
  combined_seurat,
  features,
  alpha = 0.6,
  reduction = "umap.learn.scvi",
  split.by = "group"
)

#用Nebulosa包，只能一个组，可以看代表基因在哪个细胞，适合低表达突出重点，参数有feature和降维--------------------------
library(Nebulosa)
#install.packages("Nebulosa")
plot_density(combined_seurat, "Spp1",reduction = "tsne.scvi")




#RNA速率分化，导出到python-----------------------------
setwd("/home/rstudio/work/23-12-9_52WFru_singlecell_GSE218300")

library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)

DefaultAssay(combined_seurat) <- "RNA"
save.dir  <- "SeuratObjects"
SaveH5Seurat(combined_seurat, filename = paste0(save.dir,
    "/", get_time(),"combined_seurat-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden_all_joined-tsne.h5Seurat"))

#这部分作为导出和转换，可以看看其他的scanpy分析
View(maartenutils::gen_file_overview("SeuratObjects"))
Convert("SeuratObjects/2024-03-08-18-02combined_seurat-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden_all_joined-tsne.h5Seurat", dest = "h5ad")






#R中的速率分析这部分暂时运行不了，需要从BAM生成splice，unsplice----------------------------------------
library(Seurat)
BiocManager::install("pcaMethods")
#sudo apt-get install libboost-all-dev
library(devtools)
#install_github("velocyto-team/velocyto.R",dependencies = )
library(velocyto.R)
library(SeuratWrappers)

bm <- RunVelocity(object = bm, deltaT = 1, kCells = 25, fit.quantile = 0.02)
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = bm)))
names(x = ident.colors) <- levels(x = bm)
cell.colors <- ident.colors[Idents(object = bm)]
names(x = cell.colors) <- colnames(x = bm)
show.velocity.on.embedding.cor(emb = Embeddings(object = bm, reduction = "umap"), vel = Tool(object = bm, 
    slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
    cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
    do.par = FALSE, cell.border.alpha = 0.1)

    