library(Seurat)
library(SeuratWrappers)
library(tidyverse)
setwd("/home/rstudio/work/23-12-9_52WFru_singlecell_GSE218300")
View(maartenutils::gen_file_overview("SeuratObjects"))
load("SeuratObjects/2024-03-08combined_seurat-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden_all_joined-tsne2024-03-08-13-38.RData")

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

#细胞比例分析==================================

library(Seurat)
library(ggplot2)
library(tidyverse)

# Extract cell type annotations from the metadata


cell_types <- combined_seurat@meta.data$cell_type

colnames(combined_seurat@meta.data)
combined_seurat@meta.data[[RNA_snn_scvi_leiden_res.2]]
library(mclust)  # for adjustedRandIndex
library(aricode) # for NMI
#install.packages("aricode")
# 提取聚类结果和已知细胞类型注释
metadata <- combined_seurat@meta.data
colnames(metadata)

adjustedRandIndex(metadata$SCT_snn_harmony_leiden_res.0.8,
metadata$MouseRNA_SingleR.labels)


cluster_cols <- grep("_clusters$|_res\\.", colnames(metadata), value = TRUE)

ari_scores <- sapply(cluster_cols, function(col) {
  adjustedRandIndex(metadata[[col]], metadata$MouseRNA_SingleR.labels)
})

ari_df <- data.frame(Cluster = cluster_cols, ARI = ari_scores)

ari_df <- ari_df[order(-ari_df$ARI), ]

library(ggplot2)
ggplot(ari_df, aes(x = reorder(Cluster, ARI), y = ARI)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 0.7) +
  coord_flip() +
  xlab("Clustering") +
  ylab("Adjusted Rand Index") +
  ggtitle("Comparison of Clustering Results with MouseRNA_SingleR.labels") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10))





# 计算 ARI
ari <- adjustedRandIndex(clustering, reference)

# 计算 NMI
nmi <- NMI(clustering, reference)

clusterings <- colnames(combined_seurat@meta.data[, -c(1:5, 12, 13, 18:22)])

references <- c("MouseRNA_SingleR.labels", "ImmGen_SingleR.labels","Blueprint_SingleR.labels",
"ImmGen_fine_SingleR.labels","HPA_SingleR.labels")


# Define the clusterings and references to evaluate
clusterings <- colnames(combined_seurat@meta.data[, -c(1:5, 12, 13, 18:22)])
references <- c("MouseRNA_SingleR.labels", "ImmGen_SingleR.labels",
                "Blueprint_SingleR.labels", "ImmGen_fine_SingleR.labels",
                "HPA_SingleR.labels")


compare_clusterings <- function(metadata, ref_col) {
  # 提取所有聚类结果的列名
  cluster_cols <- grep("_clusters$|_res\\.", colnames(metadata), value = TRUE)
  
  # 计算每个聚类结果与参考分组的ARI
  ari_scores <- sapply(cluster_cols, function(col) {
    adjustedRandIndex(metadata[[col]], metadata[[ref_col]])
  })
  
  # 将ARI得分与对应的聚类结果名称组合为数据框
  ari_df <- data.frame(Cluster = cluster_cols, ARI = ari_scores)
  
  # 按ARI得分降序排列数据框
  ari_df <- ari_df[order(-ari_df$ARI), ]
  
  # 导出CSV文件
  write.csv(ari_df, file = paste0(get_time(), ref_col, "_clustering_comparison.csv"), row.names = FALSE)
  
  # 绘制条形图
  library(ggplot2)
  p <- ggplot(ari_df, aes(x = reorder(Cluster, ARI), y = ARI)) +
    geom_bar(stat = "identity", fill = "skyblue", width = 0.7) +
    coord_flip() +
    xlab("Clustering") +
    ylab("Adjusted Rand Index") +
    ggtitle(paste("Comparison of Clustering Results with", ref_col)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 10))
  
  # 保存图像文件
  ggsave(paste0(get_time(),ref_col, "_clustering_comparison.pdf"), plot = p, ,units = "in",width = 10, height = 20, dpi = 1200)
}

# 应用函数到参考基因列表
ref_cols <- c("MouseRNA_SingleR.labels", "ImmGen_SingleR.labels",
              "Blueprint_SingleR.labels", "ImmGen_fine_SingleR.labels",
              "HPA_SingleR.labels")

for (ref_col in ref_cols) {
  compare_clusterings(metadata, ref_col)
}

#分配聚类结果名称，采用slmrse0.2-scvi==========================================
colnames(combined_seurat@meta.data)
# "RNA_snn_scvi_slm_res.0.2"
Idents(combined_seurat)  <-  "RNA_snn_scvi_slm_res.0.2"
table(Idents(combined_seurat))
table(combined_seurat@meta.data$ImmGen_SingleR.labels)
table(combined_seurat@meta.data$MouseRNA_SingleR.labels)

combined_seurat  <- RenameIdents(combined_seurat,
              "0" = "Endothelial cells")
table(Idents(combined_seurat))
combined_seurat  <- RenameIdents(combined_seurat,
              "1" = "Macrophages")

combined_seurat  <- RenameIdents(combined_seurat,
              "2" = "Hepatocytes-1",
              "3" = "Fibroblasts-1",
              "4" = "Hepatocytes-2",
              "5" = "Monocytes-1",
              "6" = "Hepatocytes-3",
              "7" = "Monocytes-2",
              "8" = "Neutrophils",
              "9" = "T cells",
              "10" = "B cells-1",
              "11" = "B cells-2",
              "12" = "Erythrocytes-1",
              "13" = "Fibroblasts-2",
              "14" = "Erythrocytes-2",
              "15" = "Monocytes/B cells",
              "17" = "Macrophages-2",

              "18" = "Macrophages-3")
table(Idents(combined_seurat))
head(Idents(combined_seurat))

RidgePlot(combined_seurat,features = "Lcn2",split.by  = "group")
VlnPlot(combined_seurat,features = "Spp1",split.by  = "group")
DotPlot(combined_seurat,features = "Lcn2")
FeaturePlot(combined_seurat,reduction = "tsne.scvi",features = "Lcn2",split.by  = "group")


DimPlot(
  combined_seurat,
  dims = c(1, 2),
  cells = NULL,
  cols = NULL,
  pt.size = NULL,
  reduction = "tsne.scvi",
  split.by = "group",
  shape.by = NULL,
  order = NULL,
  shuffle = FALSE,
  seed = 1,
  label = FALSE,
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
