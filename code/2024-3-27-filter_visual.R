load(
  "~/work/23-12-9_52WFru_singlecell_GSE218300/SeuratObjects/2024-03-26-18-41combined_seurat_filter_filter_combinedN6_reselute_clusters_slm.RData"
)
# failed seurat disk to anndata---------------------------
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(SeuratDisk)
# Load Packages
library(Seurat)
library(rliger)
library(scCustomize)
library(qs)
library(scran)
library(scater)
sce <- as.SingleCellExperiment(combined_seurat_filter_filter)

sce
combined_seurat_filter_filter



plotReducedDim(sce, dimred="FT_UMAP.SCVI",colour_by = "ft_RNA_snn_scvi_slm_res.2")
SingleCellExperiment::reducedDimNames(sce)
SingleCellExperiment::colData(sce)
colnames(SingleCellExperiment::colData(sce))

#分配聚类结果名称，采用ft_RNA_snn_scvi_louv_res.0.5==========================================
colnames(combined_seurat_filter@meta.data)
# "RNA_snn_scvi_slm_res.0.2"
Idents(combined_seurat_filter)  <-  "ft_RNA_snn_scvi_louv_res.0.5"
table(Idents(combined_seurat_filter))
table(combined_seurat_filter@meta.data$ImmGen_SingleR.labels)
table(combined_seurat_filter@meta.data$MouseRNA_SingleR.labels)

# combined_seurat_filter  <- RenameIdents(combined_seurat_filter,
#               "0" = "Endothelial cells",
#               "8" = "Endothelial cells",
#               "15"= "Endothelial cells")
table(Idents(combined_seurat_filter))
Idents(combined_seurat_filter)  <-  "ft_RNA_snn_scvi_louv_res.0.5"
combined_seurat_filter  <- RenameIdents(combined_seurat_filter,
                                "0" = "Endothelial cells",
                                "1" = "Mononucl. phagocytes",
                                 "2" = "Mononucl. phagocytes",
                                 "3" = "Mesenchymal cells",
                                 "4" = "Hepatocytes",
                                 "5" = "Hepatocytes",
                                 "6" = "Hepatocytes",
                                 "7" = "Cholangiocyte-like",
                                 "8" = "Endothelial cells",
                                 "9" = "Mononucl. phagocytes",
                                 "10" = "Neutrophils",
                                 "11" = "T cells",
                                 "12" = "Mononucl. phagocytes",
                                 "13" = "Mononucl. phagocytes",
                                 "14" = "B cells",
                                 "15" = "Endothelial cells",
                                 "16" = "Mononucl. phagocytes",
                                 "17" = "T cells",
                                 
                                 "18" = "B cells",
                                 "19" = "Erythrocytes",
                                  "20" = "Mesenchymal cells",
                                  "21" = "NK cells",
                                  "22" = "Mononucl. phagocytes",
                                  "23" = "Mononucl. phagocytes",
                                  "24" = "Mononucl. phagocytes",
                                  "25" = "Mononucl. phagocytes",
                                  "26" = "Mononucl. phagocytes",
                                  "27" = "Cholangiocyte-like",
                                  "28" = "Mesenchymal cells",
                                  "29" = "Erythrocytes",
                                  "30" = "Basophils",
                                  "31" = "Mononucl. phagocytes")

# combined_seurat_filter  <- RenameIdents(combined_seurat_filter,
#                                 "0" = "Endothelial cells",
#                                 "1" = "Macrophages",
#                                  "2" = "Macrophages",
#                                  "3" = "Mesenchymal cells",
#                                  "4" = "Hepatocytes",
#                                  "5" = "Hepatocytes",
#                                  "6" = "Hepatocytes",
#                                  "7" = "Cholangiocyte-like",
#                                  "8" = "Endothelial cells",
#                                  "9" = "Monocytes",
#                                  "10" = "Neutrophils",
#                                  "11" = "T cells",
#                                  "12" = "Monocytes",
#                                  "13" = "Monocytes",
#                                  "14" = "B cells",
#                                  "15" = "Endothelial cells",
#                                  "16" = "Monocytes",
#                                  "17" = "T cells",
                                 
#                                  "18" = "B cells",
#                                  "19" = "Erythrocytes",
#                                   "20" = "Mesenchymal cells",
#                                   "21" = "NK cells",
#                                   "22" = "Monocytes",
#                                   "23" = "Monocytes",
#                                   "24" = "Macrophages",
#                                   "25" = "Macrophages",
#                                   "26" = "Monocytes",
#                                   "27" = "Cholangiocyte-like",
#                                   "28" = "Mesenchymal cells",
#                                   "29" = "Erythrocytes",
#                                   "30" = "Basophils",
#                                   "31" = "Macrophages")

combined_seurat_filter
table(Idents(combined_seurat_filter))
head(Idents(combined_seurat_filter))

RidgePlot(combined_seurat_filter,features = "Lcn2",split.by  = "group")
VlnPlot(combined_seurat_filter,features = "Spp1",split.by  = "group")
DotPlot(combined_seurat_filter,features = "Lcn2")
FeaturePlot(combined_seurat_filter,reduction = "ft_tsne.scvi",features = "Lcn2",split.by  = "group")
FeaturePlot(combined_seurat_filter,reduction = "ft_tsne.scvi",features = "Lcn2",split.by  = "group")

p1<- DimPlot(
  combined_seurat_filter,
  dims = c(1, 2),
  cells = NULL,
  cols = NULL,
  pt.size = NULL,
  reduction = "ft_tsne.scvi",
  split.by = NULL,
  shape.by = NULL,
  order = NULL,
  shuffle = FALSE,
  seed = 1,
  label = TRUE,
  label.size = 4,
  label.color = "black",
  label.box = FALSE,
  repel = TRUE,
  alpha = 0.2,
  cells.highlight = NULL,
  cols.highlight = "#DE2D26",
  sizes.highlight = 1,
  na.value = "grey50",
  ncol = 2,
  combine = TRUE,
  raster = NULL,
  raster.dpi = c(512, 512)
)+ NoLegend()
p1


colnames(combined_seurat_filter@meta.data)

cell_props <- prop.table(table(combined_seurat_filter@active.ident, combined_seurat_filter$group), margin = 2)
library(ggplot2)
cell_props_df <- as.data.frame(cell_props)

p2 <- ggplot(cell_props_df, aes(x = Var2, y = Freq, fill = Var1)) + 
  geom_bar(stat = "identity") +
  labs(x = "Group", y = "Proportion", fill = "Cell Type") +
  theme_classic()
library(patchwork)
p1 + p2 + plot_layout(widths = c(4,1))

chisq.test(table(combined_seurat_filter@active.ident, combined_seurat_filter$group))
save.dir  <- "SeuratObjects"

save(
  combined_seurat_filter,
  file = paste0(
    save.dir,"/",get_time(),
    "combined_seurat_filter_combinedN6_reselute_clusters_slm_annotated",
    ".RData"
  )
)
library("loupeR")
create_loupe_from_seurat(
  combined_seurat_filter,
  output_name = paste0(save.dir,"/",get_time(),
    "combined_seurat_filter_combinedN6_reselute_clusters_slm_annotated"
  )
)

