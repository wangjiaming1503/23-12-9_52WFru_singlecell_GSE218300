
# --------------------------自动注释----------------------------
load(
  "~/work/23-12-9_52WFru_singlecell_GSE218300/SeuratObjects/combined_seurat_2024-02-24.RData"
)
library(Seurat)
library(celldex)
ref_blu <- BlueprintEncodeData()
ref_m <- MouseRNAseqData()
ref_immgen <- ImmGenData()

# # 用SingleR进行细胞注释


combined_seurat[["RNA"]] <- JoinLayers(combined_seurat[["RNA"]])
save(combined_seurat, file = "SeuratObjects/combined_seurat_tsne_joined_2024-02-24.RData")
sce_combined_seurat <- as.SingleCellExperiment(combined_seurat)
save(sce_combined_seurat, file = "SeuratObjects/sce_combined_seurat_2024-02-24.RData")
library(SingleR)


pred_m <-
  SingleR(
    test = sce_combined_seurat,
    ref = ref_m,
    labels = ref_m$label.main
  )
table(pred_m$labels)
plotScoreHeatmap(pred_m)
combined_seurat[["MouseRNA_SingleR.labels"]] <- pred_m$labels
save(combined_seurat, file = "SeuratObjects/combined_seurat_anotated_mouse_v1_tsne_joined_2024-02-24.RData")
DimPlot(
  combined_seurat,
  group.by = "MouseRNA_SingleR.labels",
  reduction = "umap",
  label = TRUE
)

pred_blu <-
  SingleR(
    test = sce_combined_seurat,
    ref = ref_blu,
    labels = ref_blu$label.main
  )
table(pred_blu$labels)
plotScoreHeatmap(pred_blu)
combined_seurat[["Blueprint_SingleR.labels"]] <- pred_blu$labels
DimPlot(
  combined_seurat,
  group.by = "Blueprint_SingleR.labels",
  reduction = "umap",
  label = TRUE
)
save(combined_seurat, file = "SeuratObjects/combined_seurat_anotated_mouse_blu_v2_tsne_joined_2024-02-24.RData")

pred_immgen <-
  SingleR(
    test = sce_combined_seurat,
    ref = ref_immgen,
    labels = ref_immgen$label.main
  )
table(pred_immgen$labels)
plotScoreHeatmap(pred_immgen)
combined_seurat[["ImmGen_SingleR.labels"]] <- pred_immgen$labels
DimPlot(
  combined_seurat,
  group.by = "ImmGen_SingleR.labels",
  reduction = "umap",
  label = TRUE
)
save(combined_seurat, file = "SeuratObjects/combined_seurat_anotated_mouse_blu_immgen_v3_tsne_joined_2024-02-24.RData")

ref_hpa <- celldex::HumanPrimaryCellAtlasDat()
pred_hpa <-
  SingleR(
    test = sce_combined_seurat,
    ref = ref_hpa,
    labels = ref_hpa$label.main
  )
table(pred_hpa$labels)
plotScoreHeatmap(pred_hpa)
combined_seurat[["HPA_SingleR.labels"]] <- pred_hpa$labels
DimPlot(
  combined_seurat,
  group.by = "SingleR.labels_HPA",
  reduction = "umap",
  label = TRUE
)
save(combined_seurat, file = "SeuratObjects/combined_seurat_anotated_mouse_blu_immgen_hpa_v4_tsne_joined_2024-02-24.RData")
library("loupeR")
create_loupe_from_seurat(combined_seurat, output_name = "combined_seurat_combined_umap_N4--SCT-N3-cca-N3-louvain_slm_muti")

pred_immgen_fine <-
  SingleR(
    test = sce_combined_seurat,
    ref = ref_immgen,
    labels = ref_immgen$label.fine
  )
table(pred_immgen_fine$labels)
plotScoreHeatmap(pred_immgen_fine)
combined_seurat[["ImmGen_fine_SingleR.labels"]] <-
  pred_immgen_fine$labels
DimPlot(combined_seurat, group.by = "SingleR.labels_ImmGen_fine", reduction = "umap") + NoLegend()


# save(combined_seurat, file = "SeuratObjects/combined_seurat_anotated_mouse_blu_immgen_hpa_v5_fine_tsne_joined_2024-02-24.RData")
save(
  combined_seurat,
  file = paste0(
    "SeuratObjects/combined_seurat_combined_umap_N4--SCT-N3-cca-N3-louvain_slm_muti",
    format(Sys.Date(), "%Y-%m-%d"),
    ".RData"
  )
)

library("loupeR")
create_loupe_from_seurat(
  combined_seurat,
  output_name = paste0(
    "SeuratObjects/combined_seurat_combined_umap_N4--SCT-N3-cca-N3-louvain_slm_muti_",
    format(Sys.Date(), "%Y-%m-%d")
  )
)


# 首先，确保annotation文件夹存在，如果不存在则创建
if (!dir.exists("annotation")) {
  dir.create("annotation")
}

# 列出工作环境中所有以'pred'开头的对象的名称
pred_objects <- ls(pattern = "^pred")

# 使用save()函数保存这些对象
# 假设我们想要将文件保存为'pred_objects.RData'
save(list = pred_objects, file = "annotation/pred_objects.RData")

ref_objects <- ls(pattern = "^ref")

# 使用save()函数保存这些对象
# 假设我们想要将文件保存为'ref_objects.RData'
save(list = ref_objects, file = "annotation/ref_objects.RData")
