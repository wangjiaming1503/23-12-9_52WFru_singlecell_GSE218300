library(Seurat)
library(tidyverse)
load("~/work/23-12-9_52WFru_singlecell_GSE218300/SeuratObjects/Seurat_GSM6738781_WD_4538.RData")
load("~/work/23-12-9_52WFru_singlecell_GSE218300/SeuratObjects/Seurat_GSM6738785_WD_6929.RData")
load("~/work/23-12-9_52WFru_singlecell_GSE218300/SeuratObjects/Seurat_GSM6738786_WD_6930.RData")
load("~/work/23-12-9_52WFru_singlecell_GSE218300/SeuratObjects/Seurat_GSM6738791_Chow_6927.RData")
load("~/work/23-12-9_52WFru_singlecell_GSE218300/SeuratObjects/Seurat_GSM6738792_Chow_6928.RData")
ls()

pbmc <- GSM6738781_WD_4538
GSM6738781_WD_4538 <- ifnb
rm(pbmc)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
GSM6738781_WD_4538[["percent.mt"]] <- PercentageFeatureSet(GSM6738781_WD_4538, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(GSM6738781_WD_4538, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(GSM6738781_WD_4538, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(GSM6738781_WD_4538, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# 这个wd 4538的MT已经去掉过了

GSM6738781_WD_4538 <- subset(GSM6738781_WD_4538, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

GSM6738781_WD_4538 <- NormalizeData(GSM6738781_WD_4538, normalization.method = "LogNormalize", scale.factor = 10000)
GSM6738781_WD_4538 <- NormalizeData(GSM6738781_WD_4538)

GSM6738781_WD_4538 <- FindVariableFeatures(GSM6738781_WD_4538, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(GSM6738781_WD_4538), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(GSM6738781_WD_4538)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


GSM6738781_WD_4538 <- ScaleData(GSM6738781_WD_4538, features = rownames(GSM6738781_WD_4538))

GSM6738781_WD_4538 <- RunPCA(GSM6738781_WD_4538, features = VariableFeatures(object = GSM6738781_WD_4538))

# Examine and visualize PCA results a few different ways
print(GSM6738781_WD_4538[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(GSM6738781_WD_4538, dims = 1:2, reduction = "pca")
DimPlot(GSM6738781_WD_4538, reduction = "pca") + NoLegend()
DimHeatmap(GSM6738781_WD_4538, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(GSM6738781_WD_4538, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(GSM6738781_WD_4538)
DimHeatmap(GSM6738781_WD_4538, dims = 1:19, cells = 500, balanced = TRUE)
# 建议用户在选择该参数时宁可偏高 鼓励用户使用不同数量的 PC（10、15，甚至 50！）重复下游分析
GSM6738781_WD_4538 <- FindNeighbors(GSM6738781_WD_4538, dims = 1:15)
GSM6738781_WD_4538 <- FindClusters(GSM6738781_WD_4538, resolution = 0.5)

# 对于较大的数据集，最佳分辨率通常会增加。可以使用该Idents()函数找到簇。
head(Idents(GSM6738781_WD_4538), 5)


GSM6738781_WD_4538 <- RunUMAP(GSM6738781_WD_4538, dims = 1:15)

DimPlot(GSM6738781_WD_4538, reduction = "umap")
DimPlot(GSM6738781_WD_4538, reduction = "umap", label = TRUE)

save(GSM6738781_WD_4538, file = "~/work/23-12-9_52WFru_singlecell_GSE218300/SeuratObjects/Seurat_GSM6738781_WD_4538.RData")
# 
devtools::install_github("immunogenomics/presto")
# # 本地zip安装
# devtools::install_local("presto-master.zip")

# find all markers of cluster 2
cluster2.markers <- FindMarkers(GSM6738781_WD_4538, ident.1 = 2)

head(GSM6738781_WD_4538, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(GSM6738781_WD_4538, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(GSM6738781_WD_4538, only.pos = TRUE)

pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# ROC 测试返回任何单个标记的“分类能力”（范围从 0 - 随机，到 1 - 完美）
cluster0.markers <- FindMarkers(GSM6738781_WD_4538, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(GSM6738781_WD_4538, features = c("Igfbp7", "Lcn2"))

# you can plot raw counts as well
VlnPlot(GSM6738781_WD_4538, features = c("Igfbp7", "Gpihbp1"), slot = "counts", log = TRUE)
# PC_ 1 
# Positive:  Fabp4, Sparc, Plpp1, Bgn, Ctla2a 
# Negative:  Cd52, Lyz2, Ctsc, Csf1r, Ctsh 
# PC_ 2 
# Positive:  Fabp4, Tmsb4x, Plpp1, Cyp4b1, Klf2 
# Negative:  Krt18, Ttr, Alb, Cp, Serpina1c 
# PC_ 3 
# Positive:  Fxyd3, Tspan8, Krt19, Mmp7, Cldn7 
# Negative:  Apoc3, Uox, Serpina1a, Hpd, Adh1 
# PC_ 4 
# Positive:  Il18bp, Ctsb, Cd5l, Adgre1, Igf1 
# Negative:  Alox5ap, Pglyrp1, Coro1a, Retnlg, Gsr 
# PC_ 5 
# Positive:  Cldn7, Mmp7, Dsg2, Sox9, Wfdc15b 
# Negative:  Serping1, C1s1, Efemp1, Igfbp6, Prss23 
FeaturePlot(GSM6738781_WD_4538, features = c("Igfbp7", "Gpihbp1", "Lcn2","Il18bp", "Ctsb", "Cd5l"))

pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(GSM6738781_WD_4538, features = top10$gene) + NoLegend()

# # 以下为手动标注
new.cluster.ids <- c("Naive CD4 T-test", "CD14+ Mono-test", "Memory CD4 T-test", "B-test", "CD8 T-test", "FCGR3A+ Mono-test",
                     "NK-test", "DC-test", "Platelet","Microphage M1-test","Microphage M2-test")

# 数量暂时对不上
names(new.cluster.ids) <- levels(GSM6738781_WD_4538)
GSM6738781_WD_4538 <- RenameIdents(GSM6738781_WD_4538, new.cluster.ids)
DimPlot(GSM6738781_WD_4538, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

library(ggplot2)
plot <- DimPlot(GSM6738781_WD_4538, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))

ggsave(filename = "pbmc3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)

plot

# -----------------------------探索可视化方法----------------

# SeuratData::InstallData("ifnb")

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
# pbmc3k.final <- LoadData("pbmc3k", type = "pbmc3k.final")
# 此处为随机给予分组，创造了groups槽，后续合并后可以根据名字来分组
GSM6738781_WD_4538$groups <- sample(c("group1", "group2"), size = ncol(GSM6738781_WD_4538), replace = TRUE)


features <- c("Igfbp7", "Gpihbp1", "Lcn2","Il18bp", "Ctsb", "Cd5l")

GSM6738781_WD_4538

# Ridge plots - from ggridges. Visualize single cell expression distributions in each cluster
RidgePlot(GSM6738781_WD_4538, features = features, ncol = 2)

# Violin plot - Visualize single cell expression distributions in each cluster
VlnPlot(GSM6738781_WD_4538, features = features)

# Feature plot - visualize feature expression in low-dimensional space
FeaturePlot(GSM6738781_WD_4538, features = features)

# Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
DotPlot(GSM6738781_WD_4538, features = features) + RotatedAxis()

# Single cell heatmap of feature expression
DoHeatmap(subset(GSM6738781_WD_4538, downsample = 100), features = features, size = 3)

# Plot a legend to map colors to expression levels
FeaturePlot(GSM6738781_WD_4538, features = "Adgre1")
FeaturePlot(GSM6738781_WD_4538, features = "Mrc1")

# Adjust the contrast in the plot
FeaturePlot(GSM6738781_WD_4538, features = "Adgre1", min.cutoff = 1, max.cutoff = 3)

# Calculate feature-specific contrast levels based on quantiles of non-zero expression.
# Particularly useful when plotting multiple markers
FeaturePlot(GSM6738781_WD_4538, features = c("Adgre1", "Mrc1"), min.cutoff = "q10", max.cutoff = "q90")

# Visualize co-expression of two features simultaneously
FeaturePlot(GSM6738781_WD_4538, features = c("Adgre1", "Mrc1"), blend = TRUE)

# Split visualization to view expression by groups (replaces FeatureHeatmap)
FeaturePlot(GSM6738781_WD_4538, features = c("Adgre1", "Mrc1","Lcn2"), split.by = "group")

# ------------Updated and expanded visualization functions----------
# Violin plots can also be split on some variable. Simply add the splitting variable to object
# metadata and pass it to the split.by argument
VlnPlot(GSM6738781_WD_4538, features = c("Adgre1", "Mrc1"), split.by = "groups")

# SplitDotPlotGG has been replaced with the `split.by` parameter for DotPlot
DotPlot(GSM6738781_WD_4538, features = features, split.by = "groups") + RotatedAxis()

# DimPlot replaces TSNEPlot, PCAPlot, etc. In addition, it will plot either 'umap', 'tsne', or
# 'pca' by default, in that order
DimPlot(GSM6738781_WD_4538)

GSM6738781_WD_4538.no.umap <- GSM6738781_WD_4538
GSM6738781_WD_4538.no.umap[["umap"]] <- NULL
DimPlot(GSM6738781_WD_4538.no.umap) + RotatedAxis()

# DoHeatmap now shows a grouping bar, splitting the heatmap into groups or clusters. This can
# be changed with the `group.by` parameter
DoHeatmap(GSM6738781_WD_4538, features = VariableFeatures(GSM6738781_WD_4538)[1:100], cells = 1:500, size = 4,
          angle = 90) + NoLegend()

baseplot <- DimPlot(GSM6738781_WD_4538, reduction = "umap")
# Add custom labels and titles
baseplot + labs(title = "Clustering of 2,700 PBMCs")

# Use community-created themes, overwriting the default Seurat-applied theme Install ggmin
Sys.setenv(http_proxy = "http://192.168.0.105:7890")
Sys.setenv(https_proxy = "http://192.168.0.105:7890")

# remotes::install_github('sjessa/ggmin')

baseplot + ggmin::theme_powerpoint()

# Seurat also provides several built-in themes, such as DarkTheme; for more details see
# ?SeuratTheme
baseplot + DarkTheme()

# Chain themes together
baseplot + FontSize(x.title = 20, y.title = 20) + NoLegend()

# Include additional data to display alongside cell names by passing in a data frame of
# information.  Works well when using FetchData
plot <- FeaturePlot(GSM6738781_WD_4538, features = "Adgre1")
HoverLocator(plot = plot, information = FetchData(GSM6738781_WD_4538, vars = c("ident", "PC_1", "nFeature_RNA")))

# 手动选择 2024年2月21日11:50:04todo
# GSM6738781_WD_4538 <- RenameIdents(GSM6738781_WD_4538, DC = "CD14+ Mono")
plot <- DimPlot(GSM6738781_WD_4538, reduction = "umap")
select.cells <- CellSelector(plot = plot)

head(select.cells)

Idents(GSM6738781_WD_4538, cells = select.cells) <- "NewCells"

# Now, we find markers that are specific to the new cells, and find clear DC markers
newcells.markers <- FindMarkers(GSM6738781_WD_4538, ident.1 = "NewCells", min.diff.pct = 0.3,
                                only.pos = TRUE)
head(newcells.markers)

# -----组合图
plot1 <- DimPlot(GSM6738781_WD_4538)
# Create scatter plot with the Pearson correlation value as the title
plot2 <- FeatureScatter(GSM6738781_WD_4538, feature1 = "Adgre1", feature2 = "Mrc1")
# Combine two plots
plot1 + plot2

# Remove the legend from all plots
(plot1 + plot2) & NoLegend()

# Seurat Command List----------------------
# Standard Seurat workflow---------------
pbmc <- GSM6738781_WD_4538
pbmc <- NormalizeData(object = pbmc)
pbmc <- FindVariableFeatures(object = pbmc)
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(object = pbmc)
pbmc <- FindNeighbors(object = pbmc, dims = 1:30)
pbmc <- FindClusters(object = pbmc)
pbmc <- RunUMAP(object = pbmc, dims = 1:30)
p1 <- DimPlot(object = pbmc, reduction = "umap")

# SCtransform version---------------
pbmc <- SCTransform(object = pbmc)
pbmc <- RunPCA(object = pbmc)
pbmc <- FindNeighbors(object = pbmc, dims = 1:30)
pbmc <- FindClusters(object = pbmc)
pbmc <- RunUMAP(object = pbmc, dims = 1:30)
p2 <-DimPlot(object = pbmc, reduction = "umap")
p1 + p2
# note that you can chain multiple commands together with %>%
pbmc <- SCTransform(pbmc) %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:30)

# Seurat Object Data Access----------------------------
# Cell, feature, and layer names--------------------

# Get cell and feature names, and total numbers We show multiple ways to get the same output cell names
colnames(pbmc)
Cells(pbmc)

# feature names
Features(pbmc)
rownames(pbmc)

# number of cells/features
num_cells <- ncol(pbmc)
num_features <- nrow(pbmc)

# List of object layers
Layers(pbmc)

# # working with multimodal objects list assays
# Assays(cbmc)
# 
# # Assay-specific features (genes/ADT)
# Features(cbmc[["RNA"]])
# Features(cbmc[["ADT"]])

# Variable feature names
VariableFeatures(pbmc)

# Identity class labels

# ---------------Setting and retrieving cell identities------------------------------------------------

# Set identity classes to an existing column in meta data
Idents(object = pbmc) <- "seurat_annotations"

# View cell identities, get summary table
Idents(pbmc)
table(Idents(pbmc))

# Set identity to CD4 T cells for all cells
Idents(pbmc) <- "CD4 T cells"

# Set for a selected group of cells
pbmc.cells <- Cells(pbmc)
Idents(object = pbmc, cells = pbmc.cells[1:10]) <- "CD4 T cells"

# Get cell identity classes
Idents(object = pbmc)
levels(x = pbmc)

# Stash cell identity classes in metadata
pbmc[["old.ident"]] <- Idents(object = pbmc)
pbmc <- StashIdent(object = pbmc, save.name = "old.ident")

# Rename identity classes
pbmc <- RenameIdents(object = pbmc, `Naive CD4 T-test` = "T Helper cells")

# ------------------------Cell metadata----------------

# View metadata data frame, stored in object@meta.data
pbmc[[]]

# Retrieve specific values from the metadata
pbmc$nCount_RNA
pbmc[[c("seurat_clusters", "nFeature_RNA")]]

# Add metadata, see ?AddMetaData
random_group_labels <- sample(x = c("g1", "g2"), size = ncol(x = pbmc), replace = TRUE)
pbmc$groups <- random_group_labels

#---------------- Expression data (stored as layers in Seurat v5)-----------------------
str(pbmc)
# Retrieve data in an expression matrix RNA counts matrix
pbmc[["RNA"]]$counts

# Alternate accessor function with the same result
LayerData(pbmc, assay = "RNA", layer = "counts")

# GetAssayData from Seurat v4 is still supported
GetAssayData(object = pbmc, assay = "RNA", slot = "counts")

# ADT counts matrix (multimodal object)
cbmc[["ADT"]]$counts

# Set expression data assume new.data is a new expression matrix
pbmc[["RNA"]]$counts <- new.data

# Alternate setter function with the same result
LayerData(pbmc, assay = "RNA", layer = "counts") <- new.data

# SetAssayData from Seurat v4 is still supported
pbmc <- SetAssayData(object = pbmc, slot = "counts", new.data = new.data)

# Get cell embeddings and feature loadings stored on pbmc[['pca']]@cell.embeddings
Embeddings(pbmc, reduction = "pca")

# stored in pbmc[['pca]]@feature.loadings
Loadings(pbmc, reduction = "pca")

# Create custom dimensional reduction loadings matrix is optional

# new_reduction <- CreateDimReducObject(embeddings = new.embeddings, loadings = new.loadings, key = "custom_pca")
# pbmc[["custom_pca"]] <- new_reduction

# FetchData can access anything from expression matrices, cell embeddings, or metadata Use the previously listed
# commands to access entire matrices Use FetchData to access individual/small groups of variables
FetchData(object = pbmc, vars = c("PC_1", "nFeature_RNA", "MS4A1"), layer = "counts")

# Subsetting and merging

# Subset Seurat objects

# Get cell identity classes
Idents(object = pbmc)
levels(x = pbmc)


# Subset Seurat object based on identity class, also see ?SubsetData
subset(x = pbmc, idents = "B-test")
subset(x = pbmc, idents = c("Platelet", "CD8 T-test"), invert = TRUE)

# Subset on the expression level of a gene/feature
subset(x = pbmc, subset = Lcn2 > 2.5)

# Subset on a combination of criteria
subset(x = pbmc, subset = Lcn2 > 2.5 & PC_1 > 5)
subset(x = pbmc, subset = Lcn2 > 2.5, idents = "B-test")

# Subset on a value in the object meta data
subset(x = pbmc, subset = groups == "g1")

# Downsample the number of cells per identity class
subset(x = pbmc, downsample = 100)

# Split layers--------------------------
# In Seurat v5, users can now split in object directly into different layers keeps expression data in one object, but
# splits multiple samples into layers can proceed directly to integration workflow after splitting layers
InstallData("ifnb")
rm(ifnb)
ifnb
ifnb <- UpdateSeuratObject(object = ifnb)
ifnb
ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)
Layers(ifnb)
ifnb@assays
# If desired, for example after intergation, the layers can be joined together again
ifnb <- JoinLayers(ifnb)
Layers(ifnb)

# In line with prior workflows, you can also into split your object into a list of multiple objects based on a metadata
# column creates a list of two objects
ifnb_list <- SplitObject(ifnb, split.by = "stim")
ifnb_list$CTRL
ifnb_list$STIM

# Merge Seurat objects--------------------------

# Merge two Seurat objects
# merged_obj <- merge(x = ifnb_list$CTRL, y = ifnb_list$STIM)
merged_obj[["RNA"]] <- JoinLayers(merged_obj)

# Example to merge more than two Seurat objects
merge(x = pbmc1, y = list(pbmc2, pbmc3))


merged_obj <- merge(x = ifnb_list$CTRL, y = ifnb_list$STIM)
merged_obj <- NormalizeData(merged_obj)
merged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj)
merged_obj <- IntegrateLayers(object = merged_obj, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca",
                              verbose = FALSE)

# # now that integration is complete, rejoin layers
# merged_obj[["RNA"]] <- JoinLayers(merged_obj)
# 
# ? JoinLayers

# Pseudobulk analysis
# Group cells together, based on multiple categories

# pseudobulk cells only by cell type
bulk <- AggregateExpression(ifnb, group.by = "seurat_annotations", return.seurat = TRUE)
Cells(bulk)

# pseudobulk cells by stimulation condition AND cell type
bulk <- AggregateExpression(ifnb, group.by = c("stim", "seurat_annotations"), return.seurat = TRUE)
Cells(bulk)

# pseudobulk cells by stimulation condition AND cell type AND donor
bulk <- AggregateExpression(ifnb, group.by = c("stim", "seurat_annotations", "donor_id"), return.seurat = TRUE)
Cells(bulk)

# Visualization in Seurat----------------------------
# Dimensional reduction plot
DimPlot(object = pbmc, reduction = "pca")

# Dimensional reduction plot, with cells colored by a quantitative feature Defaults to UMAP if available
FeaturePlot(object = pbmc, features = "Lcn2")

# Scatter plot across single cells
FeatureScatter(object = pbmc, feature1 = "Lcn2", feature2 = "PC_1")

# Scatter plot across individual features, repleaces CellPlot
# CellScatter(object = pbmc, cell1 = "AGTCTACTAGGGTG", cell2 = "CACAGATGGTTTCT")

VariableFeaturePlot(object = pbmc)

# Violin and Ridge plots
VlnPlot(object = pbmc, features = c("Lcn2"))
RidgePlot(object = pbmc, feature = c("Lcn2"))


# Heatmaps (visualize scale.data slot)
DimHeatmap(object = pbmc, reduction = "pca", cells = 200)

# standard workflow
heatmap_markers <- c("Lcn2","Il6")
# rm(pbmc)
# pbmc3k.final <- LoadData("pbmc3k", type = "pbmc3k.final")
pbmc <- ScaleData(pbmc, features = c("Lcn2","Il6"))
DoHeatmap(object = pbmc, features = heatmap_markers)

# sctransform workflow
pbmc <- GetResidual(pbmc, features = heatmap_markers)
DoHeatmap(object = pbmc, features = heatmap_markers)

# heatmap with maximum of 100 cells per group
# DoHeatmap(pbmc, heatmap_markers, cells = subset(pbmc, downsample = 100))

