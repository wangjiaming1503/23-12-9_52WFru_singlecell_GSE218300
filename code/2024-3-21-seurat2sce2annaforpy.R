library(Seurat)
library(tidyverse)
library(SingleCellExperiment)
combined_seurat <-
  readRDS(
    "~/work/23-12-9_52WFru_singlecell_GSE218300/SeuratObjects/2024-03-15-19-26_combined_seurat_anotated.rds"
  )

sce <- as.SingleCellExperiment(combined_seurat)

library(zellkonverter)

zellkonverter::writeH5AD(sce, "SeuratObjects/2024-03-21-00-15combined_seurat_anotated.h5ad", verbose = TRUE)
View(maartenutils::gen_file_overview("SeuratObjects"))

#python完成了过滤，进行双细胞标记
View(maartenutils::gen_file_overview("SeuratObjects"))
library(zellkonverter)
sce <- zellkonverter::readH5AD("./SeuratObjects/2024-03-21-00-25anotated-alltsne-combined_seurat-quality_control.h5ad", verbose = TRUE)

library(scDblFinder)
sce = scDblFinder(sce)
#Training model...
#iter=0, 5864 cells excluded from training.
#iter=1, 7735 cells excluded from training.
#iter=2, 7659 cells excluded from training.
#Threshold found:0.016
#6058 (17.6%) doublets called

zellkonverter::writeH5AD(sce, "SeuratObjects/2024-03-21-00-15combined_seurat_anotated_double.h5ad", verbose = TRUE)

#完成了部分normalization导入到R中
#2024-03-21-00-58anotated-alltsne-combined_seurat-quality_control_doublet_normalized.h5ad
sce <- zellkonverter::readH5AD("./SeuratObjects/2024-03-21-00-58anotated-alltsne-combined_seurat-quality_control_doublet_normalized.h5ad", verbose = TRUE)

library(scry)
sce <- devianceFeatureSelection(sce,assay="counts", sorted=TRUE)
plot(rowData(sce)$binomial_deviance, type="l", xlab="ranked genes",
     ylab="binomial deviance", main="Feature Selection with Deviance")
abline(v=2000, lty=2, col="red")
# 获取binomial deviance值
binomial_deviance <- rowData(sce)$binomial_deviance

# 对binomial_deviance向量进行排序,选择前4000个高度偏差的基因
idx <- order(binomial_deviance, decreasing = TRUE)[1:2000]
highly_deviant <- rep(FALSE, nrow(sce))
highly_deviant[idx] <- TRUE

# 将highly_deviant和binomial_deviance保存到adata的rowData中
rowData(sce)$highly_deviant <- highly_deviant  
rowData(sce)$binomial_deviance <- binomial_deviance

zellkonverter::writeH5AD(sce, "SeuratObjects/2024-03-21-01-09anotated-alltsne-combined_seurat-quality_control_doublet_normalized_feature_selected.h5ad", verbose = TRUE)

sce <- zellkonverter::readH5AD("./SeuratObjects/2024-03-21-01-23anotated-alltsne-combined_seurat-quality_control_doublet_normalized_featureselect_pca_tsne-umap-leiden.h5ad", verbose = TRUE)
sce
combined_seurat
summary(colnames(sce))
summary(colnames(combined_seurat))
metadata(sce)

#py to R v1 -----------------------------------------------------------------------------
#这里应该是自动补全了，试下不在anadata中subset。
#2024-03-21-11-06anotated-alltsne-combined_seurat-quality_control_notsubseted.h5ad
library(zellkonverter)
sce <- zellkonverter::readH5AD("./SeuratObjects/2024-03-21-11-06anotated-alltsne-combined_seurat-quality_control_notsubseted.h5ad", verbose = TRUE)

library(scDblFinder)
library(BiocParallel)

sce = scDblFinder(sce,samples="orig.ident", BPPARAM=MulticoreParam(3))
doublet_score = sce$scDblFinder.score
doublet_class = sce$scDblFinder.class
table(doublet_class)
combined_seurat[['doublet_score']] = doublet_score
combined_seurat[['doublet_class']] = doublet_class
combined_seurat@meta.data
str(combined_seurat,max.level = 2)
colnames(combined_seurat@meta.data)

str(sce,max.level = 3)

summary(combined_seurat@meta.data)
colnames(colData(sce))

colnames(colData(sce) [,183:203])

combined_seurat[['outlier']] = sce$outlier
combined_seurat[['mt_outlier']] = sce$mt_outlier
#combined_seurat[['doublet_score']] = sce$doublet_score
#combined_seurat[['doublet_class']] = sce$doublet_class
# > colnames(colData(sce) [,183:203])
#  [1] "ident"                      "n_genes_by_counts"         
#  [3] "log1p_n_genes_by_counts"    "total_counts"              
#  [5] "log1p_total_counts"         "pct_counts_in_top_20_genes"
#  [7] "total_counts_mt"            "log1p_total_counts_mt"     
#  [9] "pct_counts_mt"              "total_counts_ribo"         
# [11] "log1p_total_counts_ribo"    "pct_counts_ribo"           
# [13] "total_counts_hb"            "log1p_total_counts_hb"     
# [15] "pct_counts_hb"              "outlier"                   
# [17] "mt_outlier"                 "scDblFinder.class"         
# [19] "scDblFinder.score"          "scDblFinder.weighted"      
# [21] "scDblFinder.cxds_score"  

#把上面的都加到combined_seurat中
combined_seurat[['n_genes_by_counts']] = colData(sce)$n_genes_by_counts
combined_seurat[['log1p_n_genes_by_counts']] = colData(sce)$log1p_n_genes_by_counts
combined_seurat[['total_counts']] = colData(sce)$total_counts
combined_seurat[['log1p_total_counts']] = colData(sce)$log1p_total_counts
combined_seurat[['pct_counts_in_top_20_genes']] = colData(sce)$pct_counts_in_top_20_genes
combined_seurat[['total_counts_mt']] = colData(sce)$total_counts_mt
combined_seurat[['log1p_total_counts_mt']] = colData(sce)$log1p_total_counts_mt
combined_seurat[['pct_counts_mt']] = colData(sce)$pct_counts_mt
combined_seurat[['total_counts_ribo']] = colData(sce)$total_counts_ribo
combined_seurat[['log1p_total_counts_ribo']] = colData(sce)$log1p_total_counts_ribo
combined_seurat[['pct_counts_ribo']] = colData(sce)$pct_counts_ribo
combined_seurat[['total_counts_hb']] = colData(sce)$total_counts_hb
combined_seurat[['log1p_total_counts_hb']] = colData(sce)$log1p_total_counts_hb
combined_seurat[['pct_counts_hb']] = colData(sce)$pct_counts_hb
combined_seurat[['outlier']] = colData(sce)$outlier
combined_seurat[['mt_outlier']] = colData(sce)$mt_outlier
combined_seurat[['doublet_score']] = colData(sce)$scDblFinder.score
combined_seurat[['doublet_class']] = colData(sce)$scDblFinder.class
combined_seurat[['scDblFinder.weighted']] = colData(sce)$scDblFinder.weighted
combined_seurat[['scDblFinder.cxds_score']] = colData(sce)$scDblFinder.cxds_score

colnames(combined_seurat@meta.data)
save(combined_seurat, file = "SeuratObjects/2024-03-26-12-54anotated-alltsne-combined_seurat-quality_control_doublet_normalized_addmeta.RData")
load("SeuratObjects/2024-03-21-11-23anotated-alltsne-combined_seurat-quality_control_doublet_normalized_addmeta.RData")
library("loupeR")
create_loupe_from_seurat(
  combined_seurat,
  output_name = paste0(
    "SeuratObjects/2024-03-21-11-23anotated-alltsne-combined_seurat-quality_control_doublet_normalized_addmeta"
  )
)
# 将 outlier 转换为因子型
combined_seurat@meta.data$outlier <- factor(combined_seurat@meta.data$outlier, levels = c("TRUE", "FALSE"))

# 将 mt_outlier 转换为因子型 
combined_seurat@meta.data$mt_outlier <- factor(combined_seurat@meta.data$mt_outlier, levels = c("TRUE", "FALSE"))


#str(combined_seurat@meta.data)

create_loupe_from_seurat(
  combined_seurat,
  output_name = paste0(
    "SeuratObjects/2024-03-26-13-00anotated-alltsne-combined_seurat-quality_control_doublet_normalized_addmeta"
  )
)

summary(combined_seurat@meta.data)
summary(combined_seurat@meta.data$outlier)
combined_seurat_filter <-
  subset(combined_seurat, subset = outlier == "FALSE" &
    mt_outlier == "FALSE" & doublet_class == "singlet")
combined_seurat_filter

saveRDS(combined_seurat_filter, file = "SeuratObjects/2024-03-26-13-03anotated-alltsne-combined_seurat-quality_control_doublet_normalized_addmeta_filter.rds")
create_loupe_from_seurat(
  combined_seurat_filter,
  output_name = paste0(
    "SeuratObjects/2024-03-21-11-55anotated-alltsne-combined_seurat-quality_control_doublet_normalized_addmeta_filter"
  )
)
combined_seurat@meta.data$outlier
#把combined_seurat@meta.data$outlier里面的NULL改为TRUE，或者不是FALSE的都改为TRUE
# combined_seurat@meta.data$outlier[is.null(combined_seurat@meta.data$outlier) | combined_seurat@meta.data$outlier != FALSE] <- TRUE
# str(combined_seurat@meta.data)

combined_seurat_filter <- combined_seurat

combined_seurat  <- combined_seurat_filter

DefaultAssay(combined_seurat_filter) <- "SCT"

combined_seurat

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