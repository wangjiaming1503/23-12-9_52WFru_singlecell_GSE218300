library(Seurat)
library(tidyverse)
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
