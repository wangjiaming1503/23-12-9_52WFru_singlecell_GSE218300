library(scry)

#安装包，用bioconductor
#BiocManager::install("scry")
library(zellkonverter)

# 加载已经归一化的数据集
sce <- zellkonverter::readH5AD("./SeuratObjects/2024-03-20-19-05anotated-alltsne-combined_seurat-quality_control_doublet_normalized.h5ad")
sce

# 使用deviance进行特征选择
sce <- devianceFeatureSelection(sce,assay="X")
str(sce,max.level=3)
names(sce)
col(sce)
# 获取binomial deviance值
binomial_deviance <- rowData(sce)$binomial_deviance

# 对binomial_deviance向量进行排序,选择前4000个高度偏差的基因
idx <- order(binomial_deviance, decreasing = TRUE)[1:4000]
highly_deviant <- rep(FALSE, nrow(adata))
highly_deviant[idx] <- TRUE

# 将highly_deviant和binomial_deviance保存到adata的rowData中
rowData(adata)$highly_deviant <- highly_deviant  
rowData(adata)$binomial_deviance <- binomial_deviance

# 保存结果
saveH5AD(adata, "s4d8_feature_selection.h5ad")