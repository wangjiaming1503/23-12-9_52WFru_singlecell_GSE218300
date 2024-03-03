library(Seurat)
library(patchwork)
library(ggplot2)



load("~/work/23-12-9_52WFru_singlecell_GSE218300/SeuratObjects/seurat_list.RData")
seurat.list


# 假设你有一个包含多个Seurat对象的列表，名为seurat.list
# 其中每个对象代表不同处理条件下的一组细胞

# 合并Seurat对象
combined_seurat <- merge(seurat.list[[1]], y = seurat.list[-1])


# 把meta.data 中的orig.ident的_ _ 两个中间的部分提取,加入到group列
combined_seurat$group <- gsub(".*_([^_]+)_.*", "\\1", combined_seurat$orig.ident)


  
  combined_seurat[["percent.mt"]] <- PercentageFeatureSet(combined_seurat, pattern = "^MT-")
  
  
  
  # Visualize QC metrics as a violin plot
  P1 <- VlnPlot(combined_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  combined_seurat <- subset(combined_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  P2 <- VlnPlot(combined_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  P1 / P2
  # 加上时间保存图片，QC
  
  ggsave(paste0("violinplot.png",Sys.Date(),".png"))

  # 标准化数据
  combined_seurat <- NormalizeData(combined_seurat)
  
  # 查找变量特征
  combined_seurat <- FindVariableFeatures(combined_seurat)
  # 进行线性降维，比如PCA
  # combined_seurat <- ScaleData(combined_seurat)
  
  combined_seurat <- SCTransform(object = combined_seurat)
  
  combined_seurat <- RunPCA(combined_seurat)
  
  # 进行非线性降维，比如UMAP
  combined_seurat <- RunUMAP(combined_seurat, dims = 1:20)
  
  # 可视化UMAP结果
  p3<- DimPlot(combined_seurat, reduction = "umap", group.by = c("orig.ident","group"))
  
  p3
  
  combined_seurat <- IntegrateLayers(object = combined_seurat,
                                     method = CCAIntegration, 
                                     normalization.method = "SCT", 
                                     verbose = F
                            )
  
  
  combined_seurat <- FindNeighbors(combined_seurat, reduction = "integrated.dr", dims = 1:30)
  
  combined_seurat <- FindClusters(combined_seurat, resolution = 0.6)
  
  combined_seurat <- RunUMAP(combined_seurat, dims = 1:30, reduction = "integrated.dr")
  
  p4<- DimPlot(combined_seurat, reduction = "umap", group.by = c("orig.ident","group"))
  p4
  p3 / p4
  
  # 加上日期
  ggsave(paste0("umap_merged_intergrated_",Sys.Date(),".png"), p4, width = 10, height = 5)
  ggsave(paste0("umap_merged_intergrated_compair",Sys.Date(),".png"), p3 / p4, width = 10, height = 10)
  
  
  # 加上日期时间具体到小时分,文件夹在SeuratObjects下
  save(combined_seurat, file = paste0("combined_seurat_",Sys.Date(),".RData"))
  
  Idents(combined_seurat) <- "seurat_annotations"
  
  
  # nk.markers <- FindConservedMarkers(combined_seurat, ident.1 = "NK", grouping.var = "group", verbose = FALSE)
  # 
  # head(nk.markers)



