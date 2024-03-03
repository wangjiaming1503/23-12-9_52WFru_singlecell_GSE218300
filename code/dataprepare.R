
##============= 数据读取保存====================================
library(Seurat)

# 设置包含所有数据文件夹的父目录路径
parent.dir <- "~/work/23-12-9_52WFru_singlecell_GSE218300/data"

# 获取所有子目录名称
folders <- list.dirs(path = parent.dir, full.names = TRUE, recursive = FALSE)

# 初始化一个空列表来存储Seurat对象
seurat.list <- list()

# 循环遍历每个子目录并导入数据
for (folder in folders) {
  # 提取文件夹名称作为Seurat对象的名称
  folder.name <- basename(folder)
  
  # 检查路径是否存在
  if (dir.exists(folder)) {
    # 导入数据
    seurat.data <- Read10X(data.dir = folder)
    
    # 创建Seurat对象
    seurat.obj <- CreateSeuratObject(counts = seurat.data, project = folder.name)
    
    # 将Seurat对象添加到列表中，以文件夹名称作为键
    seurat.list[[folder.name]] <- seurat.obj
  } else {
    warning(paste("Data directory does not exist:", folder))
  }
}
# 删除不再需要的变量
rm(seurat.data, seurat.obj)
# 检查导入的Seurat对象列表
print(names(seurat.list))


# ================函数化不留变量在外==============================

library(Seurat)

# 定义一个函数来读取数据并创建Seurat对象
read_and_create_seurat <- function(folder) {
  # 导入数据
  seurat.data <- Read10X(data.dir = folder)
  
  # 创建Seurat对象
  seurat.obj <- CreateSeuratObject(counts = seurat.data, project = basename(folder))
  
  return(seurat.obj)
}

# 设置包含所有数据文件夹的父目录路径
parent.dir <- "~/work/2023-12-9_52WFru_singlecell_GSE218300/data"

# 获取所有子目录名称
folders <- list.dirs(path = parent.dir, full.names = TRUE, recursive = FALSE)

# 初始化一个空列表来存储Seurat对象
seurat.list <- list()

# 循环遍历每个子目录并导入数据
for (folder in folders) {
  # 检查路径是否存在
  if (dir.exists(folder)) {
    # 使用函数读取数据并创建Seurat对象
    seurat.list[[basename(folder)]] <- read_and_create_seurat(folder)
  } else {
    warning(paste("Data directory does not exist:", folder))
  }
}

# 检查导入的Seurat对象列表
print(names(seurat.list))

# 创建一个新的文件夹来存放Seurat对象
save.dir <- "SeuratObjects"
dir.create(save.dir)

# 分别保存列表中的每个Seurat对象
for (name in names(seurat.list)) {
  # 创建文件名，包含路径
  file.name <- file.path(save.dir, paste0("Seurat_", name, ".RData"))
  
  # 保存Seurat对象
  save(list = name, file = file.name, envir = list2env(seurat.list))
}

# 如果你也想保存整个列表，可以这样做
save.file <- file.path(save.dir, "seurat_list.RData")
save(seurat.list, file = save.file)

# 把seruat.list中的所有seurat对象放到环境中
list2env(seurat.list, envir = .GlobalEnv)


# 合并seurat对象

for (i in 1:length(seurat.list)) {
  seurat.list[[i]] <- NormalizeData(seurat.list[[i]], verbose = FALSE)
  seurat.list[[i]] <- FindVariableFeatures(
    seurat.list[[i]], selection.method = "vst",
    nfeatures = 2000, verbose = FALSE
  )
}


# 查找集成的锚点
anchors <- FindIntegrationAnchors(object.list = seurat.list, dims = 1:30)
saveRDS(anchors, file = "./SeuratObjects/anchors.rds")
# 整合数据
integrated.data <- IntegrateData(anchorset = anchors, dims = 1:30)
saveRDS(integrated.data, file = "./SeuratObjects/integrated_data.rds")



# 在整合数据上运行标准分析工作流
integrated.data <- ScaleData(integrated.data)
integrated.data <- RunPCA(integrated.data)
integrated.data <- FindNeighbors(integrated.data, dims = 1:30)
integrated.data <- FindClusters(integrated.data)
integrated.data <- RunUMAP(integrated.data, dims = 1:30)

saveRDS(integrated.data, file = "./SeuratObjects/integrated_data_afaz_v1.rds")

integrated.data@meta.data
# 把meta.data 中的orig.ident的_ _ 两个中间的部分提取,加入到group列
integrated.data$group <- gsub(".*_([^_]+)_.*", "\\1", integrated.data$orig.ident)


# 可视化
DimPlot(integrated.data, reduction = "umap", group.by = c("group","seurat_clusters"))
# 
# # ------------Azimuth示例--------------------# --------------------
# 
# library(Seurat)
# library(Azimuth)
# library(SeuratData)
# library(patchwork)
# # 
# # # integrated.data <- RunAzimuth(integrated.data, reference = "pbmcref")
# ifnb <- RunAzimuth(ifnb, reference = "pbmcref")
# 
# Error in `object[[assay]]`:
#   ! 'SCT' not found in this Seurat object
# 
# Run `rlang::last_trace()` to see where the error occurred.
# > rlang::last_trace()
# <error/rlang_error>
#   Error in `object[[assay]]`:
#   ! 'SCT' not found in this Seurat object
# 
# ---
#   Backtrace:
#   ▆
# 1. ├─Azimuth::RunAzimuth(ifnb, reference = "pbmcref")
# 2. └─Azimuth:::RunAzimuth.Seurat(ifnb, reference = "pbmcref")
# 3.   └─Seurat::PercentageFeatureSet(...)
# 4.     ├─features %||% ...
# 5.     ├─base::grep(...)
# 6.     │ └─base::is.factor(x)
# 7.     ├─base::rownames(x = object[[assay]][layer])
# 8.     ├─object[[assay]]
# 9.     └─SeuratObject:::`[[.Seurat`(object, assay)
# Run rlang::last_trace(drop = FALSE) to see 5 hidden frames.
# # ---
# # Error in `GetAssayData()`:
# #   ! GetAssayData doesn't work for multiple layers in v5 assay.
# # 
# # Backtrace:
# #     ▆
# #  1. ├─Azimuth::RunAzimuth(integrated.data, reference = "pbmcref")
# #  2. └─Azimuth:::RunAzimuth.Seurat(integrated.data, reference = "pbmcref")
# #  3.   └─Azimuth::ConvertGeneNames(...)
# #  4.     ├─SeuratObject::GetAssayData(object = object[["RNA"]], slot = "counts")
# #  5.     └─SeuratObject:::GetAssayData.StdAssay(object = object[["RNA"]], slot = "counts")
# #  6.       └─rlang::abort(...)
# 
# FeaturePlot(integrated.data, features = "CCR7")
# # <error/rlang_error>
# #   Error in `FeaturePlot()`:
# #   ! None of the requested features were found: CCR7 in slot  data
# 
# 

# 
# InstallData("pbmcsca")
# 
# # returns a Seurat object named pbmcsca
# pbmcsca <- LoadData("pbmcsca")
# 
# # The RunAzimuth function can take a Seurat object as input
# pbmcsca <- RunAzimuth(pbmcsca, reference = "pbmcref")
# 
# 
# p1 <- DimPlot(pbmcsca, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) + NoLegend()
# p2 <- DimPlot(pbmcsca, group.by = "Method")
# p1 + p2
# 
# pbmcsca <- NormalizeData(pbmcsca)
# Idents(pbmcsca) <- "predicted.celltype.l2"
# 
# p1 <- FeaturePlot(pbmcsca, features = "CCR7")
# p2 <- FeaturePlot(pbmcsca, features = "FCGR3A")
# p3 <- VlnPlot(pbmcsca, features = "AXL", group.by = "predicted.celltype.l2", idents = c("ASDC",
#                                                                                         "pDC", "cDC1", "cDC2"))
# p4 <- FeaturePlot(pbmcsca, features = "predictionscorecelltypel2_Treg")
# 
# p1 + p2 + p3 + p4 + plot_layout(ncol = 2)
# 
# # 查看reference数据,似乎不是很适合小鼠肝脏
# available_data <- AvailableData()
# available_data[grep("Azimuth", available_data[, 3]), 1:3]
# 
# # Dataset Version                        Summary
# # adiposeref.SeuratData         adiposeref   1.0.0     Azimuth Reference: adipose
# # bonemarrowref.SeuratData   bonemarrowref   1.0.0  Azimuth Reference: bonemarrow
# # fetusref.SeuratData             fetusref   1.0.0       Azimuth Reference: fetus
# # heartref.SeuratData             heartref   1.0.0       Azimuth Reference: heart
# # humancortexref.SeuratData humancortexref   1.0.0 Azimuth Reference: humancortex
# # kidneyref.SeuratData           kidneyref   1.0.1      Azimuth Reference: kidney
# # lungref.SeuratData               lungref   2.0.0        Azimuth Reference: lung
# # mousecortexref.SeuratData mousecortexref   1.0.0 Azimuth Reference: mousecortex
# # pancreasref.SeuratData       pancreasref   1.0.0    Azimuth Reference: pancreas
# # pbmcref.SeuratData               pbmcref   1.0.0        Azimuth Reference: pbmc
# # tonsilref.SeuratData           tonsilref   2.0.0      Azimuth Reference: tonsil
# 
# # --------------------
pbmc <- integrated.data
pbmc <- SCTransform(object = pbmc)
pbmc <- RunPCA(object = pbmc)
pbmc <- FindNeighbors(object = pbmc, dims = 1:30)
pbmc <- FindClusters(object = pbmc)
pbmc <- RunUMAP(object = pbmc, dims = 1:30)
DimPlot(object = pbmc, reduction = "umap")
DimPlot(pbmc, reduction = "umap", group.by = c("group","seurat_clusters"))
DimPlot(pbmc, reduction = "umap", group.by = c("orig.ident","seurat_clusters"))

# saveRDS(pbmc, file = "./SeuratObjects/integrated_data_afaz_SCT_v1.rds")
save(pbmc,file = "./SeuratObjects/integrated_data_afaz_SCT_v1.Rdata")

ifnb <- pbmc

ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])
ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$group)

# integrate datasets
ifnb <- IntegrateLayers(object = ifnb, method = CCAIntegration, normalization.method = "SCT", verbose = F)
# re-join layers after integration
ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])

# run standard anlaysis workflow
ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb)
ifnb <- ScaleData(ifnb)
ifnb <- RunPCA(ifnb) 

ifnb <- FindNeighbors(ifnb, dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 1)

ifnb <- RunUMAP(ifnb, dims = 1:30)

DimPlot(pbmc, reduction = "umap", group.by = c("group","seurat_clusters"))
DimPlot(ifnb, reduction = "umap", group.by = c("group","seurat_clusters"))
DimPlot(ifnb, reduction = "umap", group.by = c("orig.ident","seurat_clusters"))

save(ifnb,file = "./SeuratObjects/integrated_data_afaz_SCT_IntegrateLayers_v2.Rdata")
ifnb

FeaturePlot(ifnb, features = c("Adgre1", "Mrc1","Lcn2"), split.by = "group")


Layers(ifnb)
