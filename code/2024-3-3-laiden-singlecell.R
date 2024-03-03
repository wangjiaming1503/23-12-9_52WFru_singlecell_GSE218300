setwd("/home/rstudio/work/23-12-9_52WFru_singlecell_GSE218300")
.libPaths("/home/rstudio/R/x86_64-pc-linux-gnu-library/4.3")
.libPaths()
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
Sys.setenv(http_proxy = "http://172.27.255.100:7890")
Sys.setenv(https_proxy = "http://172.27.255.100:7890")
load(
  "SeuratObjects/combined_seurat_combined_umap_N4-N3-anatoedd-cluster-sct-harmony-ladden-diff-resolution-0-8_2_v3-2024-03-03.RData"
)

load(
  "SeuratObjects/combined_seurat_combined_umap_N4-N3-anatoedd-cluster-sct-harmony-ladden-diff-resolution-0-8v1-22024-03-02.RData"
)


options(future.globals.maxSize = 1e9)
library(reticulate)
Reselute_clusters_leiden <- function(combined_seurat, resolution) {
  combined_seurat <-
    FindClusters(combined_seurat,
                 resolution = resolution,
                 # cluster.name = "sct_harmony_2_clusters",
                 algorithm = 4)
  return(combined_seurat)
}
c(seq(0.5, 1.9, by = 0.1))[-4]
combined_seurat <-
  Reselute_clusters_leiden(combined_seurat, c(seq(0.5, 1.9, by = 0.1))[-4])

Reselute_clusters_slm <- function(combined_seurat, resolution) {
  combined_seurat <-
    FindClusters(combined_seurat,
                 resolution = resolution,
                 # cluster.name = "sct_harmony_2_clusters",
                 algorithm = 3)
  return(combined_seurat)
}


Reselute_clusters_louvain_muti <-
  function(combined_seurat, resolution) {
    combined_seurat <-
      FindClusters(combined_seurat,
                   resolution = resolution,
                   # cluster.name = "sct_harmony_2_clusters",
                   algorithm = 2)
    return(combined_seurat)
  }

Reselute_clusters_louvain <- function(combined_seurat, resolution) {
  combined_seurat <-
    FindClusters(combined_seurat,
                 resolution = resolution,
                 # cluster.name = "sct_harmony_2_clusters",
                 algorithm = 1)
  return(combined_seurat)
}

#0.1-0.5分辨率
resolution <- seq(0.1, 0.4, by = 0.1)
resolution

combined_seurat <-
  Reselute_clusters_slm(combined_seurat, resolution)

colnames(combined_seurat@meta.data)
str(combined_seurat@meta.data)
# # 通过gusb函数替换含有"SCT_snn_res"的列名
colnames(combined_seurat@meta.data) <-
  gsub("SCT_snn_res",
       "SCT_snn_harmony_slm_res",
       colnames(combined_seurat@meta.data))

colnames(combined_seurat@meta.data)


combined_seurat <-
  Reselute_clusters_louvain_muti(combined_seurat, resolution)

colnames(combined_seurat@meta.data)

colnames(combined_seurat@meta.data)[91:94] <-
  gsub(
    "SCT_snn_harmony_slm",
    "SCT_snn_harmony_muti_louvain",
    colnames(combined_seurat@meta.data[91:94])
  )
colnames(combined_seurat@meta.data)


# # 通过gusb函数替换含有"SCT_snn_res"的列名
colnames(combined_seurat@meta.data) <-
  gsub(
    "SCT_snn_res",
    "SCT_snn_harmony_muti_louvain_res",
    colnames(combined_seurat@meta.data)
  )

combined_seurat <-
  Reselute_clusters_louvain(combined_seurat, resolution)
colnames(combined_seurat@meta.data)

# # 通过gusb函数替换含有"SCT_snn_res"的列名
colnames(combined_seurat@meta.data) <-
  gsub(
    "SCT_snn_res",
    "SCT_snn_harmony_Louvain_res",
    colnames(combined_seurat@meta.data)
  )

combined_seurat <-
  Reselute_clusters_leiden(combined_seurat, resolution)

colnames(combined_seurat@meta.data)
# # 通过gusb函数替换含有"SCT_snn_res"的列名
colnames(combined_seurat@meta.data) <-
  gsub("SCT_snn_res",
       "SCT_snn_harmony_leiden_res",
       colnames(combined_seurat@meta.data))

save(
  combined_seurat,
  file = paste0(
    save.dir,
    "/combined_seurat_combined_umap_N4-N3-anatoedd-cluster-sct-harmony-ladden-diff-resolution-all_0_1-2-before-leiden-v1-",
    Sys.Date(),
    ".RData"
  )
)



# #前处理把全部加上harmony,去掉部分harmony----------------------
colnames(combined_seurat@meta.data) <-
  gsub("SCT_snn_harmony",
       "SCT_snn",
       colnames(combined_seurat@meta.data))

# 确认替换成功
colnames(combined_seurat@meta.data)

# #把全部加上harmony,去掉部分harmony-换回来
colnames(combined_seurat@meta.data) <-
  gsub("SCT_snn",
       "SCT_snn_harmony",
       colnames(combined_seurat@meta.data))

# 确认替换成功
colnames(combined_seurat@meta.data)


use_condaenv("/home/rstudio/work/anaconda3/envs/leiden")
combined_seurat <-
  FindClusters(combined_seurat,
               resolution = 2,
               # cluster.name = "sct_harmony_2_clusters",
               algorithm = 4)




save.dir <- "SeuratObjects"

# 通过gusb函数替换含有"SCT_snn_res"的列名
colnames(combined_seurat@meta.data) <-
  gsub("SCT_snn_harmony_leiden_res",
       "SCT_snn_res",
       colnames(combined_seurat@meta.data))

# 通过gusb函数替换含有"SCT_snn_res"的列名
colnames(combined_seurat@meta.data) <-
  gsub("SCT_snn_leiden",
       "SCT_snn_harmony_leiden",
       colnames(combined_seurat@meta.data))

# # 通过gusb函数替换含有"SCT_snn_res"的列名
colnames(combined_seurat@meta.data) <-
  gsub("SCT_snn_",
       "SCT_snn_harmony_",
       colnames(combined_seurat@meta.data))





save(
  combined_seurat,
  file = paste0(
    save.dir,
    "/combined_seurat_combined_umap_N4-N3-anatoedd-cluster-sct-harmony-ladden-diff-resolution-all_v4-",
    Sys.Date(),
    ".RData"
  )
)

load(
  paste0(
    save.dir,
    "/combined_seurat_combined_umap_N4-N3-anatoedd-cluster-sct-harmony-ladden-diff-resolution-all_v4-",
    Sys.Date(),
    ".RData"
  )
)

library("loupeR")
create_loupe_from_seurat(
  combined_seurat,
  output_name =  paste0(
    save.dir,
    "/combined_seurat_combined_umap_N4-N3-anatoedd-cluster-sct-harmony-ladden-diff-resolution-v10-8_2_v2-",
    format(Sys.Date(), "%Y-%m-%d")
  )
)

# relevel combined_seurat@meta.data$SCT_snn_leiden_res.0.8 按照数字大小重新排序
str(combined_seurat@meta.data)

str(combined_seurat@meta.data$SCT_snn_leiden_res.0.8)
# Factor w/ 34 levels "1","10","11",..: 25 1 1 1 1 8 1 8 12 5 ...
library(dplyr)
library(forcats)
library(Seurat)

# 假设 seurat_object 是你的 Seurat 对象
# 使用 AddMetaData 更新元数据，例如 SCT_snn_leiden_res.0.8

# 首先，提取你想要操作的元数据列
metadata <-
  FetchData(combined_seurat, vars = "SCT_snn_leiden_res.0.8")
combined_seurat$SCT_snn_leiden_res.0.8 %>% str()
# 接下来，按照需求修改这个元数据列
# 例如，根据数值大小对因子进行排序
metadata$SCT_snn_leiden_res.0.8 <-
  factor(metadata$SCT_snn_leiden_res.0.8,
         levels = sort(as.numeric(
           levels(metadata$SCT_snn_leiden_res.0.8)
         )))

# 最后，使用 AddMetaData 将修改后的元数据添加回 Seurat 对象
combined_seurat <- AddMetaData(combined_seurat, metadata = metadata)
# 检查更改后的因子结构
combined_seurat$SCT_snn_leiden_res.0.8 %>% str()
combined_seurat$SCT_snn_leiden_res.0.8

# List of resolution values to loop over
res_values <-
  c(0.5,
    0.6,
    0.7,
    0.8,
    0.9,
    1,
    1.1,
    1.2,
    1.3,
    1.4,
    1.5,
    1.6,
    1.7,
    1.8,
    1.9,
    2)

# Loop over each resolution value
for (res in res_values) {
  col_name <- paste0("SCT_snn_harmony_leiden_res.", res)
  
  # Check if the column exists in the metadata
  if (col_name %in% names(combined_seurat@meta.data)) {
    # Extract the factor and store as a character vector
    metadata_col <-
      FetchData(combined_seurat, vars = col_name)[, col_name]
    
    # Convert the factor levels to numeric, sort them, and convert back to a factor with the new levels
    combined_seurat <- Seurat::AddMetaData(
      combined_seurat,
      metadata = factor(metadata_col, levels = sort(as.numeric(
        levels(metadata_col)
      ))),
      col.name = col_name
    )
  }
}

# After running this loop, the factor levels for each resolution column should be numerically sorted.


save(
  combined_seurat,
  file = paste0(
    save.dir,
    "/combined_seurat_combined_umap_N4-N3-anatoedd-cluster-sct-harmony-ladden-diff-resolution-all_sorted-v6-",
    Sys.Date(),
    ".RData"
  )
)


library("loupeR")
create_loupe_from_seurat(
  combined_seurat,
  output_name =  paste0(
    save.dir,
    "/combined_seurat_combined_umap_N4-N3-anatoedd-cluster-sct-harmony-ladden-diff-resolution-all_0_1-2-after-leiden-v2",
    format(Sys.Date(), "%Y-%m-%d")
  )
)
colnames(combined_seurat@meta.data)

load(
  "~/work/23-12-9_52WFru_singlecell_GSE218300/SeuratObjects/combined_seurat_combined_umap_N4-N3-anatoedd-cluster-sct-harmony-ladden-diff-resolution-all_0_1-2-after-leiden-v2-2024-03-03.RData"
)

tree_cluster <-
  function(combined_seurat,
           prefix,
           width = 25,
           height = 16) {
    library(clustree)
    # install.packages("clustree")
    p_cluster <-
      clustree(combined_seurat, prefix = prefix)
    
    # # install.packages("ploty")
    # install.packages("plotly")
    #
    # library(plotly)
    #
    # ggplotly(p_cluster)
    
    # width = 25
    # height = 16
    ggsave(
      width = width,
      height = height,
      paste0(
        prefix,
        "width_",
        width,
        "_",
        "height",
        height,
        "_",
        Sys.Date(),
        ".pdf"
      )
      
    )
    return(p_cluster)
  }

p_cluster_leiden <-
  tree_cluster(combined_seurat,
               prefix = "SCT_snn_harmony_leiden_res.",
               width = 25,
               height = 20)


colnames(combined_seurat@meta.data)
# 更新聚类，增加了低分辨率数据
p_cluster_muti_louvain <-
  tree_cluster(combined_seurat,
               prefix = "SCT_snn_harmony_muti_louvain_res.",
               width = 25,
               height = 20)

p_cluster_slm <-
  tree_cluster(combined_seurat,
               prefix = "SCT_snn_harmony_slm_res.",
               width = 25,
               height = 20)

p_cluster_louvain <-
  tree_cluster(combined_seurat,
               prefix = "SCT_snn_harmony_Louvain_res.",
               width = 25,
               height = 20)

# 保存到plotobject文件夹下，命名为cluster_muti_res.RData
save(
  list = ls(pattern = "p_cluster"),
  file = paste0(
    "plotobject/cluster_slm_louvain_muti_leiden_0_1-2",
    Sys.Date(),
    ".RData"
  )
)


# ------------scvi-----------------------------------
setwd("/home/rstudio/work/23-12-9_52WFru_singlecell_GSE218300")
.libPaths("/home/rstudio/R/x86_64-pc-linux-gnu-library/4.3")
.libPaths()
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
Sys.setenv(http_proxy = "http://172.27.255.100:7890")
Sys.setenv(https_proxy = "http://172.27.255.100:7890")
load(
  "SeuratObjects/combined_seurat_combined_umap_N4-N3-anatoedd-cluster-sct-harmony-ladden-diff-resolution-0-8_2_v3-2024-03-03.RData"
)
library(reticulate)


use_condaenv("/home/rstudio/work/root_conda/anaconda3/envs/scvi-env")
py_config()

combined_seurat <- IntegrateLayers(
  object = combined_seurat,
  method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "/home/rstudio/work/root_conda/anaconda3/envs/scvi-env",
  verbose = FALSE
)
save.dir <- "SeuratObjects"