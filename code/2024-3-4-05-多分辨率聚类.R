
# 这是后面加进去的，不同的cluster算法如何避免记住参数，把其包装出来-------------------------------------------
Reselute_clusters_slm <- function(combined_seurat, resolution) {
  combined_seurat <-
    FindClusters(combined_seurat,
      resolution = resolution,
      # cluster.name = "sct_harmony_2_clusters",
      algorithm = 3
    )
  return(combined_seurat)
}

Reselute_clusters_louvain_muti <-
  function(combined_seurat, resolution) {
    combined_seurat <-
      FindClusters(combined_seurat,
        resolution = resolution,
        # cluster.name = "sct_harmony_2_clusters",
        algorithm = 2
      )
    return(combined_seurat)
  }

Reselute_clusters_leiden <- function(combined_seurat, resolution) {
  combined_seurat <-
    FindClusters(combined_seurat,
      resolution = resolution,
      # cluster.name = "sct_harmony_2_clusters",
      algorithm = 4
    )
  return(combined_seurat)
}

seq(0.5, 2, by = 0.1)
combined_seurat <-
  Reselute_clusters_slm(combined_seurat, seq(0.5, 2, by = 0.1))


# 通过gusb函数替换含有"SCT_snn_res"的列名，这里先检查再确定怎么换
colnames(combined_seurat@meta.data) <-
  gsub(
    "SCT_snn_res",
    "SCT_snn_slm_res",
    colnames(combined_seurat@meta.data)
  )

# 确认替换成功
colnames(combined_seurat@meta.data)
save(
  combined_seurat,
  file = paste0(
    save.dir,
    "/combined_seurat_combined_umap_N4--SCT-N3-cca-N2-louvain_slm",
    Sys.Date(),
    ".RData"
  )
)
combined_seurat <-
  Reselute_clusters_louvain_muti(combined_seurat, seq(0.5, 2, by = 0.1))

# 通过gusb函数替换含有"SCT_snn_res"的列名
colnames(combined_seurat@meta.data) <-
  gsub(
    "SCT_snn_res",
    "SCT_snn_muti_louvain_res",
    colnames(combined_seurat@meta.data)
  )

# 确认替换成功
colnames(combined_seurat@meta.data)


# 这是没考虑分辨率的时候
save(
  combined_seurat,
  file = paste0(
    save.dir,
    "/combined_seurat_combined_umap_N4--SCT-N3-cca-N3-louvain_slm_muti",
    Sys.Date(),
    ".RData"
  )
)

# 第一次试leiden，最高分辨率，崩了
combined_seurat <-
  FindClusters(combined_seurat,
    resolution = 2,
    # cluster.name = "sct_harmony_2_clusters",
    algorithm = 4
  )


save(
  combined_seurat,
  file = paste0(
    save.dir,
    "/combined_seurat_combined_umap_N4-N3-anatoedd-cluster-sct-harmony-ladden-diff-resolution-v1",
    Sys.Date(),
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

combined_seurat <- Reselute_clusters_leiden(combined_seurat, 0.5)


colnames(combined_seurat@meta.data)

# 保存了harmony的SCTumap，不过可以看看有啥参数，这个可能是设定了dims
combined_seurat <-
  RunUMAP(
    combined_seurat,
    reduction = "integrated.sct.harmony",
    dims = 1:30,
    reduction.name = "umap.2.sct.harmony"
  )

# 通过gusb函数替换含有"SCT_snn_res"的列名
colnames(combined_seurat@meta.data) <-
  gsub(
    "SCT_snn_res",
    "SCT_snn_Louvain_res",
    colnames(combined_seurat@meta.data)
  )

# 确认替换成功
colnames(combined_seurat@meta.data)


# 把不同分辨率的画树，宽度是按照50个细胞减半，高度是按照分辨率的数量，0.5-2是15个，
# 这样不会重叠，如果文件存在会报错，想一想有啥方法，也许把后缀改为时间好一点
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

colnames(combined_seurat@meta.data)

p_cluster_muti_louvain <-
  tree_cluster(combined_seurat,
    prefix = "SCT_snn_muti_louvain_res.",
    width = 25,
    height = 16
  )

p_cluster_slm <-
  tree_cluster(combined_seurat,
    prefix = "SCT_snn_slm_res.",
    width = 25,
    height = 16
  )

p_cluster_louvain <-
  tree_cluster(combined_seurat,
    prefix = "SCT_snn_Louvain_res.",
    width = 25,
    height = 16
  )

# 列出p_cluster开头的对象名，保存他们
p_cluster_object <- ls(pattern = "p_cluster")
p_cluster_object

# 建立文件夹
dir.create("plotobject")
# 保存到plotobject文件夹下，命名为cluster_muti_res.RData
save(
  list = ls(pattern = "p_cluster"),
  file = paste0(
    "plotobject/cluster_muti_louvain_res",
    Sys.Date(),
    ".RData"
  ),
  compress = TRUE
)

load(gzfile("plotobject/cluster_muti_louvain_res2024-03-02.RData.gz"))

