setwd("/home/rstudio/work/23-12-9_52WFru_singlecell_GSE218300")
.libPaths("/home/rstudio/R/x86_64-pc-linux-gnu-library/4.3_host")
.libPaths()
library(Seurat)
library(SeuratWrappers)
options(future.globals.maxSize = 1e9)
library(sceasy)
library(reticulate)
library(anndata)
library(ggplot2)
get_time <- function(my.tz = "Asia/Shanghai") {
  # 获取当前时间，并按照指定时区和格式化字符串
  datetime <- format(Sys.time(), "%Y-%m-%d-%H-%M", tz = my.tz)

  # 返回格式化的日期时间字符串
  return(datetime)
}
get_time()

save.dir <- "SeuratObjects"




load("SeuratObjects/combined_seurat_combined_umap_N4_harmony_cca_rcpa-mnn-keepsct-SCVI-intergraed-umap-add-learn-N3-slm-louvain-muti_joined2024-03-06-00-27.RData")
