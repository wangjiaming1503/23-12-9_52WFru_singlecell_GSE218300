setwd("/home/rstudio/work/23-12-9_52WFru_singlecell_GSE218300")
.libPaths()
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)
# 数据读入-------------------------------------------------------------
# Define a function to read data and create Seurat objects
read_and_create_seurat <- function(folder) {
  seurat.data <- Read10X(data.dir = folder)
  seurat.obj <-
    CreateSeuratObject(counts = seurat.data, project = basename(folder))
  return(seurat.obj)
}

# Set the parent directory path where all data folders are located
parent.dir <- "~/work/23-12-9_52WFru_singlecell_GSE218300/data"

# Get the names of all subdirectories
folders <-
  list.dirs(
    path = parent.dir,
    full.names = TRUE,
    recursive = FALSE
  )

# Initialize an empty list to store Seurat objects
seurat.list <- list()

# Loop through each subdirectory to import data
for (folder in folders) {
  if (dir.exists(folder)) {
    seurat.list[[basename(folder)]] <- read_and_create_seurat(folder)
  } else {
    warning(paste("Data directory does not exist:", folder))
  }
}

# Print the names of imported Seurat objects
print(names(seurat.list))

# Create a new folder to store Seurat objects
save.dir <- "SeuratObjects"
dir.create(save.dir)

# Save each Seurat object from the list to a separate file
for (name in names(seurat.list)) {
  file.name <- file.path(save.dir, paste0("Seurat_", name, ".RData"))
  save(
    list = name,
    file = file.name,
    envir = list2env(seurat.list)
  )
}

# Optionally, save the entire list of Seurat objects
save.file <- file.path(save.dir, "seurat_list.RData")
save(seurat.list, file = save.file)
# 完成数据读取----------数据过滤------------------===================================================
# Load the list of Seurat objects
save.dir <- "SeuratObjects"
load(file.path(save.dir, "seurat_list.RData"))

# Merge Seurat objects
combined_seurat <- merge(seurat.list[[1]], y = seurat.list[-1])

# Extract the middle part of the 'orig.ident' in meta.data and add it to a new column 'group'
combined_seurat$group <-
  gsub(".*_([^_]+)_.*", "\\1", combined_seurat$orig.ident)

# Calculate the percentage of mitochondrial genes
combined_seurat[["percent.mt"]] <-
  PercentageFeatureSet(combined_seurat, pattern = "^MT-")

# Visualize QC metrics
P1 <-
  VlnPlot(
    combined_seurat,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3
  )
combined_seurat <-
  subset(combined_seurat, subset = nFeature_RNA > 200 &
    nFeature_RNA < 2500 & percent.mt < 5)
P2 <-
  VlnPlot(
    combined_seurat,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3
  )
plot_violin <- P1 / P2
plot_violin
# Save violin plot with current date in the filename
ggsave(
  paste0("violinplot_", Sys.Date(), ".png"),
  plot_violin,
  width = 10,
  height = 10
)

save(
  combined_seurat,
  file = paste0(
    save.dir,
    "/combined_seurat_af_raw_merge_fillter_",
    Sys.Date(),
    ".RData"
  )
)