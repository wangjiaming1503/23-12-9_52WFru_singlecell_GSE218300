load("~/work/23-12-9_52WFru_singlecell_GSE218300/SeuratObjects/combined_seurat_anotated_mouse_blu_immgen_hpa_v5_fine_tsne_joined_2024-02-24.RData")

library(Seurat)
library(reticulate)
library(sceasy)

sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)

ls()

combined_seurat[["RNA"]] <- JoinLayers(combined_seurat[["RNA"]])
# adata <- convertFormat(combined_seurat, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
# Error in ncol(df) : 
#   no slot of name "meta.features" for this object of class "Assay5"

print(adata) # Note generally in Python, dataset conventions are obs x var