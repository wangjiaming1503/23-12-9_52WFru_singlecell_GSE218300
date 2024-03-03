library("loupeR")
create_loupe_from_seurat(
  combined_seurat,
  output_name = paste0(
    "SeuratObjects/combined_seurat_combined_umap_N4--SCT-N3-cca-N3-louvain_slm_muti_",
    format(Sys.Date(), "%Y-%m-%d")
  )
)
