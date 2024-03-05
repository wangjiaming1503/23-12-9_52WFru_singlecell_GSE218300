
library(httpgd)
DimPlot(combined_seurat,
  group.by = c("SCT_snn_harmony_slm_res.2"),
  label = TRUE
)
colnames(combined_seurat@meta.data)
dev.off()
combined_seurat@commands
# str(combined_seurat)
# Prepare the Seurat object for differential expression analysis
combined_seurat <- PrepSCTFindMarkers(combined_seurat)
Idents(combined_seurat)

# After the preparation, you can try running FindAllMarkers again
markers <-
  FindAllMarkers(
    object = combined_seurat,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )

cluster1_markers <-
  FindMarkers(
    object = combined_seurat,
    ident.1 = 1,
    min.pct = 0.25
  )

cluster2.markers <-
  FindMarkers(combined_seurat, ident.1 = 2, ident.2 = c(7))
head(cluster2.markers)


markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

top10
write.csv(top10, file = "top10.csv")
combined_seurat@meta.data
View(top10)
dev.off()
pdf("FeaturePlot.pdf")
FeaturePlot(
  combined_seurat,
  features = c("Adgre1", "Mrc1", "Lcn2"),
  split.by = c("group")
)
dev.off()

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 3) %>%
  ungroup() -> top3


DoHeatmap(combined_seurat, features = top3$gene) + NoLegend()

ggsave("heatmap.png", width = 10, height = 10)