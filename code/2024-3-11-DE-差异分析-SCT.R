# Load the SCT data
library(readr)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(EnhancedVolcano)
setwd("/home/rstudio/work/23-12-9_52WFru_singlecell_GSE218300")


detable_MouseRNA_SCT <-
  read_csv(
    "2024年3月27日active_cluster Features.csv",
    comment = "#"
  )

#打包函数化火山图重构===========================--------------------------------------------------------------------------




volcano_plot_func <- function(detable_data,
                              cell_type,
                              highlight_gene,
                              FeatureName = "FeatureName",
                              selectLab = NA,
                              up_color = "red3",
                              down_color = "navyblue",
                              nosig_color = "grey50") {
  # 创建自定义颜色向量
  custom_colors <-
    ifelse(
      detable_data[[paste0(cell_type, " P-Value")]] < 0.05 &
        detable_data[[paste0(cell_type, " Log2 Fold Change")]] > 0.5,
      up_color,
      ifelse(
        detable_data[[paste0(cell_type, " P-Value")]] < 0.05 &
          detable_data[[paste0(cell_type, " Log2 Fold Change")]] < -0.5,
        down_color,
        nosig_color
      )
    )
  custom_colors[is.na(custom_colors)] <- nosig_color
  
  names(custom_colors)[custom_colors == up_color] <- 'high'
  names(custom_colors)[custom_colors == down_color] <- 'low'
  names(custom_colors)[custom_colors == nosig_color] <- 'nosig'
  
  # 创建火山图
  volcano_plot <- EnhancedVolcano(
    toptable = detable_data,
    lab = detable_data[[FeatureName]],
    x = paste0(cell_type, " Log2 Fold Change"),
    y = paste0(cell_type, " P-Value"),
    pCutoff = 0.05,
    selectLab = selectLab,
    drawConnectors = FALSE,
    boxedLabels = FALSE,
    parseLabels = TRUE,
    labCol = 'black',
    labFace = 'bold',
    title = cell_type,
    subtitle = NULL,
    ylim = c(0, -log10(min(detable_data[[paste0(cell_type, " P-Value")]][detable_data[[paste0(cell_type, " P-Value")]] > 0]))),
    xlim = c(min(detable_data[[paste0(cell_type, " Log2 Fold Change")]]), max(detable_data[[paste0(cell_type, " Log2 Fold Change")]])),
    caption = NULL,
    FCcutoff = 0.5,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    legendPosition = 'none',
    pointSize = 1,
    colAlpha = 0.1,
    colCustom = custom_colors
  )
  
  # 突出显示指定的基因
  volcano_plot <- volcano_plot +
    geom_point(
      data = subset(detable_data, FeatureName %in% highlight_gene),
      aes(x = .data[[paste0(cell_type, " Log2 Fold Change")]], y = -log10(.data[[paste0(cell_type, " P-Value")]])),
      size = 2.5,
      alpha = 1,
      colour = ifelse(
        subset(detable_data, FeatureName %in% highlight_gene)[[paste0(cell_type, " Log2 Fold Change")]] > 0,
        up_color,
        down_color
      )
    ) +
    xlab(expression(log[2] * "(fold change)")) +
    ylab(expression(-lg * italic(" P"))) +
    geom_text_repel(
      data = subset(detable_data, FeatureName %in% highlight_gene),
      aes(
        x = .data[[paste0(cell_type, " Log2 Fold Change")]],
        y = -log10(.data[[paste0(cell_type, " P-Value")]]),
        label = FeatureName
      ),
      fontface = "italic",
      nudge_y = 0.5
    )
  
  return(volcano_plot)
  print(volcano_plot)
}

colnames(detable_MouseRNA_SCT)

volcano_plot_func(
  detable_data = detable_MouseRNA_SCT,
  highlight_gene = c("Lcn2"),
  cell_type = "Hepatocytes"
)

volcano_plot_func(
  detable_data = detable_MouseRNA_SCT,
  highlight_gene = c("Lcn2", "Lyz2"),
  cell_type = "Granulocytes"
)

#save.image(
#   "/home/rstudio/work/23-12-9_52WFru_singlecell_GSE218300/Rdat
#     a_time/2024年3月11日14-loupe-下游差异火山图分析.RData"
# )
load("/home/rstudio/work/23-12-9_52WFru_singlecell_GSE218300/Rdata_time/2024年3月11日14-loupe-下游差异火山图分析.RData")


volcano_plot_func(
  detable_data = detable_MouseRNA_SCT,
  highlight_gene = c("Lcn2", "Lyz2","Ly6g"),
  cell_type = "Granulocytes"
)

volcano_plot_func(
  detable_data = detable_MouseRNA_SCT,
  highlight_gene = c("Mmp2"),
  cell_type = "Fibroblasts"
)

volcano_plot_func(
  detable_data = detable_MouseRNA_SCT,
  highlight_gene = c("Saa1","Saa2"),
  cell_type = "Fibroblasts"
)

volcano_plot_func(
  detable_data = detable_MouseRNA_SCT,
  highlight_gene = c("Lcn2", "Spp1","Saa1","Saa2"),
  cell_type = "Hepatocytes"
)
#打开绘图
library(httpgd)
hgd()
hgd_view()

#monocle3===========================================
library(monocle3)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)


summary(combined_seurat@reductions)
colnames(combined_seurat@meta.data)


cds <- as.cell_data_set(
  combined_seurat,
  default.reduction = "tsne.sct.harmony"
)
cds <- reduce_dimension(cds)

cds <- cluster_cells(cds)

p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

DimPlot(combined_seurat, group.by = c("MouseRNA_SingleR.labels"),reduction = "tsne.sct.harmony")
str(cds)

save(cds,file = paste0("SeuratObjects/",get_time(),"cds.Rdata"))
