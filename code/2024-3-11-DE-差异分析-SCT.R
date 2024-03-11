library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(EnhancedVolcano)

# Load the SCT data
library(readr)
detable_MouseRNA_SCT <-
  read_csv(
    "SeuratObjects/2024-3-11-MouseRNA_SingleR.labels Features-remove-1-chow-for-granulocyte-SCT-all.csv",
    comment = "#"
  )
View(detable_MouseRNA_SCT)
colnames(detable_MouseRNA_SCT)

EnhancedVolcano(
  toptable = detable_MouseRNA_SCT,
  lab = detable_MouseRNA_SCT$FeatureName,
  x = "Hepatocytes Log2 Fold Change",
  y = "Hepatocytes P-Value",pCutoff = 0.05,
  title = NULL,
  subtitle = NULL,
  selectLab = "Lcn2"
)

#claude-------------------------------
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(EnhancedVolcano)

# Load the SCT data
library(readr)
detable_MouseRNA_SCT <-
  read_csv(
    "SeuratObjects/2024-3-11-MouseRNA_SingleR.labels Features-remove-1-chow-for-granulocyte-SCT-all.csv",
    comment = "#"
  )
View(detable_MouseRNA_SCT)
colnames(detable_MouseRNA_SCT)

lab_italics <- paste0("italic('", detable_MouseRNA_SCT$FeatureName, "')")
selectLab_italics = paste0("italic('", "Lcn2", "')")

point_sizes <- ifelse(detable_MouseRNA_SCT$FeatureName == 'Lcn2', 2.5, 1)
alpha_values <- ifelse(detable_MouseRNA_SCT$FeatureName == 'Lcn2', 1, 0.2)

custom_colors <- ifelse(detable_MouseRNA_SCT$`Hepatocytes P-Value` < 0.05 & detable_MouseRNA_SCT$`Hepatocytes Log2 Fold Change` > 0.5, 
                        'red3', 
                        ifelse(detable_MouseRNA_SCT$`Hepatocytes P-Value` < 0.05 & detable_MouseRNA_SCT$`Hepatocytes Log2 Fold Change` < -0.5, 
                               'navyblue', 
                               'grey50'))
custom_colors[is.na(custom_colors)] <- 'grey50'

names(custom_colors)[custom_colors == "red3"] <- 'high'
names(custom_colors)[custom_colors == "navyblue"] <- 'low'
names(custom_colors)[custom_colors == "grey50"] <- 'nosig'

volcano_plot <- EnhancedVolcano(
  toptable = detable_MouseRNA_SCT,
  lab = lab_italics,
  x = "Hepatocytes Log2 Fold Change",
  y = "Hepatocytes P-Value",
  pCutoff = 0.05,
  selectLab = NA, #配合下方ggrpel，如果NULL则全部画出来，不用ggrpel则为selectLab_italics
  drawConnectors = FALSE,
  boxedLabels = FALSE,
  parseLabels = TRUE,
  labCol = 'black',
  labFace = 'bold',
  title = "Hepatocytes",
  subtitle = NULL,
  ylim = c(0, -log10(min(detable_MouseRNA_SCT$`Hepatocytes P-Value`[detable_MouseRNA_SCT$`Hepatocytes P-Value` > 0]))),
  xlim = c(min(detable_MouseRNA_SCT$`Hepatocytes Log2 Fold Change`), max(detable_MouseRNA_SCT$`Hepatocytes Log2 Fold Change`)),
  caption = NULL,
  FCcutoff = 0.5,
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  legendPosition = 'none',
  #pointSize = point_sizes,
  #colAlpha = alpha_values,
  pointSize = 1,  # 设置默认点大小
  colAlpha = 0.2,  # 设置默认透明度
  colCustom = custom_colors
)
#str(volcano_plot)
volcano_plot <- volcano_plot +
  xlab(expression(log[2]*"(fold change)")) +
  ylab(expression(-lg*italic(" P")))


print(volcano_plot)
highlight_gene <- "Lcn2"

volcano_plot <- volcano_plot +
  xlab(expression(log[2]*"(fold change)")) +
  ylab(expression(-lg*italic(" P"))) +
  geom_text_repel(
    data = subset(detable_MouseRNA_SCT, FeatureName %in% c("Lcn2")),
    aes(x = `Hepatocytes Log2 Fold Change`, y = -log10(`Hepatocytes P-Value`), label = FeatureName),
    fontface = "italic",
    nudge_y = 0.5
  )
print(volcano_plot)

volcano_plot <- volcano_plot +
  geom_point(
    data = subset(detable_MouseRNA_SCT, FeatureName %in% highlight_gene),
    aes(x = `Hepatocytes Log2 Fold Change`, y = -log10(`Hepatocytes P-Value`)),
    size = 2.5,  # 设置突出显示的点的大小
    alpha = 1 , # 设置突出显示的点的透明度
    colour = custom_colors
  ) +
  xlab(expression(log[2]*"(fold change)")) +
  ylab(expression(-lg*italic(" P"))) +
  geom_text_repel(
    data = subset(detable_MouseRNA_SCT, FeatureName %in% highlight_gene),
    aes(x = `Hepatocytes Log2 Fold Change`, y = -log10(`Hepatocytes P-Value`), label = FeatureName),
    fontface = "italic",
    nudge_y = 0.5
  )
print(volcano_plot)
detable_MouseRNA_SCT[["FeatureName"]]

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
  custom_colors <- ifelse(detable_data[[paste0(cell_type, " P-Value")]] < 0.05 & detable_data[[paste0(cell_type, " Log2 Fold Change")]] > 0.5, 
                          up_color, 
                          ifelse(detable_data[[paste0(cell_type, " P-Value")]] < 0.05 & detable_data[[paste0(cell_type, " Log2 Fold Change")]] < -0.5, 
                                 down_color, 
                                 nosig_color))
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
    colAlpha = 0.2,
    colCustom = custom_colors
  )
  
  # 突出显示指定的基因
  volcano_plot <- volcano_plot +
    geom_point(
      data = subset(detable_data, FeatureName %in% highlight_gene),
      aes(x = .data[[paste0(cell_type, " Log2 Fold Change")]], y = -log10(.data[[paste0(cell_type, " P-Value")]])),
      size = 2.5,
      alpha = 1,
      colour = ifelse(subset(detable_data, FeatureName %in% highlight_gene)[[paste0(cell_type, " Log2 Fold Change")]] > 0, up_color, down_color)
    ) +
    xlab(expression(log[2]*"(fold change)")) +
    ylab(expression(-lg*italic(" P"))) +
    geom_text_repel(
      data = subset(detable_data, FeatureName %in% highlight_gene),
      aes(x = .data[[paste0(cell_type, " Log2 Fold Change")]], y = -log10(.data[[paste0(cell_type, " P-Value")]]), label = FeatureName),
      fontface = "italic",
      nudge_y = 0.5
    )
  
  return(volcano_plot)
  print(volcano_plot)
}


colnames(detable_MouseRNA_SCT)

volcano_plot_func(detable_data = detable_MouseRNA_SCT,
                  highlight_gene = c("Lcn2","Spp1"),cell_type = "Hepatocytes")

volcano_plot_func(detable_data = detable_MouseRNA_SCT,
                  highlight_gene = c("Lcn2","Lyz2"),cell_type = "Granulocytes")

save.image("/home/rstudio/work/23-12-9_52WFru_singlecell_GSE218300/Rdat
    a_time/2024年3月11日14-loupe-下游差异火山图分析.RData")
