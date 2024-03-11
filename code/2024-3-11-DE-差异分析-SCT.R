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
  selectLab = selectLab_italics,
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
  pointSize = point_sizes,
  colAlpha = alpha_values,
  colCustom = custom_colors
)
str(volcano_plot)
volcano_plot <- volcano_plot +
  xlab(expression(log[2]*"(fold change)")) +
  ylab(expression(-lg*italic(" P")))


print(volcano_plot)

#导出 600 600 250% 