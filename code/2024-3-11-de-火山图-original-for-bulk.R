# library(clusterProfiler)
# library(org.Mm.eg.db)
# BiocManager::install("EnhancedVolcano")
# devtools::install_github('kevinblighe/EnhancedVolcano')


library(EnhancedVolcano)
library(readr)
## ----setup2, message=FALSE, eval=TRUE---------------------------------------------------
library(limma)
# library(Glimma)
library(edgeR)
# library(Mus.musculus)
# library(htmlwidgets)
load("~/bioinfo/23-02-12肝脏/2023年2月13日 4V5 有权重版/4v5 weitht.RData")
treatedvsuntreated <- topTable(efit, coef=1, n=Inf)

# write.csv(treatedvsuntreated,"treatedvsuntreated.csv")
lab_italics <- paste0("italic('", treatedvsuntreated$SYMBOL, "')")#基因变成斜体
selectLab_italics = paste0(
  "italic('",
  c('Lcn2'),
  "')")
# 默认带背景
EnhancedVolcano(treatedvsuntreated ,
                lab = lab_italics ,
                x = 'logFC',
                y = 'P.Value',
                pCutoff = 	
                  0.05,
                 selectLab =selectLab_italics,
                drawConnectors = TRUE,boxedLabels = FALSE,parseLabels = TRUE,
                labCol = 'black',
                labFace = 'bold',
                title = NULL,
                subtitle = NULL,ylim = c(0, 12),xlim = c(-8, 8),
                legendLabSize = 10,caption = NULL
                ,FCcutoff = 0.5
                  )
# 默认配色，去掉背景
EnhancedVolcano(treatedvsuntreated,
                lab = lab_italics,
                x = 'logFC',
                y = 'P.Value',
                pCutoff = 0.05,
                selectLab = selectLab_italics,
                drawConnectors = TRUE,
                boxedLabels = FALSE,
                parseLabels = TRUE,
                labCol = 'black',
                labFace = 'bold',
                title = NULL,
                subtitle = NULL,
                ylim = c(0, 12),
                # xlim = c(-8, 8),
                legendLabSize = 0,
                caption = NULL,
                FCcutoff = 0.5,gridlines.major = FALSE,
                gridlines.minor = FALSE,legendPosition = 'none') 
??EnhancedVolcano



EnhancedVolcano(treatedvsuntreated,
                lab = lab_italics,
                x = 'logFC',
                y = 'P.Value',
                pCutoff = 0.05,
                selectLab = selectLab_italics,
                drawConnectors = TRUE,
                boxedLabels = FALSE,
                parseLabels = TRUE,
                labCol = 'black',
                labFace = 'bold',
                title = NULL,
                subtitle = NULL,
                ylim = c(0, 10.5),
                xlim = c(-8, 6),
                # legendLabSize = 0,
                caption = NULL,
                FCcutoff = 0.5,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                legendPosition = 'none',
                pointSize = 0.5, # decreasing point size
                colAlpha = 0.5, # adding some transparency
                col =c("grey50","grey50", "grey50", "red3")
)

point_sizes <- ifelse(treatedvsuntreated$SYMBOL == 'Lcn2', 2.5, 1)
alpha_values <- ifelse(treatedvsuntreated$SYMBOL == 'Lcn2', 1, 0.2)



library(ggplot2)

volcano_plot <- EnhancedVolcano(treatedvsuntreated,
                                lab = lab_italics,
                                x = 'logFC',
                                y = 'P.Value',
                                pCutoff = 0.05,
                                selectLab = selectLab_italics,
                                drawConnectors = FALSE,
                                boxedLabels = FALSE,
                                parseLabels = TRUE,
                                labCol = 'black',
                                labFace = 'bold',
                                title = NULL,
                                subtitle = NULL,
                                ylim = c(0, 10.5),
                                xlim = c(-8, 6),
                                # legendLabSize = 0,
                                caption = NULL,
                                FCcutoff = 0.5,
                                gridlines.major = FALSE,
                                gridlines.minor = FALSE,
                                legendPosition = 'none',
                                pointSize = point_sizes, # use the size vector here
                                colAlpha = alpha_values, # use the alpha vector here
                                col =c("grey50","grey50", "grey50", "red3")
)

# Change x and y axis labels
volcano_plot <- volcano_plot +
  xlab(expression(log[2]*"(fold change)")) +
  ylab(expression(lg*italic(" P")))

print(volcano_plot)
# 红蓝上下调配色

# 创建自定义颜色向量
custom_colors <- ifelse(treatedvsuntreated$P.Value < 0.05 & treatedvsuntreated$logFC > 0.5, 
                        'red3', 
                        ifelse(treatedvsuntreated$P.Value < 0.05 & treatedvsuntreated$logFC < -0.5, 
                               'navyblue', 
                               'grey50'))
custom_colors[is.na(custom_colors)]<- 'grey50'

# 这里需要名字，不然会报错
names(custom_colors)[custom_colors == "red3"] <- 'high'
names(custom_colors)[custom_colors == "navyblue"] <- 'low'
names(custom_colors)[custom_colors == "grey50"] <- 'nosig'
# 在绘制火山图的函数中使用这个自定义颜色向量
volcano_plot <- EnhancedVolcano(treatedvsuntreated,
                                lab = lab_italics,
                                x = 'logFC',
                                y = 'P.Value',
                                pCutoff = 0.05,
                                selectLab = selectLab_italics,
                                drawConnectors = FALSE,
                                boxedLabels = FALSE,
                                parseLabels = TRUE,
                                labCol = 'black',
                                labFace = 'bold',
                                title = NULL,
                                subtitle = NULL,
                                ylim = c(0, 10.5),
                                xlim = c(-8, 6),
                                caption = NULL,
                                FCcutoff = 0.5,
                                gridlines.major = FALSE,
                                gridlines.minor = FALSE,
                                legendPosition = 'none',
                                pointSize = point_sizes, # 使用这里的大小向量
                                colAlpha = alpha_values, # 使用这里的alpha向量
                                colCustom  = custom_colors  # 使用这里的自定义颜色向量
)
# Change x and y axis labels
volcano_plot <- volcano_plot +
  xlab(expression(log[2]*"(fold change)")) +
  ylab(expression(-lg*italic(" P")))

print(volcano_plot)

EnhancedVolcano(treatedvsuntreated ,
                lab = treatedvsuntreated$SYMBOL ,
                x = 'logFC',
                y = 'P.Value',pCutoff = 	
                  0.10,selectLab= c("Fgf15"))
