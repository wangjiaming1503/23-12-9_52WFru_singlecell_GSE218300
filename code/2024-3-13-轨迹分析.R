# claude
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(SingleCellExperiment)

View(maartenutils::gen_file_overview("SeuratObjects"))
View(maartenutils::gen_file_overview("sceobject"))
load(
  "SeuratObjects/2024-03-11-00-49combined_seurat-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden_all_joined-tsne-anotated-alltsne.Rdata"
)


load("./sceobject/2024-03-13-15-44miloobj-all_joined-tsne-anotated-alltsne.Rdata")
load("sceobject/2024-03-13-14-52-SCE-combined_seurat-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden_all_joined-tsne-anotated-alltsne.Rdata")

# 将 Seurat 对象转换为 SingleCellExperiment 对象
sce <- as.SingleCellExperiment(combined_seurat)
str(sce, max.level = 2)
# 现在可以将 SingleCellExperiment 对象传递给 Milo
milo_obj <- Milo(sce)

# 使用miloR比较亚群丰度差异
# BiocManager::install("miloR")
library(miloR)

save(
  sce,
  file = paste0(
    "sceobject/",
    get_time(),
    "-SCE-combined_seurat-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden_all_joined-tsne-anotated-alltsne.Rdata"
  )
)


load(
  "./sceobject/2024-03-13-14-52-SCE-combined_seurat-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden_all_joined-tsne-anotated-alltsne.Rdata"
)


sce$celltype.group <-
  paste(sce$MouseRNA_SingleR.labels, sce$group, sep = "_")



colnames(colData(sce))
milo.meta  <- colData(sce)
milo_obj <- buildGraph(milo_obj, k = 20, d = 30)
milo_obj <- makeNhoods(
  milo_obj,
  k = 20,
  d = 30,
  refined = TRUE,
  prop = 0.2
)
milo_obj <- calcNhoodDistance(milo_obj, d = 30)
milo_obj <-
  countCells(milo_obj, samples = "celltype.group", meta.data = colData(sce))
milo.obj  <- milo_obj
save(
  milo_obj,
  file = paste0(
    "./sceobject/",
    get_time(),
    "miloobj-all_joined-tsne-anotated-alltsne.Rdata"
  )
)
#需要把分组插入标注中
milo.design <-
  as.data.frame(xtabs( ~ celltype.group + group, data = milo.meta))
milo.design <- milo.design[milo.design$Freq > 0,]
rownames(milo.design) <- milo.design$celltype.group
milo.design <- milo.design[colnames(nhoodCounts(milo.obj)), ]

milo.res <-
  testNhoods(milo.obj,
             design =  ~ celltype.group + group,
             design.df = milo.design)
head(milo.res)





milo_obj <-
  countCells(milo_obj, samples = "MouseRNA_SingleR.labels", meta.data = colData(sce))
da_results <-
  testNhoods(milo_obj, design = ~ MouseRNA_SingleR.labels + group)
#milo_obj <- calcNhoodDistance(milo_obj, d=30)



milo_obj <-
  countCells(milo_obj, samples = "orig.ident", meta.data = colData(sce))


da_results <-
  testNhoods(milo_obj,
             design = ~ orig.ident + group ,
             data = colData(sce))
milo_obj <- buildNhoodGraph(milo_obj)





# 使用pySCENIC进行转录因子和靶基因分析
library(reticulate)
use_virtualenv("path/to/pySCENIC_env")
importpySCENIC()
scenicOptions <- initializeScenic(org = "mmu")
genesKept <-
  geneFiltering(
    as.matrix(sc_obj@assays$RNA@counts),
    scenicOptions,
    minCountsPerGene = 3 * .01 * ncol(sc_obj@assays$RNA)
  )
exprMat_filtered <- as.matrix(sc_obj@assays$RNA@counts)[genesKept,]
runCorrelation(exprMat_filtered, scenicOptions)
runGenie3(exprMat_filtered, scenicOptions)
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)

# 使用scVelo推断mononuclear phagocyte的转换
library(scvelo)
ldata <- sc_obj@assays$RNA@data
spliced <- ldata[grepl("^Exon", rownames(ldata)),]
unspliced <- ldata[grepl("^Intron", rownames(ldata)),]
ambiguous <- ldata[grepl("^Ambiguous", rownames(ldata)),]
cell.dist <- as.dist(1 - armaCor(t(ldata)))
emb <- Rtsne::Rtsne(as.matrix(cell.dist))$Y
sc_obj@reductions$tsne <- emb
scv <- readVelocity(sc_obj, spliced, unspliced, ambiguous)
scv <- gene.relative.velocity.estimates(scv)
scv <- velocity_graph(scv)
cell.dist <- as.dist(1 - armaCor(t(ldata)))
show.velocity.on.embedding.cor(
  scv,
  n = 200,
  scale = "sqrt",
  cell.colors = ac(traj.col, alpha = 0.5),
  cex = 0.8,
  arrow.scale = 3,
  show.grid.flow = TRUE,
  min.grid.cell.mass = 0.5,
  grid.n = 40,
  arrow.lwd = 1,
  do.par = FALSE,
  cell.border.alpha = 0.1
)

# 使用Destiny和Slingshot进行轨迹分析
library(scater)
library(destiny)
library(slingshot)
dm <- DiffusionMap(sc_obj)
sc_obj@reducedDims$dm <- dm@eigenvectors

clusters <- sc_obj@active.ident
sds <-
  getLineages(
    sc_obj,
    reducedDim = "dm",
    clusterLabels = clusters,
    start.clus = "MP"
  )

# 使用NicheNet预测ligand-receptor相互作用
library(nichenetr)
lr_network <-
  readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
weighted_networks <-
  readRDS(url(
    "https://zenodo.org/record/3260758/files/weighted_networks.rds"
  ))
ligand_target_matrix <-
  readRDS(url(
    "https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"
  ))
real_links <- lr_network %>%
  filter(database == "omnipath") %>%
  pull(interaction_name)
ligands <- weighted_networks$lr_sig %>%
  pull(from) %>%
  unique()
receptors <- weighted_networks$lr_sig %>%
  pull(to) %>%
  unique()

### 对特定细胞群使用ligand_activity_heatmap看配体表达
active_ligands <- ligand_activity_heatmap(
  geneset = marker_genes,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33,
  ligands = ligands,
  color_low = "whitesmoke",
  color_high = "red"
)
active_ligands_df <- active_ligands$ligand_activity %>%
  as.data.frame() %>%
  t()

### 对特定细胞群使用receptor_activity_heatmap看受体表达
active_receptors <- receptor_activity_heatmap(
  geneset = marker_genes,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33,
  receptors = receptors,
  color_low = "whitesmoke",
  color_high = "red"
)
