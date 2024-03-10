library(Seurat)
library(SeuratWrappers)
options(future.globals.maxSize = 1e9)

load("SeuratObjects/2024-03-10-19-38combined_seurat-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden_all_joined-tsne-anotated.Rdata")


combined_seurat@reductions
summary(combined_seurat@reductions)
combined_seurat <- RunTSNE(
  combined_seurat,
  reduction = "integrated.sct.cca",
  cells = NULL,
  dims = 1:30,
  features = NULL,
  seed.use = 1,
  tsne.method = "Rtsne",
  dim.embed = 2,
  distance.matrix = NULL,
  reduction.name = "tsne.sct.cca",
  reduction.key = "sct.cca.tSNE_"
)
#setting key from sct.cca.tSNE_ to sctccatSNE_
combined_seurat <- RunTSNE(
  combined_seurat,
  reduction = "integrated.sct.rpca",
  dims = 1:30,
  reduction.name = "tsne.sct.rpca"
)
DefaultAssay(combined_seurat)


save(combined_seurat,"./SeuratObjects/",get_time()"combined_seurat-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden_all_joined-tsne-anotated-alltsne")
library("loupeR")
create_loupe_from_seurat(
  combined_seurat,
  output_name = paste0(
    "SeuratObjects/",get_time(),"combined_seurat-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden_all_joined-tsne-anotated-alltsne",
    format(Sys.Date(), "%Y-%m-%d")
  )
)
DefaultAssay(combined_seurat) <- "SCT"
create_loupe_from_seurat(
  combined_seurat,
  output_name = paste0(
    "SeuratObjects/",get_time(),"combined_seurat-SCVI-intergraed-umap-add-learn-N4-slm-louvain-muti-leiden_all_joined-tsne-anotated-alltsne-SCT",
    format(Sys.Date(), "%Y-%m-%d")
  )
)


#3dtsne===========================================================================================
# Install plot_ly
install.packages('plotly')

# Load plot_ly
library(plotly)

# Construct a dataframe using data from your pre-clustered Seurat v3.0.0 object
# Here 'seurat_clusters' is list of numeric cluster identities, you can find it here: yourseuratobject[["seurat_cluster"]], 
# or yourseuratobject$seurat_clusters, where 'yourseuratobject' is a Seurat object created with Seurat v3.1.1 (works for v3.0.0 as well)
yourseuratobject <- ThisIsWhateverYourSeuratObjectIs

# Re-run tSNE s that you have accurate calculations for all tSNE(s)
yourseuratobject <- RunTSNE(yourseuratobject,
                        reduction.use = "pca",
                        dims.use = 1:20,
                        dim.embed = 3)

# Extract tSNE information from Seurat Object
tsne_1 <- yourseuratobject[["tsne"]]@cell.embeddings[,1]
tsne_2 <- yourseuratobject[["tsne"]]@cell.embeddings[,2]
tsne_3 <- yourseuratobject[["tsne"]]@cell.embeddings[,3]

# Prepare a dataframe for cell plotting
plot.data <- FetchData(object = yourseuratobject, vars = c("tSNE_1", "tSNE_2", "tSNE_3", "seurat_clusters"))

# Make a column of row name identities (these will be your cell/barcode names)
plot.data$label <- paste(rownames(plot.data))

# Plot your data, in this example my Seurat object had 21 clusters (0-20)
plot_ly(data = plotting.data, 
        x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
        color = ~seurat_clusters, 
        colors = c("lightseagreen",
                   "gray50",
                   "darkgreen",
                   "red4",
                   "red",
                   "turquoise4",
                   "black",
                   "yellow4",
                   "royalblue1",
                   "lightcyan3",
                   "peachpuff3",
                   "khaki3",
                   "gray20",
                   "orange2",
                   "royalblue4",
                   "yellow3",
                   "gray80",
                   "darkorchid1",
                   "lawngreen",
                   "plum2",
                   "darkmagenta"),
        type = "scatter3d", 
        mode = "markers", 
        marker = list(size = 5, width=2), # controls size of points
        text=~label, #This is that extra column we made earlier for which we will use
        hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names

# Say you wanto make a gene-expression 3D plot, where you can plot gene expression against a color scale
# Here using the same seurat object as above, we extract gene expression information for beta-actin 'ACTB'
# create a dataframe
plotting.data <- FetchData(object = yourseuratobject, vars = c("tSNE_1", "tSNE_2", "tSNE_3", "ACTB"))

# Say you want change the scale, so that every cell having an expression >1 will be one color
# Basically, you are re-adjusting the scale here, so that any cell having a certain expression 
# will light up on your 3D plot

# First make another column in your dataframe, where all values above 1 are re-assigned a value of 1
# This information is stored in the 'changed' column of your dataframe
plotting.data$changed <- ifelse(test = plotting.data$ACTB <1, yes = plotting.data$ACTB, no = 1)

# Add the label column, so that now the column has 'cellname-its expression value'
plotting.data$label <- paste(rownames(plotting.data)," - ", plotting.data$ACTB, sep="")

# Plot your data, in this example my Seurat object had 21 clusters (0-20), and cells express a gene called ACTB
plot_ly(data = plotting.data, 
        x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
        color = ~changed, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
        opacity = .5,
        colors = c('darkgreen', 'red'), 
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 5, width=2), 
        text=~label,
        hoverinfo="text"
)

