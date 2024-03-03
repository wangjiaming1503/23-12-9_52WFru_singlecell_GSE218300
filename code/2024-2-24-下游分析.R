pbmc <- combined_seurat
rm(pbmc)
# Get cell and feature names, and total numbers We show multiple ways to get the same output cell names
colnames(pbmc)
Cells(pbmc)

# feature names
Features(pbmc)
rownames(pbmc)

# number of cells/features
num_cells <- ncol(pbmc)
num_features <- nrow(pbmc)

# List of object layers
Layers(pbmc)


# Variable feature names
VariableFeatures(pbmc)

# Setting and retrieving cell identities

# Set identity classes to an existing column in meta data
Idents(object = pbmc) <- "seurat_annotations"

# View cell identities, get summary table
Idents(pbmc)
table(Idents(pbmc))

# Set identity to CD4 T cells for all cells
Idents(pbmc) <- "CD4 T cells"

# Set for a selected group of cells
pbmc.cells <- Cells(pbmc)
Idents(object = pbmc, cells = pbmc.cells[1:10]) <- "CD4 T cells"

# Get cell identity classes
Idents(object = pbmc)
levels(x = pbmc)

# Stash cell identity classes in metadata
pbmc[["old.ident"]] <- Idents(object = pbmc)
pbmc <- StashIdent(object = pbmc, save.name = "old.ident")

# Rename identity classes
pbmc <- RenameIdents(object = pbmc, `CD4 T cells` = "T Helper cells")

# View metadata data frame, stored in object@meta.data
pbmc[[]]

# Retrieve specific values from the metadata
pbmc$nCount_RNA
pbmc[[c("percent.mito", "nFeature_RNA")]]

# Add metadata, see ?AddMetaData
random_group_labels <- sample(x = c("g1", "g2"), size = ncol(x = pbmc), replace = TRUE)
pbmc$groups <- random_group_labels

# Get cell embeddings and feature loadings stored on pbmc[['pca']]@cell.embeddings
Embeddings(pbmc, reduction = "pca")

# stored in pbmc[['pca]]@feature.loadings
Loadings(pbmc, reduction = "pca")