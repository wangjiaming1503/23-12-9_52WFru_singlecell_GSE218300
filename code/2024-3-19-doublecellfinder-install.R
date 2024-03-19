
#环境安装，DblFinder------------------------------
# BiocManager::install("scater")
# BiocManager::install("scDblFinder")
# BiocManager::install("BiocParallel")
# BiocManager::install("GenomeInfoDb")
# BiocManager::valid()
# BiocManager::install(c(
#     "BiocVersion", "brew", "brio", "bslib", "callr", "cli", "commonmark",
#     "cpp11", "curl", "desc", "digest", "EGSEA", "EGSEA123", "evaluate",
#     "fansi", "gert", "ggsci", "glue", "htmltools", "htmlwidgets", "httpuv",
#     "httr2", "jsonlite", "knitr", "later", "lifecycle", "littler", "pkgbuild",
#     "pkgload", "processx", "ps", "psych", "ragg", "Rcpp", "remotes", "rlang",
#     "rmarkdown", "roxygen2", "rprojroot", "sass", "shiny", "stringi",
#     "stringr", "systemfonts", "testthat", "tinytex", "usethis", "vctrs",
#     "waldo", "withr", "xfun", "xml2", "yaml", "zip"
#   ), update = TRUE, ask = FALSE, force = TRUE)

#anadata to sce-------------------------------------------------
# BiocManager::install("zellkonverter")

View(maartenutils::gen_file_overview("SeuratObjects"))
library(zellkonverter)

sce <- zellkonverter::readH5AD("./SeuratObjects/2024-03-19-10-45anotated-alltsne-combined_seurat-quality_control.h5ad", verbose = TRUE)




library(reticulate)
save("sce", file = "./SeuratObjects/2024-03-19-12-40anotated-alltsne-combined_seurat-quality_control_sce.RData")
#use_condaenv("singlecell", required = TRUE)

#anadata read---------------------------------------------------
anndata <- reticulate::import("anndata")
adata<- anndata$read_h5ad("./SeuratObjects/2024-03-19-10-45anotated-alltsne-combined_seurat-quality_control.h5ad")
#explore adata---------------------------------------------------
head(adata$obs)
str(adata)

#doublet detection----------------------------------------------
set.seed(123)

sce = scDblFinder(sce)

doublet_score = sce$scDblFinder.score
doublet_class = sce$scDblFinder.class
adata$obs["scDblFinder_score"] = doublet_score
adata$obs["scDblFinder_class"] = doublet_class
#anada write---------------------------------------------------
adata$write_h5ad("./SeuratObjects/2024-03-19-19-05anotated-alltsne-combined_seurat-quality_control_doublet.h5ad")
#sce to h5ad write-----------------------------------------------------
zellkonverter::writeH5AD(sce, "./SeuratObjects/2024-03-19-19-09anotated-alltsne-combined_seurat-quality_control_doublet_scetrans.h5ad", verbose = TRUE)
#save sce
save("sce", file = "./SeuratObjects/2024-03-19-19-07anotated-alltsne-combined_seurat-quality_control_sce_doublet.RData")


