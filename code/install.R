# Enter commands in R (or R studio, if installed)

Sys.setenv(http_proxy = "http://192.168.0.105:7890")
Sys.setenv(https_proxy = "http://192.168.0.105:7890")
install.packages("remotes")
remotes::install_github('satijalab/seurat-wrappers')

install.packages('Seurat')
library(Seurat)
setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
install.packages(c("BPCells", "presto", "glmGamPoi"))
install.packages('Signac')
remotes::install_github("satijalab/seurat-data", quiet = TRUE)
remotes::install_github("satijalab/azimuth", quiet = TRUE)
remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE)

remotes::install_github("10XGenomics/loupeR")
loupeR::setup()