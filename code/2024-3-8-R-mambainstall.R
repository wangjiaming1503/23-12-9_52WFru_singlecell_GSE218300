.libPaths("/home/rstudio/work/miniforge3/envs/r4.3.2/lib/R/library")
.libPaths()

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# Configure BioCManager to use Posit Package Manager:
options(BioC_mirror = "https://packagemanager.posit.co/bioconductor")
options(BIOCONDUCTOR_CONFIG_FILE = "https://packagemanager.posit.co/bioconductor/config.yaml")

# Configure a CRAN snapshot compatible with Bioconductor 3.18:
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"))
chooseBioCmirror()
library(BiocManager)
BiocManager::install(c("simpleSingleCell","RNAseq123","rnaseqGene"))
BiocManager::install(c("gert","RNAseq123","usethis","BiocWorkflowTools"))
BiocManager::install("RNAseq123")
#mamba install r-gert
pkgs <- c(
  "liftOver",
  "rnaseqGene",
  "methylationArrayAnalysis",
  "RnaSeqGeneEdgeRQL",
  "RNAseq123",
  "GeoMxWorkflows",
  "TCGAWorkflow",
  "annotation",
  "cytofWorkflow",
  "maEndToEnd",
  "simpleSingleCell",
  "rnaseqDTU",
  "arrays",
  "sequencing",
  "ExpressionNormalizationWorkflow",
  "variants",
  "chipseqDB",
  "BP4RNAseq",
  "generegulation",
  "highthroughputassays",
  "ExpHunterSuite",
  "recountWorkflow",
  "BiocMetaWorkflow",
  "CAGEWorkflow",
  "csawUsersGuide",
  "EGSEA123",
  "fluentGenomics",
  "SingscoreAMLMutations"
)

c("arrangements","MKinfer","methylationArrayAnalysis","curl")


BiocManager::install(c("arrangements","MKinfer","methylationArrayAnalysis","curl"), ask = FALSE, update = TRUE)

install.packages("arrangements")