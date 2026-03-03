options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)

  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(GO.db)

  library(yaml)
  library(circlize)
  library(RColorBrewer)
})

set.seed(0)

source("R/utils_io.R")
source("R/utils_go.R")
source("R/utils_plot.R")