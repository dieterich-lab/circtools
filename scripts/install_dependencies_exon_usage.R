#!/usr/bin/env Rscript

# set mirrors
source("https://bioconductor.statistik.tu-dortmund.de/biocLite.R")
options(repos = c(CRAN = "https://cran.uni-muenster.de/"))

# we need devtools as requirement for circtest
pkgs <- c(
    "devtools",
    "ballgown",
    "edgeR",
    "biomaRt",
    "ggbio",
    "ggfortify",
    "openxlsx",
    "GenomicRanges",
    "GenomicFeatures"
)

# check if devtools is already installed
pkgs <- pkgs[!pkgs %in% installed.packages()[,1]]
if (length(pkgs) > 0)
  biocLite(pkgs)
