#if (as.numeric(R.Version()$minor) *.1 + as.numeric(R.Version()$major) < 3.4)
#  stop("Please update R to the version 3.4 or later")

source("https://bioconductor.org/biocLite.R")

options(repos = c(CRAN = "http://cran.rstudio.com"))

pkgs <- c(
  #"AnnotationFilter",
  "Rsamtools",
  "htmltools",
  "DECIPHER",
  "S4Vectors",
  "IRanges",
  "Biostrings",
  "BSgenome",
  "GenomicFeatures",
  "ensembldb"
)

gitpkgs <- c(
  "AnnotationFilter"
)

pkgs <- pkgs[!pkgs %in% installed.packages()[,1]]
if (length(pkgs) > 0)
  biocLite(pkgs)

install.packages("devtools")
library(devtools)

#install_github("rstats-db/RSQLite", ref = "v1.1-15")

install_github("dieterich-lab/circtools",
             subdir = "R/circtools",
             ref = "master")