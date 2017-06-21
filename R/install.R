source("https://bioconductor.org/biocLite.R")

pkgs <- c(
  "AnnotationFilter",
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

pkgs <- pkgs[!pkgs %in% installed.packages()[,1]]
if (length(pkgs) > 0)
  biocLite(pkgs)

install.packages("devtools")
library(devtools)
install_github("dieterich-lab/circtools",
             subdir = "R/circtools",
             ref = "r-dev")