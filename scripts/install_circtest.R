#!/usr/bin/env Rscript

# set mirrors
source("https://bioconductor.statistik.tu-dortmund.de/biocLite.R")
options(repos = c(CRAN = "https://cran.uni-muenster.de/"))

# we need devtools as requirement for circtest
pkgs <- c(
  "devtools"

)

# check if devtools is already installed
pkgs <- pkgs[!pkgs %in% installed.packages()[,1]]
if (length(pkgs) > 0)
  biocLite(pkgs)

# load devtools library
library(devtools)

# install circtest from the dieterich lab github page from master branch
install_github("dieterich-lab/CircTest", ref = "master")