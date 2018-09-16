#!/usr/bin/env Rscript

# Copyright (C) 2017 Tobias Jakobi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


# set mirrors
source("https://bioconductor.statistik.tu-dortmund.de/biocLite.R")
options(repos = c(CRAN = "https://cran.uni-muenster.de/"))
biocLite()

pkgs <- c(
    "aod",
    "ballgown",
    "biomaRt",
    "data.table",
    "edgeR",
    "GenomicFeatures",
    "GenomicRanges",
    "ggbio",
    "ggfortify",
    "ggplot2",
    "ggrepel",
    "gridExtra",
    "openxlsx",
    "plyr",
    "reshape2",
    "devtools",
    "kableExtra",
    "formattable",
    "dplyr",
    "RColorBrewer",
    "colortools",
    "EnsDb.Rnorvegicus.v79",
    "EnsDb.Mmusculus.v79",
    "EnsDb.Hsapiens.v86"
)

# check if devtools is already installed
pkgs <- pkgs[!pkgs %in% installed.packages()[,1]]
if (length(pkgs) > 0)
  biocLite(pkgs)

# load devtools library
library(devtools)

# install CircTest from the Dieterich Lab GitHub page from master branch
install_github("dieterich-lab/CircTest", ref = "master")

# install primex from the Dieterich Lab GitHub page from master branch
install_github("dieterich-lab/primex", ref = "master")
