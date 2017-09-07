source("https://bioconductor.org/biocLite.R")
options(repos = c(CRAN = "http://cran.rstudio.com"))

installRSQL <- FALSE
if (requireNamespace("RSQLite")) {
  if (packageVersion("RSQLite") > "1.1.15") {
    answer <- "!"
    while (!answer %in% c("y", "n")) {
      cat(
        paste0(
          "The current version of circtools can work ",
          "only with RSQLite version <= 1.1.5\n",
          "Your version is ", packageVersion("RSQLite"), "\n",
          "Would you like to install the 1.1.15 one? [y/n]:  "
        )
      )
      answer <- readLines(con = "stdin", n = 1)
      if (answer %in% c("n"))
        quit()
      if (answer %in% c("y"))
        installRSQL <- TRUE
    }
  }
}

if (!requireNamespace("devtools")) {
  install.packages("devtools")
}

if (!requireNamespace("devtools")) {
  install.packages("devtools")
}

if (installRSQL)
  devtools::install_github("rstats-db/RSQLite", ref = "v1.1-15")

pkgs <- c(
  "Rsamtools",
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
             ref = "master")