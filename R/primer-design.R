library(circtools)

ensPackage <- "EnsDb.Hsapiens.v86"

if (!is.null(ensPackage)) {
  a <- requireNamespace(ensPackage, quietly = TRUE)
  if (!a) {
    cat(
      paste0(
        "The specified package is not available. \n",
        "Please check if the name is correct: ",
        ensPackage, "\n"
      )
    )
    cat(
      paste0(
        "\nIn order to install the missing Bioconductor package, \n",
        "please run in the R prompt:\n",
        "source(\"https://bioconductor.org/biocLite.R\")\n",
        "biocLite(\"",
        ensPackage,
        "\""
      )
    )
  }
}
