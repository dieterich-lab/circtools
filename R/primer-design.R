library(circtools)

ensPackage <- "EnsDb.Hsapiens.v86"
circFile <- "circs.tab"
typeExons <- "all"
bsgPackage <- "BSgenome.Hsapiens.NCBI.GRCh38"
reportFile <- "report.html"

# check for installed annotation
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
  db <- get(ensPackage, envir = asNamespace(ensPackage))
} else {
  stop(paste(
    "The gtf support is not yet implemented. \n",
    "Please construct an ensembldb object from the gtf file."))
}

tableToGR <- function(table,
                         columns = list(
                           chr    = 1,
                           start  = 2,
                           end    = 3,
                           strand = 4
                         )) {
  gr <- GenomicRanges::GRanges(
    seqnames = table[, columns$chr],
    ranges = IRanges::IRanges(start = table[, columns$start],
                              end   = table[, columns$end]),
    strand = table[, columns$strand]
  )
  gr$sjId <- paste0(table[, columns$chr], ":",
                    table[, columns$start], "-",
                    table[, columns$end])
  gr
}

circs <- read.table(file = circFile, header =TRUE,sep = " ")
circData <- CircData(db, tableToGR(circs, list(chr=1, start = 2, end = 3, strand = 5)))


if (!is.null(bsgPackage)) {
  a <- requireNamespace(ensPackage, quietly = TRUE)
  if (!a) {
    cat(
      paste0(
        "The specified package is not available. \n",
        "Please check if the name is correct: ",
        bsgPackage,
        "\n"
      )
    )
    cat(
      paste0(
        "\nIn order to install the missing Bioconductor package, \n",
        "please run in the R prompt:\n",
        "source(\"https://bioconductor.org/biocLite.R\")\n",
        "biocLite(\"",
        bsgPackage,
        "\""
      )
    )
  }
  bsg <- get(bsgPackage, envir = asNamespace(ensPackage))
} else {
  # check for fasta
}

exSeqs <- getExonSeqs(circData = circData, bsg = bsg, type = typeExons)
if (!is.null(reportFile))
  reportCircs(exSeq = exSeqs, file = reportFile)

primers <- designPrimers(exSeq = exSeqs, db = db, bsg = bsg)
