clArgs <- commandArgs(trailingOnly = TRUE)

default <- list(
  # input
  circFile    = NULL, # "circs.tsv",
  ensPackage  = NULL, # "EnsDb.Hsapiens.v86",
  bsgPackage  = NULL, # "BSgenome.Hsapiens.NCBI.GRCh38",
  # processing 
  typeExons   = "longest",
  # output
  # where to write an HTML report with exon sequences
  reportFile  = "report.html",
  primerFile  = "primers.tsv",
  productFile = "products.tsv",
  # where to save the result of the designPrimers R function as an RDS object
  rdsFile     = "rds-",
  # separator in input and output files
  sep         = "\t"
)

parseArgs <- function(args, default) {
  args <- matrix(args, ncol = 2, byrow = TRUE)
  # remove "--"
  args[, 1] <- substring(args[, 1], 3)
  extraArgs <- setdiff(args[,1], names(default))
  if (length(extraArgs) > 0)
    stop("Wrong arguments: ", paste(extraArgs, collapse = ", "), "!")
  default[args[,1]] <- args[,2]
  default
}

default <- parseArgs(clArgs, default)

invisible(lapply(names(default), function(x) {
  assign(x, default[[x]], envir = globalenv())
}))

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

suppressPackageStartupMessages(library(circtools))
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

cat("Reading the circs...")
circs <- read.table(file = circFile, sep = sep)
circData <- CircData(db, tableToGR(circs, list(
  chr = 1,
  start = 2,
  end = 3,
  strand = 4
)))
cat("OK!\n")


if (!is.null(bsgPackage)) {
  a <- requireNamespace(bsgPackage, quietly = TRUE)
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
  bsg <- get(bsgPackage, envir = asNamespace(bsgPackage))
} else {
  # check for fasta
}

cat("Quering exon sequences....")
exSeqs <- getExonSeqs(circData = circData, bsg = bsg, type = "shortest")
if (!is.null(reportFile))
  reportCircs(exSeq = exSeqs, file = reportFile)
cat("OK!\n")

cat("Primer design...")
res <- designPrimers(exSeq = exSeqs, db = db, bsg = bsg)
cat("OK!\n")

cat("writting.....")
primerList <- unlist(res$primers, recursive = FALSE)
primerList <- do.call(rbind, lapply(primerList,
                                    function(x) {
                                      x <- as.data.frame(x)
                                      x[order(x$start), ]
                                    }))
write.table(
  primerList,
  file = primerFile,
  sep = sep,
  row.names = FALSE,
  col.names = TRUE
)
write.table(
  do.call(rbind, res$products),
  file = productFile,
  sep = sep,
  row.names = FALSE,
  col.names = TRUE
)

if (!is.null(rdsFile)) 
  saveRDS(res, rdsFile)
cat("OK!\n")