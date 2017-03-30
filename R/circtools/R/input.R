
getExons <- function(txdb, circsGR) {
  extxdb <- exons(txdb, columns = c("EXONNAME", "GENEID", "TXNAME"))
  circExonsMap <- list(# t because we want it ordered as circs
    rightSide = t(findOverlaps(extxdb, circsGR, type = "start")),
    leftSide = t(findOverlaps(extxdb, circsGR, type = "end")))
  sjExonCoords <- lapply(circExonsMap, function(x) extxdb[x@to])
  geneKey <- "GENEID"
  genes <- unlist(mcols(sjExonCoords$rightSide)[[geneKey]])
  exByGene <- lapply(sjExonCoords, function(x) split(x, genes))
  exByGene
}

getIntersectingTx <- function(exByGene) {
  intersectingTx <-  lapply(
    names(exByGene$leftSide),
    function(gene) {
      intersect(unlist(mcols(exByGene[[1]][[gene]])$TXNAME),
                unlist(mcols(exByGene[[2]][[gene]])$TXNAME))
    })
  names(intersectingTx) <- names(exByGene$leftSide)
  intersectingTx
}

getTx <- function(txdb, genes) {
  txCoords <- transcripts(txdb,
                          columns = c("GENEID", "TXNAME"),
                          filter = list(gene_id = genes))
  txCoordsByGene <- split(txCoords, unlist(mcols(txCoords)$GENEID))
  txCoordsByGene
}

getExonsSeq <- function() {
  
}

getCounts <- function() {
  
}

getPrimersCoord <- function() {
  
}

getCircCoords <- function(table) {
  circIR <- IRanges(start = table[, 2],
                    end = table[, 3],
                    names = table[, 1])
  grCircs <- GRanges(seqnames = table[, 1],
                     ranges = circIR,
                     strand = table[, 4])
  mcols(grCircs)$CIRCID <- apply(table[, 1:3],
                                 FUN = paste,
                                 MARGIN = 1,
                                 collapse = '-')
  grCircs
}

getGeneIds <- function() {
  
}

library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb2 <- TxDb.Hsapiens.UCSC.hg19.knownGene
annotationDir <- "/home/alex/prj/mcf7_rna_metabolism/annotation/"
txdb <- makeTxDbFromGFF(paste0(annotationDir,"GRCh38.84.gtf"))
saveDb(txdb, paste0(annotationDir, "txdb"))
circsCoords <- read.table(
  "/home/alex/prj/mcf7_rna_metabolism/data/CircCoordinates",
  header = TRUE,
  sep = "\t",
  nrows = 20
)
grCircs <- getCircCoords(circsCoords[, c("Chr", "Start", "End", "Strand")])
exByGene <- getExons(txdb,grCircs)
intersectingTx <- getIntersectingTx(exByGene)
txByGene <- getTx(txdb, names(intersectingTx))

