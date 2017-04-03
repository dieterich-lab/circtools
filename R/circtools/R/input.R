
getSjExons <- function(txdb, circsGR) {
  extxdb <- exons(txdb, columns = c("EXONNAME", "GENEID", "TXNAME"))
  circExonsMap <- list(# t because we want it ordered as circs
    rightSide = t(findOverlaps(extxdb, circsGR, type = "start")),
    leftSide = t(findOverlaps(extxdb, circsGR, type = "end")))
  sjExonCoords <- lapply(circExonsMap, function(x) extxdb[x@to])
  geneKey <- "GENEID"
  genes <- unlist(mcols(sjExonCoords$rightSide)[[geneKey]])
  sjExByGene <- lapply(sjExonCoords, function(x) split(x, genes))
  sjExByGene
}

getTxByGenes <- function(txdb, genes) {
  suppressMessages(
    txByGeneMap <- mapIds(
      txdb,
      column = "TXNAME",
      key = genes,
      keytype = "GENEID",
      multiVals = "CharacterList"
    )
  )
  txByGeneMap
}

getIntersectingTx <- function(exByGene, type=c("all", "<", ">")) {
  #TODO: different output types
  intersectingTx <-  lapply(
    names(exByGene$leftSide),
    function(gene) {
      intersect(unlist(mcols(exByGene[[1]][[gene]])$TXNAME),
                unlist(mcols(exByGene[[2]][[gene]])$TXNAME))
    })
  names(intersectingTx) <- names(exByGene$leftSide)
  intersectingTx
}

#' @import GenomicFeatures
getTx <- function(txdb, genes) {
  txCoords <-exons(txdb,
                          columns = c("GENEID", "TXNAME", "EXONNAME"),
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

#' Create circular RNA coordinates GRanges object
#'
#' @param chr 
#' @param start 
#' @param end 
#' @param strand 
#' @param ids 
#'
#' @return a GRanges object
#' @export
#' @import GenomicRanges IRanges
#' @examples
getCircCoords <- function(chr, start, end, strand, ids = NULL) {  
  circIR <- IRanges(start = start, 
                    end = end,
                    names =  chr)
  grCircs <- GRanges(seqnames = chr,
                     ranges = circIR,
                     strand = strand)
  if (is.null(ids)) {
    mcols(grCircs)$CIRCID <- apply(
      cbind(chr, start, end),
      FUN = paste,
      MARGIN = 1,
      collapse = '-'
    )
  } else {
    mcols(grCircs)$CIRCID <- ids
  }
  grCircs
}

getGeneIds <- function() {
  
}

