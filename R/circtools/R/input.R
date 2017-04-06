#' Create lists for exons around the splice-junctions
#'
#' @param db a TxDB or ensembldb object 
#'
#' @param circsGR a GRanges object with the splice junction coordinates
#' 
#' @details The db object must have a key "GENEID" and columns
#' "EXONID" and "TXNAME", which are used for creation of the output object.
#' 
#' @return a list of two GRangesLists (by gene) objects: 
#' for the left and right side splice  junctions exons
#'
#' @importFrom GenomicFeatures exons
#' @importFrom GenomicRanges findOverlaps
getSjExons <- function(db, circsGR) {
  exdb <- exons(db, columns = c("EXONID", "GENEID", "TXNAME"))
  circExonsMap <- list(# t because we want it ordered as circs
    rightSide = t(findOverlaps(exdb, circsGR, type = "start")),
    leftSide = t(findOverlaps(exdb, circsGR, type = "end")))
  sjExonCoords <- lapply(circExonsMap, function(x) exdb[x@to])
  geneKey <- "GENEID"
  genes <- unlist(mcols(sjExonCoords$rightSide)[[geneKey]])
  sjExByGene <- lapply(sjExonCoords, function(x) split(x, genes))
  sjExByGene
}

getTxByGenes <- function(db, genes) {
  suppressMessages(
    txByGeneMap <- mapIds(
      db,
      column = "TXNAME",
      key = genes,
      keytype = "GENEID",
      multiVals = "CharacterList"
    )
  )
  txByGeneMap
}

#' @importFrom BiocGenerics intersect
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

getTx <- function(txdb, genes) {
  txCoords <- GenomicFeatures::exons(txdb,
                          columns = c("GENEID", "TXNAME", "EXONID"),
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
#' 
#' @export
#' @examples
getCircCoords <- function(chr, start, end, strand, ids = NULL) {  
  circIR <- IRanges::IRanges(start = start, 
                    end = end,
                    names =  chr)
  grCircs <- GenomicRanges::GRanges(seqnames = chr,
                     ranges = circIR,
                     strand = strand)
  if (is.null(ids)) {
    S4Vectors::mcols(grCircs)$CIRCID <- apply(
      cbind(chr, start, end),
      FUN = paste,
      MARGIN = 1,
      collapse = '-'
    )
  } else {
    S4Vectors::mcols(grCircs)$CIRCID <- ids
  }
  grCircs
}

getGeneIds <- function() {
  
}

