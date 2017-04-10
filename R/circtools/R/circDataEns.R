
#' @importFrom ensembldb exons
#' @importFrom GenomicRanges findOverlaps mcols mcols<-
#' @importFrom S4Vectors to from 
getSjExons <- function(db, circsGR) {
  exdb <- exons(db)
  circExonsMap <- list(# t because we want it ordered as circs
    rightSide = t(findOverlaps(exdb, circsGR, type = "start")),
    leftSide = t(findOverlaps(exdb, circsGR, type = "end")))
  sjExonCoords <- lapply(circExonsMap,
                         function(x) {
                           res <- exdb[to(x)]
                           mcols(res)$CIRCID <-
                             mcols(circsGR)$CIRCID[from(x)]
                           res
                         })
  sjExonCoords
}


#' Create a CircData object
#'
#' @param db 
#' @param circCoords 
#'
#' @return
#' @export
#'
#' @importFrom GenomicRanges findOverlaps mcols
#' @importFrom ensembldb genes GRangesFilter
#' @examples
CircData <- function(db, circCoords) {
  # why does not work with a list of two filters?
  sjFilter <- GRangesFilter(circCoords, "overlapping")
  sjGenes <- genes(db, filter = sjFilter)
  circsGeneHits <- findOverlaps(circCoords, sjGenes)
  
  circData <- list(
    db = db,
    sjGeneIds = mcols(sjGenes)$gene_id,
    circsGeneHits = circsGeneHits,
    circCoords = circCoords,
    sjGenes = sjGenes
  )
  circData
}

#' Plot a transcript scheme 
#'
#' @param circIds a character identifiers for the circular transcript 
#'   to be plotted.
#' @param circGene a character identifier for a gene, which transcripts are to 
#'   be plotted. If `circIds` is not set, all circular transcripts defined in
#'   `circData` will be used. The identifier is identical to the GENEID
#'   record in the EnsDb or TxDB objects.
#' @param circData an object returned by CircData
#'
#' @export
#'
plotCirc <- function(circIds,
                     circGenes,
                     circData,
                     counts = NULL) {
  if (!missing(circGenes)) {
    if (length(circGenes) > 1)
      warning(paste(
        "More than one gene id provided: the transcripts from the genes ",
              " will be plotted together."))
    if (!all(circGenes %in% circData$sjGeneIds))
      stop(paste( "Error: The provided gene id ",
          circGenes, " does not overlap with the splice junction. "))
    ind <- which(circData$sjGeneIds == circGenes)
    allCircInd <-  with(circData, from(circsGeneHits)[to(circsGeneHits) == ind])
    allCircId <- mcols(circData$circCoords)$CIRCID[allCircInd]
    if (missing(circIds)) {
      circIds <- allCircId
    } else {
      if (!all(circIds %in% allCircId))
        stop("Some provided circs do not overlap with the genes in circGenes")
    }
  } else {
    circIndex <- which(mcols(circData$circCoords)$CIRCID == circIds)
    circGenes <- with(
      circData,
      sjGenes[to(circsGeneHits)[from(circsGeneHits) == circIndex]])
    circGenes <- mcols(circGenes)$gene_id
  }
  circIndex <- which(mcols(circData$circCoords)$CIRCID == circIds)
  circs <- circData$circCoords[circIndex]
  ex <- exons(
    circData$db,
    filter = GeneidFilter(circGenes),
    columns = c("gene_id", "tx_id", "tx_name")
  )
  plotTranscripts(ex, circs = circs,counts = counts)
}

