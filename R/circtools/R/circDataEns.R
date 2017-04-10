
#' Get coordinates of splice junction exons
#'
#' @param dban ensembldb object
#'
#' @param circsGR a GRanges of circular splice junctions. 
#'   Must have `CIRCID` column with the identifiers of the circular splice 
#'   juctions.
#' @param filter a filter to use for exons output
#' 
#' @return a list with `leftSide` and `rightSide` GRanges of the 
#' intersecting exons.
#'
#' @importFrom ensembldb exons
#' @importFrom GenomicRanges findOverlaps mcols mcols<-
#' @importFrom S4Vectors to from 
getSjExons <- function(db, circsGR, filter=list()) {
  exdb <- exons(db, filter = filter)
  circExonsMap <- list(# t because we want it ordered as circs
    rightSide = t(findOverlaps(exdb, circsGR, type = "start")),
    leftSide  = t(findOverlaps(exdb, circsGR, type = "end")))
  sjExonCoords <- lapply(circExonsMap,
                         function(x) {
                           res <- exdb[to(x)]
                           mcols(res)$CIRCID <- mcols(circsGR)$CIRCID[from(x)]
                           res
                         })
  sjExonCoords
}


#' Create a CircData object
#'
#' @param db 
#' @param circCoords 
#'
#' @return a CircData object
#' @export
#'
#' @importFrom GenomicRanges findOverlaps mcols
#' @importFrom ensembldb genes GRangesFilter
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
  class(circData) <- "CircData"
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
#' @param primers a data.frame or an IRanges object for the primers
#' @param counts a data.frame with an id column (tx_id) and corresponding
#'   read counts.
#'
#' @export
#'
plotCirc <- function(circIds,
                     circGenes,
                     circData,
                     primers = NULL,
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
  plotTranscripts(ex, circs = circs, primers = primers, counts = counts)
}

#' Retrieve sequences of the exons, which have common start or end with 
#' the splice junction coordinates provided in circData
#'
#' @param circData a `CircData` object
#' @param bsg a `BSgenome` entity
#'
#' @return a list of `leftSide` and `rightSide` exons
#' @export
#'
#' @examples
#' \dontrun{
#' 
#' }
getExonSeqs <- function(circData, bsg) {
  ex <- getSjExons(db = circData$db, 
                   circsGR = circData$circCoords,
                   filter = GeneidFilter(circData$sjGeneIds))
  for (side in names(ex)) {
    GenomicRanges::mcols(ex[[side]])$seq <- BSgenome::getSeq(x = bsg,
                                                             names = ex[[side]])
  }
  ex
}