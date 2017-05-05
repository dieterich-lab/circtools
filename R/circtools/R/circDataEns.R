
#' Get coordinates of splice junction exons
#'
#' @param db an ensembldb object
#'
#' @param circsGR a GRanges of circular splice junctions. 
#'   Must have `sjId` column with the identifiers of the circular splice 
#'   juctions.
#' @param filter a filter to use for exons output
#' 
#' @return a list with `left` and `right` GRanges of the 
#' intersecting exons.
#'
#' @importFrom ensembldb exons
#' @importFrom GenomicRanges findOverlaps mcols mcols<-
#' @importFrom S4Vectors to from 
getSjExons <- function(db, circsGR, filter=list()) {
  exdb <- exons(db, filter = filter)
  circExonsMap <- list(# t because we want it ordered as circs
    right = t(findOverlaps(exdb, circsGR, type = "start")),
    left  = t(findOverlaps(exdb, circsGR, type = "end")))
  sjExonCoords <- lapply(
    names(circExonsMap),
    function(side) {
      res <- exdb[to(circExonsMap[[side]])]
      mcols(res)$sjId <- mcols(circsGR)$sjId[from(circExonsMap[[side]])]
      mcols(res)$side <- side
      res
    })
  do.call(c, sjExonCoords)
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
#' @param sjIds a character identifiers for the circular transcript 
#'   to be plotted.
#' @param circGene a character identifier for a gene, which transcripts are to 
#'   be plotted. If `sjIds` is not set, all circular transcripts defined in
#'   `circData` will be used. The identifier is identical to the GENEID
#'   record in the EnsDb or TxDB objects.
#' @param circData an object returned by CircData
#' @param primers a data.frame or an IRanges object for the primers
#' @param counts a data.frame with an id column (tx_id) and corresponding
#'   read counts.
#'
#' @export
#'
plotCirc <- function(sjIds,
                     circGenes,
                     circData,
                     primers = NULL,
                     counts = NULL, 
                     opts = NULL) {
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
    allCircId <- mcols(circData$circCoords)$sjId[allCircInd]
    if (missing(sjIds)) {
      sjIds <- allCircId
    } else {
      if (!all(sjIds %in% allCircId))
        stop("Some provided circs do not overlap with the genes in circGenes")
    }
  } else {
    circIndex <- which(mcols(circData$circCoords)$sjId == sjIds)
    circGenes <- with(
      circData,
      sjGenes[to(circsGeneHits)[from(circsGeneHits) == circIndex]])
    circGenes <- mcols(circGenes)$gene_id
  }
  circIndex <- which(mcols(circData$circCoords)$sjId == sjIds)
  circs <- circData$circCoords[circIndex]
  ex <- exons(
    circData$db,
    filter = GeneidFilter(circGenes),
    columns = c("gene_id", "tx_id", "tx_name")
  )
  plotTranscripts(
    ex,
    circs = circs,
    primers = primers,
    counts = counts,
    opts = opts
  )
}

#' Retrieve sequences of the exons, which have common start or end with 
#' the splice junction coordinates provided in circData
#'
#' @param circData a `CircData` object
#' @param bsg a `BSgenome` entity
#' @param faFile a \code{\link{FaFile}} object
#' @param type 
#'
#' @return GRangesList object with a record for every circular id
#' defined in the `circData` argument. The metadata in the list items 
#' include:
#' * `exon_id` and `gene_id`
#' * `seq`, a DNAStringSet object
#' * `side`, a character string ("left" or "right")
#' * `sjId`, a character string derived from the circData object
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' 
#' }
getExonSeqs <- function(circData, bsg, faFile, 
                        type = c('all', 'shortest', 'longest')) {
  type <- match.arg(type)
  if (!missing(faFile) && !missing(bsg))
      stop("Specify only one sequence source: BSgenome or FaFile")
  ex <- getSjExons(db = circData$db, 
                   circsGR = circData$circCoords,
                   filter = GeneidFilter(circData$sjGeneIds))
  byCirc <- GenomicRanges::split(ex, GenomicRanges::mcols(ex)$sjId)
  if (type != "all") {
    fun <- switch(
      type,
      shortest = function(x) x[which.min(width(x))],
      longest  = function(x) x[which.max(width(x))]
    )
    byCirc <- S4Vectors::endoapply(byCirc, function(gr) {
      bySide <- GenomicRanges::split(gr, gr$side)
      unlist(S4Vectors::endoapply(bySide, fun))
    })
  }
  if (!missing(bsg)) {
    S4Vectors::endoapply(byCirc, function(gr) {
      names(gr) <- NULL
      GenomicRanges::mcols(gr)$seq <- BSgenome::getSeq(x = bsg, names = gr)
      gr
    })
  } else {
    S4Vectors::endoapply(byCirc, function(gr) {
      names(gr) <- NULL
      GenomicRanges::mcols(gr)$seq <- Rsamtools::getSeq(x = faFile, gr)
      gr
    })
  }
}
