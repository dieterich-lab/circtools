
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
#' intersecting exons. 'left' means that an exon GRanges start coordinate is 
#' the same as the SJ start, 'right' is the same, but for exon GRanges ends.
#'
#' @importFrom ensembldb exons
#' @importFrom GenomicRanges findOverlaps mcols mcols<-
#' @importFrom S4Vectors to from 
getSjExons <- function(db, circsGR, filter=list()) {
  exdb <- exons(db, filter = filter)
  circExonsMap <- list(# t because we want it ordered as circs
    left = t(findOverlaps(exdb, circsGR, type = "start")),
    right = t(findOverlaps(exdb, circsGR, type = "end")))
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
#' @param db an ensembldb object
#' @param circCoords an IRanges object
#'
#' @return a CircData object
#' @export
#'
#' @importFrom GenomicRanges findOverlaps mcols
#' @importFrom ensembldb genes 
CircData <- function(db, circCoords) {
  if (requireNamespace("AnnotationFilter", quietly = TRUE)) {
    sjFilter <- AnnotationFilter::GRangesFilter(circCoords, "overlapping")
  } else {
    sjFilter <- ensembldb::GRangesFilter(circCoords, condition="overlapping")
  }
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
#' @param circData an object returned by CircData
#' @param primers a data.frame or an IRanges object for the primers
#' @param counts a data.frame with an id column (tx_id) and corresponding
#'   read counts.
#' @param circGenes a character vector or list with the gene id's, for which
#' circular transcripts must be plotted
#'   
#' @param opts a list with options for the \code{\link{plotTranscripts}} 
#' opts argument.
#' @param countThres the transcripts with the read counts exceeding the 
#' threshold value will be plotted (default: 0)
#' 
#' @details The default options in the list `opts` are the following:
#' ```
#'    normalise   = TRUE
#'    net         = TRUE
#'    primerColor = "firebrick3"
#'    exonColor   = "deepskyblue1"
#'    netColor    = "grey"
#' ```  
#' - normalise: if the coordinates should be transformed to pseudoposition 
#' for better representation
#' - net: if `normalise = TRUE`, add vertical lines to track exon and primer
#' positions
#' @export
#'
plotCirc <- function(sjIds,
                     circGenes,
                     circData,
                     primers = NULL,
                     counts = NULL, 
                     countThres = -Inf,
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
  if (requireNamespace("AnnotationFilter", quietly = TRUE)) {
    filter <- AnnotationFilter::GeneIdFilter(circGenes)
  } else {
    filter <- ensembldb::GeneidFilter(circGenes)
  }
  ex <- exons(
    circData$db,
    filter = filter, 
    columns = c("gene_id", "tx_id", "tx_name")
  )
  if (!is.null(primers)) {
    primers <-  as.data.frame(primers, row.names = NULL)
    primers$type[primers$type == 'forward'] <- 'FW'
    primers$type[primers$type == 'reverse'] <- 'RV'
    primers$id <- paste(primers$seqId, primers$type) 
  }
  if (!is.null(counts)) {
   counts <- droplevels(counts[counts$count >= countThres, ])
   ex <- ex[mcols(ex)$tx_id %in% counts$id]
  }
  plotTranscripts(
    exons   = droplevels(as.data.frame(ex, row.names = NULL)),
    circs   = as.data.frame(circs, row.names = NULL),
    primers = primers,
    counts  = counts,
    opts    = opts
  )
}

#' Retrieve sequences of the exons, which have common start or end with 
#' the splice junction coordinates provided in circData
#'
#' @param circData a `CircData` object
#' @param bsg a `BSgenome` entity
#' @param faFile a \code{\link{FaFile}} object
#' @param type a character: whether to include all, only the shortest or 
#'   the longest exons, which can form the splice junction
#'   
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
getExonSeqs <- function(circData, bsg, faFile, 
                        type = c('all', 'shortest', 'longest')) {
  type <- match.arg(type)
  if (!missing(faFile) && !missing(bsg))
      stop("Specify only one sequence source: BSgenome or FaFile")
      
  if (requireNamespace("AnnotationFilter", quietly = TRUE)) {
    filter <- AnnotationFilter::GeneIdFilter(circData$sjGeneIds)
  } else {
    filter <- ensembldb::GeneidFilter(circData$sjGeneIds)
  }
  ex <- getSjExons(db = circData$db, 
                   circsGR = circData$circCoords,
                   filter = filter)
  byCirc <- GenomicRanges::split(ex, GenomicRanges::mcols(ex)$sjId)
  if (type != "all") {
    fun <- switch(
      type,
      shortest = function(x) x[which.min(GenomicRanges::width(x))],
      longest  = function(x) x[which.max(GenomicRanges::width(x))]
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
