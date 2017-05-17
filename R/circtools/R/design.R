#' Generate sequencing of splice junction regions
#'
#' @param exSeq a GRanges object with the metadata columns:
#'   - exon_id
#'   - gene_id
#'   - sjId (splice junction id)
#'   - side ('left', 'right')
#'   - seq, a DNAStringSet
#'
#' @return a data.frame with sjId in rownames and circSeq column with
#' the splice junction sequences
#' @noRd
#'
.getCircSJSeq <- function(exSeq) {
  df <- S4Vectors::mcols(exSeq)
  sjId <- unique(S4Vectors::mcols(exSeq)$sjId)
  strand <- as.character(strand(exSeq))
  stream <- .getStream(strand, df$side)
  indByStream <- expand.grid(split(seq_along(stream), stream))
  res <-  lapply(seq_along(indByStream$downstream),
                 function(i) {
                   up <- indByStream$upstream[i]
                   down <- indByStream$downstream[i]
                   c(
                     sjId = sjId,
                     upExonId = df$exon_id[up],
                     downExonId = df$exon_id[down],
                     circSeq = paste0(df$seq[up], df$seq[down])
                   )
                 })
  as.data.frame(do.call(rbind, res), stringsAsFactors = FALSE)
}

#' Generate splice junction sequencies for a list of circs
#'
#' @param exSeqList a list of exSeq items, defined in `\link{.getCircSJSeq}`
#'
#' @return a list
#' @noRd
#'
getCircSeqFromList <- function(exSeqList) {
  seqs <- do.call(rbind, lapply(exSeqList, .getCircSJSeq))
  seqs$seqId <- row.names(seqs)
  rownames(seqs) <- NULL
  seqs
}

#' Get transcript ids which relate to the circs
#'
#' @param db ensembldb object
#' @param exSeq result of getExonSeqs
#' @param whichExons a character:
#'   - 'any': all transcripts which intersect whether with the left or right
#'       splice junction exons
#'   - 'both': only those, which intersect with the both, left and right exons.
#'
#' @return a list of transcripts ids by circ
#'
getTxOfSJExons <- function(db, exSeq, whichExons=c('any', 'both')){
  whichExons <- match.arg(whichExons)
  lapply(exSeq, function(x) {
    sides <- split(x, mcols(x)$side)
    txBySide <- lapply(sides, function(e)
      mcols(ensembldb::transcriptsByOverlaps(db, e))$tx_id)
    res <- switch(
      whichExons,
      any = union(txBySide$left, txBySide$right),
      both = intersect(txBySide$left, txBySide$right)
    )
    res
  })
}

#' Retirieve transcript seqs for the genes of circs
#'
#' @param db ensembldb object
#' @param bsg BSGenome or FaFile object 
#'   (see `x` argument in \code{\link{extractTranscriptSeqs}})
#' @param exSeq the result of \code{\link{getExonSeqs}} function
#'
#' @return a list by transcripts with their sequences (`\link[Biostrings]{DNAString}`)
#'
getTxSeqs <- function(db, bsg, exSeq) {
  geneIds <- vapply(exSeq, function(x) unique(
    S4Vectors::mcols(x)$gene_id), character(1))
  ex <- ensembldb::exonsBy(db, filter = AnnotationFilter::GeneIdFilter(geneIds))
  seqs <- GenomicFeatures::extractTranscriptSeqs(bsg, ex)
  seqs
}

#' Design primers
#'
#' @param exSeq a GRanges list with GRanges of exons on the splice junction.
#'   Every item represents an individual splice junction.
#' @param db an ensembldb object
#' @param bsg a BSgenome object
#' @param opts a named list with the design options (see details).
#'
#' @return a list with `primers` and `product` items.
#' @details  
#' 
#' The following options are available for the primer design:
#'   - minLength (0), maxLength (Inf): a min and max product length
#' @export
#'
designPrimers <- function(exSeq,
                          db,
                          bsg,
                          opts = list(minLength = 0, maxLength = Inf)) {
  # TODO choice if intersect with SJ
  seqs <- getCircSeqFromList(exSeq)
  txIntersect <- getTxOfSJExons(db, exSeq, whichExons = 'both')
  txSeqs <- getTxSeqs(db, bsg, exSeq)
  
  # lapply per sjId
  res <- lapply(names(exSeq), function(sjId) {
    ts <- txSeqs[txIntersect[[sjId]]]
    .designForSJ(seqs[seqs$sjId == sjId,], exSeq[[sjId]],
                   ts, opts = opts)
  })
  names(res) <- names(exSeq)
  list(primers = res,
       products = split(seqs, seqs$sjId))
}

.designForSJ <- function(circSeqs,
                         ex,
                         ts,
                         opts = list()) {
  # add sequences to the db and corresponding tx's
  dbConn <- RSQLite::dbConnect(RSQLite::SQLite(), ":memory:")
  stringSet <- Biostrings::DNAStringSet(circSeqs$circSeq)
  names(stringSet) <- circSeqs$seqId
  suppressWarnings({
    DECIPHER::Seqs2DB(
      c(stringSet, ts),
      "XStringSet",
      dbConn,
      identifier = c(circSeqs$seqId, names(ts)),
      verbose = FALSE
    )
  })
  # generate sequence tiles for DesignPrimers function
  tiles <- DECIPHER::TileSeqs(dbConn,
                              add2tbl = "Tiles",
                              minCoverage = 1,
                              verbose = FALSE)
  # TODO: suppress msgs
  primers <- lapply(circSeqs$seqId, function(seqId) {
    primers <- DECIPHER::DesignPrimers(
      tiles = tiles,
      identifier = seqId,
      minCoverage = 1,
      minGroupCoverage = 1,
      #numPrimerSets = 5,
      maxSearchSize = 1e4,
      worstScore = -10,
      #minProductSize = 60,
      #minEfficiency = .8,
      verbose = FALSE
    )
  })
  primers <- decipher2iranges(primers)
  names(primers) <- circSeqs$seqId
  lengthRange <- c(opts$minLength, opts$maxLength)
  bestPrimers <- lapply(names(primers), function(sjId) {
    sj <- IRanges::IRanges(
      Biostrings::width(ex$seq[ex$exon_id == circSeqs$upExonId]), width = 2)
    selectBestPrimers(p=primers[[sjId]], sj, lengthRange)
  })
  names(bestPrimers) <- names(primers)
  primersGR <- lapply(names(bestPrimers), function(seqId) {
    i <- which(seqId == circSeqs$seqId)
    upExon   <- ex[which(S4Vectors::mcols(ex)$exon_id == circSeqs$upExonId[i])]
    downExon <- ex[which(S4Vectors::mcols(ex)$exon_id == circSeqs$downExonId[i])]
    circ2genome(bestPrimers[[seqId]],
                upExon = upExon,
                downExon = downExon)
  })
  names(primersGR) <- names(bestPrimers)
  GenomicRanges::GRangesList(primersGR)
}

  

selectBestPrimers <- function(p, sj, lengthRange) {
  sjPrimers <- p[from(findOverlaps(p, sj))]
  sjPrimers <- split(sjPrimers, S4Vectors::mcols(sjPrimers)$type)
  best <- lapply(sjPrimers, function(x) {
    x[which.max(S4Vectors::mcols(x)$efficiency)]
  })
  p <- split(p, S4Vectors::mcols(p)$type)
  # forward primer on sj
  rangeForEnd <- IRanges::start(best$forward) + lengthRange - 1
  o <- IRanges::end(p$reverse) >= rangeForEnd[1] &
       IRanges::end(p$reverse) <= rangeForEnd[2]
  result <- list()
  if (any(o)) {
    result$forwardSJ <- c(best$forward,
                          p$reverse[which.max(mcols(p$reverse[o])$efficiency)])
  }
  # reverse primer on sj
  rangeForStart <- IRanges::end(best$reverse) - lengthRange + 1
  o <- IRanges::start(p$forward) >= rangeForStart[1] &
       IRanges::start(p$forward) <= rangeForStart[2]
  if (any(o)) {
    result$reverseSJ <- c(p$forward[which.max(mcols(p$forward[o])$efficiency)],
        best$reverse)
  }
  bestInd <- which.max(lapply(result, function(x) {
    sum(S4Vectors::mcols(x)$efficiency)
  }))
  result[[bestInd]]
}

#' Creates  IRanges objects 
#' 
#'  The  metadata include:
#'   - seqId
#'   - productSize
#'   - type: \['forward', 'reverse'\]
#'   - seq: the primer sequence.
#'   
#' @param primers is a list of results from 
#' \code{\link[DECIPHER]{DesignPrimers}}
#' 
#' @return  a list with an \code{\link[IRanges]{IRanges}}
#'  item for every primer pair
decipher2iranges <- function(primers) {
  res <- lapply(primers, function(p) {
    fwSelect <- p$score_forward > -Inf
    fw <- IRanges::IRanges(start = p$start_forward[fwSelect],
                  width = nchar(p$forward_primer[fwSelect,1]))
    S4Vectors::mcols(fw)$seq <- p$forward_primer[fwSelect,1]
    S4Vectors::mcols(fw)$type <- 'forward'
    S4Vectors::mcols(fw)$efficiency <- p$forward_efficiency[fwSelect,1]
    rvSelect <- p$score_reverse > -Inf
    rv <- IRanges::IRanges(start = p$start_reverse[rvSelect],
                  width = nchar(p$reverse_primer[rvSelect,1]))
    S4Vectors::mcols(rv)$seq <- p$reverse_primer[rvSelect,1]
    S4Vectors::mcols(rv)$type <- 'reverse'
    S4Vectors::mcols(rv)$efficiency <- p$reverse_efficiency[rvSelect,1]
    both <- c(fw,rv)
    S4Vectors::mcols(both)$seqId <- p$identifier[1]
    both
  })
  #names(res) <- primers$identifier
  res
}

.toGenome <- function(x, circStrand, upExon, downExon) {
  vapply(x, function(y) {
    if (circStrand == '+') {
      if (y <= IRanges::width(upExon)) {
        res <- y + IRanges::start(upExon) - 1
      } else {
        res <- y - IRanges::width(upExon) + IRanges::start(downExon) - 1
      }
      res
    } else {
      if (y <= IRanges::width(upExon)) {
        res <- IRanges::end(upExon) - y + 1
      } else {
        res <-  IRanges::end(downExon) - (y - IRanges::width(upExon)) + 1
      }
      res
    }
  }, numeric(1))
}

#' Transform from intra-circ coordinates to the genome ones
#' 
#' Works only within two adjacent splice junction exons termed
#' upstream and downatream ones.
#'
#' @param x a IRanges object
#' @param upExon a GRanges object
#' @param downExon a GRanges object
#'
#' @return GRangesList
#'
circ2genome <- function(x, upExon, downExon) { 
  circStrand <- unique(GenomicRanges::strand(upExon))
  # transform to the genome coords and order as c(min, max) using `range`
  res <-   lapply(
    list(start = IRanges::start(x), end = IRanges::end(x)),
    .toGenome,
    upExon = upExon,
    downExon = downExon,
    circStrand = circStrand
  )
  res <- do.call(rbind, Map(range, res$start, res$end))
  res <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(
      unique(GenomicRanges::seqnames(upExon)), length(x)),
    IRanges::IRanges(start = res[, 1], end = res[, 2]),
    strand = S4Vectors::Rle(unique(GenomicRanges::strand(upExon)), length(x))
  )
  S4Vectors::mcols(res) <- S4Vectors::mcols(x)
  res <- lapply(res, splitPrimer, upExon = upExon, downExon = downExon)
  do.call(c, res)
}

# We assume that primer length is > 2
splitPrimer <- function(x, upExon, downExon) {
  if (IRanges::overlapsAny(x, upExon) && IRanges::overlapsAny(x, downExon)) {
    x <- GenomicRanges::narrow(x, start = 2, end = -2)
    sp <- GenomicRanges::setdiff(c(upExon, downExon),x)
    mcols(sp) <- mcols(x)
    sp
  } else {
    x
  }
}

# get exons which intersect with the designed primers
# 
exons4primers <- function(db, primer) {
  conditions <- c('within', 'overlapping')
  ex <- lapply(conditions, function(cc) {
    ensembldb::exons(
      db,
      filter = AnnotationFilter::GRangesFilter(primer, condition = cc),
      return.type = 'data.frame',
      columns = 'exon_id'
    )
  })
  names(ex) <- conditions
}