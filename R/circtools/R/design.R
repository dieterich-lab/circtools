#' Generate sequencing of splice junction regions
#'
#' @param exSeq 
#'
#' @return a data.frame with CIRCID in rownames and circSeq column with
#' the splice junction sequences
#'
.getCircSJSeq <- function(exSeq) {
  df <- mcols(exSeq)
  circId <- unique(mcols(exSeq)$CIRCID)
  strand <- as.character(strand(exSeq))
  stream <- .getStream(strand, df$side)
  indByStream <- expand.grid(split(seq_along(stream), stream))
  res <-  lapply(seq_along(indByStream$downstream),
                 function(i) {
                   up <- indByStream$upstream[i]
                   down <- indByStream$downstream[i]
                   c(
                     CIRCID = circId,
                     upExonId = df$exon_id[up],
                     downExonId = df$exon_id[down],
                     circSeq = paste0(df$seq[down], df$seq[up])
                   )
                 })
  as.data.frame(do.call(rbind, res), stringsAsFactors = FALSE)
}

#' Generate splice junction sequencies for a list of circs
#'
#' @param exSeqList 
#'
#' @return
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
#' @param filter a character:
#'   - 'any': all transcripts which intersect whether with the left or right
#'       splice junction exons
#'   - 'both': only those, which intersect with the both, left and right exons.
#'
#' @return a list of transcripts ids by circ
#'
getTxForCircs <- function(db, exSeq, filter=c('any', 'both')){
  filter <- match.arg(filter)
  lapply(exSeq, function(x) {
    sides <- split(x, mcols(x)$side)
    txBySide <- lapply(sides, function(e)
      mcols(transcriptsByOverlaps(db, e))$tx_id)
    res <- switch(
      filter,
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
#' @return a list by transcripts with their sequences (`\link{DNAString}`)
#'
getTxSeqs <- function(db, bsg, exSeq) {
  geneIds <- vapply(exSeq, function(x) unique(mcols(x)$gene_id), character(1))
  ex <- exonsBy(db, filter = GeneidFilter(geneIds))
  seqs <- extractTranscriptSeqs(bsg, ex)
  seqs
}

designPrimers <- function(exSeq, db, bsg) {
  seqs <- getCircSeqFromList(exSeq)
  txIntersect <- getTxForCircs(db, exSeq, 'both')
  txSeqs <- getTxSeqs(db, bsg, exSeq)
  
  # lapply per circId
  res <- lapply(names(exSeq), function(circId) {
    ts <- txSeqs[txIntersect[[circId]]]
    .designForCirc(seqs[seqs$CIRCID == circId,], exSeq[[circId]], ts)
  })
  names(res) <- names(exSeq)
  list(primers = res,
       products = seqs)
}

.designForCirc <- function(circSeqs, ex, ts) {
  dbConn <- RSQLite::dbConnect(SQLite(), ":memory:")
  suppressWarnings({
    DECIPHER::Seqs2DB(
      DNAStringSet(circSeqs$circSeq),
      "XStringSet",
      dbConn,
      identifier = circSeqs$seqId,
      verbose = FALSE
    )
    # add corresponding tx's
    DECIPHER::Seqs2DB(
      ts, "XStringSet", dbConn, identifier = names(ts), verbose = FALSE)
  })
  tiles <- DECIPHER::TileSeqs(dbConn,
                              add2tbl = "Tiles",
                              minCoverage = 1,
                              verbose = FALSE)
  # design for every seqId
  # TODO: suppress msgs
  primers <- lapply(circSeqs$seqId, function(seqId) {
    primers <- DesignPrimers(
      tiles = tiles,
      identifier = seqId,
      minCoverage = 1,
      minGroupCoverage = 1,
      numPrimerSets = 1,
      maxSearchSize = 20,
      minProductSize = 60,
      minEfficiency = .9,
      verbose = FALSE
    )
  })
  primers <- decipher2iranges(primers)
  names(primers) <- circSeqs$seqId
  
  primersGR <- lapply(names(primers), function(seqId) {
    i <- which(seqId == circSeqs$seqId)
    upExon   <- ex[which(mcols(ex)$exon_id == circSeqs$upExonId[i])]
    downExon <- ex[which(mcols(ex)$exon_id == circSeqs$downExonId[i])]
    circ2genome(primers[[seqId]],
                upExon = upExon,
                downExon = downExon)
  })
  names(primersGR) <- names(primers)
  primersGR
}

#' Creates  IRanges objects with metadata on 
#'   - seqId
#'   - productSize
#'   - type: ['forward', 'reverse']
#'   - seq: the primer sequence.
#'   @param primers is a list of results from \link[DECIPHER]{DesignPrimers}
#'   @return  a list with an IRanges item for every primer pair
decipher2iranges <- function(primers) {
  res <- lapply(primers, function(p) {
    fw <- IRanges(start = p$start_forward,
                  width = nchar(p$forward_primer[1]))
    mcols(fw)$seq <- p$forward_primer[1]
    mcols(fw)$type <- 'forward'
    rv <- IRanges(start = p$start_reverse,
                  width = nchar(p$reverse_primer[1]))
    mcols(rv)$seq <- p$reverse_primer[1]
    mcols(rv)$type <- 'reverse'
    pair <- c(fw,rv)
    mcols(pair)$seqId <- p$identifier
    mcols(pair)$productSize <- p$product_size
    pair
  })
  names(res) <- primers$identifier
  res
}

.toGenome <- function(x, circStrand, upExon, downExon) {
  vapply(x, function(y) {
    if (circStrand == '+') {
      if (y <= width(upExon)) {
        res <- y + start(upExon) - 1
      } else {
        res <- y - width(upExon) + start(downExon) - 1
      }
      res
    } else {
      if (y <= width(upExon)) {
        res <- end(upExon) - y + 1
      } else {
        res <-  end(downExon) - (y - width(upExon)) + 1
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
#' @param upexon a GRanges object
#' @param downexon a GRanges object
#'
#' @return 
#'
circ2genome <- function(x, upExon, downExon) { 
  circStrand <- unique(strand(upExon))
  # transform to the genome coords and order as c(min, max) using `range`
  res <- mapply(range,
    lapply(
      list(start(x), end(x)),
      .toGenome,
      upExon = upExon,
      downExon = downExon,
      circStrand = circStrand
    )
  )
  res <- GRanges(
    seqnames = Rle(unique(seqnames(upExon)), length(x)),
    IRanges(start = res[1, ], end = res[2, ]),
    strand = Rle(unique(strand(upExon)), length(x))
  )
  mcols(res) <- mcols(x)
  res <- lapply(res, splitPrimer, upExon = upExon, downExon = downExon)
  do.call(c, res)
}

# We assume that primer length is > 2
splitPrimer <- function(x, upExon, downExon) {
  if (overlapsAny(x, upExon) && overlapsAny(x, downExon)) {
    x <- GenomicRanges::narrow(x, start = 2, end = -2)
    sp <- GenomicRanges::setdiff(c(upExon, downExon),x)
    mcols(sp) <- mcols(x)
    sp
  } else {
    x
  }
}

