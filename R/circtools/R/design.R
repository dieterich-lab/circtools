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

seqs <- getCircSeqFromList(exSeq)
txIntersect <- getTxForCircs(db, exSeq, 'both')
txSeqs <- getTxSeqs(db, bsg, exSeq)

circId <- names(txIntersect)[1]

dbConn <- dbConnect(SQLite(), ":memory:")
o <- seqs$CIRCID == circId
circDNA <- DNAStringSet(seqs$circSeq[o])
Seqs2DB(circDNA, "XStringSet", dbConn, identifier = seqs$seqId[o])
# add corresponding tx's
ts <- txSeqs[txIntersect[[circId]]]
Seqs2DB(ts, "XStringSet", dbConn, identifier = names(ts))
tiles <- TileSeqs(dbConn, add2tbl="Tiles", minCoverage=1)

seqId <- seqs$seqId[1]
primers <- DesignPrimers(
  tiles = tiles,
  identifier = seqId,
  minCoverage = 1,
  minGroupCoverage = 1,
  numPrimerSets = 1,
  maxSearchSize = 20
)
cols <- c(
'start_forward',
'start_reverse',
'forward_primer',
'reverse_primer',
'forward_efficiency',
'reverse_efficiency'
)
res <- lapply(primers, function(x) `[[`(x,1))[cols]
# include option to intersect SJ?