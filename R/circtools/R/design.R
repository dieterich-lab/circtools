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
exIntersect <- getTxForCircs(db, exSeq, 'both')
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
  maxSearchSize = 20,
  minProductSize = 60,
  minEfficiency = .9
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

decipher2genome <- function(primers, seqs) {
  primersIR <- decipher2iranges(primers)
  p <- primersIR[[1]]
  seqId <- names(primersIR)[1]
  i <- which(seqs$seqId == seqId)
  circ <- exSeq[[seqs$CIRCID[i]]]
  up <- circ[which(mcols(circ)$exon_id == seqs$upExonId[i])]
  down <- circ[which(mcols(circ)$exon_id == seqs$downExonId[i])]
    
}

decipher2iranges <- function(primers) {
  res <- lapply(seq_along(primers$identifier), function(i) {
    p <- primers[i, ]
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
    pair
  })
  names(res) <- primers$identifier
  res
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
  if (unique(strand(upExon))[1] == '+') {
    toGenome <- function(x) {
      if (x <= width(upExon)) {
        res <- x + start(upExon) - 1
      } else {
        res <- x - width(upExon) + start(downExon) - 1
      }
      res
    }
  } else {
    toGenome <- function(x) {
      if (x <= width(upExon)) {
        res <- end(upExon) - x + 1
      } else {
        res <- x - width(upExon) + end(downExon) - 1
      }
      res
    }
  }
  res <- mapply(range, vapply(start(x), toGenome, numeric(1)),
              vapply(end(x), toGenome, numeric(1)))
  res <- GRanges(seqnames = Rle(unique(seqnames(up)), length(x)),
    IRanges(start = res[1,], end = res[2,]),
    strand = Rle(unique(strand(up)), length(x)))
  mcols(res) <- mcols(x)
  res <- lapply(res, splitPrimer, upExon = upExon, downExon = downExon)
  do.call(c, res)
}

splitPrimer <- function(x, upExon, downExon) {
  if (overlapsAny(x, upExon) && overlapsAny(x, downExon)) {
    sp <- setdiff(range(c(upExon, downExon)), range(x))
    mcols(sp) <- mcols(x)
    sp
  } else {
    x
  }
}

p <- decipher2iranges(primers)
x <- p[[1]]

p <- DNAStringSet(c(primers$forward_primer[1],  primers$reverse_primer[1]))
#rp <- reverseComplement(p)
AmplifyDNA((p),
           DNAStringSet(seqs$circSeq),
           maxProductSize = 1e4,annealingTemp = 64,P = 4e-7)
AmplifyDNA((p),
           DNAStringSet(ts),
           maxProductSize = 1e4,annealingTemp = 64,P = 4e-7)
CalculateEfficiencyPCR(p[1],
                       reverseComplement(DNAStringSet(seqs$circSeq[1])),
                       #reverseComplement(p[1]),
                       temp = 64,
                       P = 4e-7,
                       ions = .225)


