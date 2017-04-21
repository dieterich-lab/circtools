
.getCircSJSeq <- function(exSeq) {
  df <- mcols(exSeq)
  strand <- as.character(strand(exSeq))
  stream <- .getStream(strand, df$side)
  indByStream <- expand.grid(split(seq_along(stream), stream))
  res <-  lapply(seq_along(indByStream$downstream),
                 function(i) {
                   up <- indByStream$upstream[i]
                   down <- indByStream$downstream[i]
                   c(
                     up_exon_id = df$exon_id[up],
                     down_exon_id = df$exon_id[down],
                     circSeq = paste0(df$seq[down], df$seq[up]),
                     seq = paste0(df$seq[up], df$seq[down])
                   )
                 })
  as.data.frame(do.call(rbind, res), stringsAsFactors = FALSE)
}

getCircSeqFromList <- function(exSeqList) {
  seqs <- do.call(rbind, lapply(exSeqList, .getCircSJSeq))
  seqs$CIRCID <- row.names(seqs)
  rownames(seqs) <- NULL
  seqs$circSeq_id <- as.character(seq_along(seqs$CIRCID))
  seqs$seq_id <- paste0(seqs$circSeq_id, "-linear")
  seqs
}

seqs <- getCircSeqFromList(exSeqAll)

dbConn <- dbConnect(SQLite(), "circdb")
circDNA <- DNAStringSet(seqs$circSeq)
Seqs2DB(circDNA, "XStringSet", dbConn, identifier = seqs$circSeq_id)
linearDNA <- DNAStringSet(seqs$seq)
Seqs2DB(circDNA, "XStringSet", dbConn, identifier = seqs$seq_id)


desc <- dbGetQuery(dbConn, "select description from Seqs")
Add2DB(data.frame(identifier=desc), dbConn)
 tiles <- TileSeqs(dbConn, add2tbl="Tiles", minCoverage=1)
 primers <- DesignPrimers(tiles, identifier="1",
 	minCoverage=1, minGroupCoverage=1)
