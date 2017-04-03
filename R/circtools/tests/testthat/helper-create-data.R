
createGene <- function(coordOffset = 0,
                       idOffset = 0) {
  exons <- data.frame(exon_start = c(0, 100, 100, 200, 300) + coordOffset)
  exons$exon_end <- exons$exon_start + c(50, 50, 70, 50, 50)
  txMap <- rbind(
    data.frame(
      tx_name = "tx_1",
      tx_id = 1L,
      exon_rank = 1:3,
      ex_id = c(1, 2, 4)
    ),
    data.frame(
      tx_name = "tx_2",
      tx_id = 2L,
      exon_rank = 1:3,
      ex_id = c(1, 3, 4)
    ),
    data.frame(
      tx_name = "tx_3",
      tx_id = 3L,
      exon_rank = 1:3,
      ex_id = c(3, 4, 5)
    )
  )
  txMap$tx_id <- as.integer(txMap$tx_id + idOffset)
  splicings <- cbind(txMap[, -4], exons[txMap$ex_id, ])
  txMap$ex_id <- as.integer(txMap$ex_id + idOffset)
  transcripts <- do.call(rbind,
                         lapply(split(splicings, splicings$tx_name),
                                function(x)
                                  data.frame(
                                    tx_start = min(x$exon_start),
                                    tx_end = max(x$exon_end),
                                    tx_strand = "+",
                                    tx_name = unique(x$tx_name),
                                    tx_id = unique(x$tx_id),
                                    tx_chrom = "chr1"
                                  )))
  splicings$tx_name <- NULL
  list(splicings = splicings, transcripts = transcripts)
}


createTestTxDB <- function() {
  gene_id <- "gene_1"
  tx <- createGene()
  genes <- data.frame(
    gene_id = gene_id,
    tx_id = tx$transcripts$tx_id
  )
  makeTxDb(tx$transcripts, tx$splicings, genes)
}

createTestCircs <- function(txdb) {
  ex <- exonsBy(txdb, by = "tx", use.names = TRUE)[[1]]
  res <- range(ex[sort(sample(length(ex), size = 2, replace = TRUE))])
  res <- as.data.frame(res)
  names(res)[1] <- "chr"
  res
}