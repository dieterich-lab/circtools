context("Create and plot CircData")

library("EnsDb.Hsapiens.v86")
db <- EnsDb.Hsapiens.v86

createCirc <- function(geneName, db){
  txs <- select(db, keys=geneName, keytype="GENENAME",columns="TXNAME")
  ex1 <- exonsBy(db, filter=TxidFilter(txs$TXNAME[1]))[[1]]
  ex2 <- exonsBy(db, filter=TxidFilter(txs$TXNAME[2]))[[1]]
  circCoords <- c(range(ex1[1:2]), range(ex2[1:2]))
  mcols(circCoords)$CIRCID <- paste0(seqnames(circCoords), ":",
                                     start(circCoords), "-", end(circCoords))
  circCoords
}

test_that("No error when create a CircData object", {
  geneName <- c("BCL6", "BCL2")
  suppressWarnings(
    circCoords <- do.call(c, lapply(geneName, createCirc, db = db)))
  circId <- mcols(circCoords)$CIRCID[1]
  expect_silent(CircData(db, circCoords))
})

test_that("No error while plot CircData object", {
  library("EnsDb.Hsapiens.v86")
  db <- EnsDb.Hsapiens.v86
  geneName <- c("BCL6" )
  suppressWarnings(
    circCoords <- do.call(c, lapply(geneName, createCirc, db = db)))
  circId <- mcols(circCoords)$CIRCID[1]
  circData <- CircData(db, circCoords)
  expect_silent(plotCirc(circId, circData = circData))
})
