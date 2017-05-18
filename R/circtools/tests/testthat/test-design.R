context("Primer design")

test_that("correct coordinate transform", {
  downExon <- GRanges(
    seqnames = 'a',
    strand = Rle('+', 1),
    ranges = IRanges(start = 1, end = 10)
  )
  upExon <- GRanges(
    seqnames = 'a',
    strand = Rle('+', 1),
    ranges = IRanges(start = 1001, end = 1010)
  )
  # inside to genome coords
  expect_equal(1004,
               circtools:::.toGenome(
                 4,
                 circStrand = '+',
                 upExon = upExon,
                 downExon = downExon
               ))
  expect_equal(4,
               circtools:::.toGenome(
                 14,
                 circStrand = '+',
                 upExon = upExon,
                 downExon = downExon
               ))
  # swap and minus strand
  expect_equal(7,
               circtools:::.toGenome(
                 4,
                 circStrand = '-',
                 upExon = downExon,
                 downExon = upExon
               ))
  expect_equal(
    1007,
    circtools:::.toGenome(
      14,
      circStrand = '-',
      upExon = downExon,
      downExon = upExon
    )
  )
  
  x <- IRanges(5,9)
  expect_equal(IRanges(1005, 1009),
               ranges(circtools:::circ2genome(x, upExon, downExon)))
  x <- IRanges(15,19)
  expect_equal(IRanges(5, 9),
               ranges(circtools:::circ2genome(x, upExon, downExon)))
  
  x <- IRanges(9,19)
  expect_equal(IRanges(c(1, 1009), c(9, 1010)),
               sort(ranges(circtools:::circ2genome(x, upExon, downExon))))
  
})

test_that("return correct exon seqs", {
  gr <- GRanges(
    seqnames = c('a', 'a'),
    ranges = IRanges(start = c(1, 10), end = c(2, 11)),
    strand = rep('+',2),
    sjId = rep('sjId',2),
    exon_id = c('e1', 'e2'),
    gene_id = rep('gene',2),
    side = c('left', 'right'),
    seq = c(DNAStringSet('AA'), DNAStringSet('GG'))
  )
  res <- .getCircSJSeq(gr) 
  expect_equal(
    res,
    data.frame(
      sjId = 'sjId',
      upExonId = 'e2',
      downExonId = 'e1',
      circSeq = 'GGAA',
      stringsAsFactors = FALSE
    )
  )
  strand(gr) <- '-'
  res <- .getCircSJSeq(gr) 
  expect_equal(
    res,
    data.frame(
      sjId = 'sjId',
      upExonId = 'e1',
      downExonId = 'e2',
      circSeq = 'AAGG',
      stringsAsFactors = FALSE
    )
  )
})
