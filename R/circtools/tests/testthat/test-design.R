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
  # primer is inside of exon
  expect_equal(downExon,
               circtools:::splitPrimer(downExon, upExon = upExon, downExon = downExon))
  # primer is splitted
  p <- range(c(upExon, downExon))
  start(p) <- 5
  end(p) <- 1005
  expectedPrimer <- c(downExon, upExon)
  start(expectedPrimer) <- c(1, 1005)
  end(expectedPrimer) <- c(5, 1010)
  
  expect_equal(expectedPrimer,
               circtools:::splitPrimer(p, upExon = upExon, downExon = downExon))
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
               ranges(circ2genome(x, upExon, downExon)))
  x <- IRanges(15,19)
  expect_equal(IRanges(5, 9),
               ranges(circ2genome(x, upExon, downExon)))
  
  x <- IRanges(9,19)
  expect_equal(IRanges(c(1, 1009), c(9, 1010)),
               sort(ranges(circ2genome(x, upExon, downExon))))
  
})