contex("Primer design")

test_that("correct coordinate transform", {
  downExon <- GRanges(seqnames = 'a',
                    strand = Rle('+', 1),
                    ranges = IRanges(start = 1, end = 10))
  upExon <- GRanges(seqnames = 'a',
                    strand = Rle('+', 1),
                    ranges = IRanges(start = 1001, end = 1010))
  # primer is inside of exon
  expect_equal(downExon,
               splitPrimer(downExon, upExon = upExon, downExon = downExon))
  # primer is splitted
  p <- range(c(upExon, downExon))
  start(p) <- start(p) + 5
  end(p) <- end(p) - 5
  expectedPrimer <- disjoin(c(range(c(upExon, downExon)),p))[c(1,3)]
  expect_equal(expectedPrimer,
    splitPrimer(p, upExon = upExon, downExon = downExon))
  # inside to genome coords
  expect_equal(1004,
  .toGenome(4, circStrand = '+',upExon = upExon, downExon = downExon))
  expect_equal(4,
  .toGenome(14, circStrand = '+',upExon = upExon, downExon = downExon))
  # swap and minus strand
  expect_equal(7,
  .toGenome(4, circStrand = '-',upExon = downExon, downExon = upExon))
  expect_equal(1007,
  .toGenome(14, circStrand = '-',upExon = downExon, downExon = upExon))
  
  
})