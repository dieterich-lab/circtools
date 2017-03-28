context("test plotting of isoforms")


test_that("function plots anything", {
  dat <- testData()
  expect_silent(
    plotTranscripts(
      exons = dat$exons,
      counts = dat$counts ,
      primers = dat$primers,
      circs = dat$circs,
      minRatio = .1
    )
  )
})
