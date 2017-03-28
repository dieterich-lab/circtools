context("test plotting of isoforms")


test_that("function plot anything", {
  dat <- testData()
  expect_output(
    plotTranscripts(
      exons = dat$exons,
      counts = dat$counts ,
      primers = dat$primers,
      circs = dat$circs,
      minRatio = .1
    )
  )
})
