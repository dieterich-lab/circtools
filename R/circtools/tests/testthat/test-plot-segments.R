context("test plotting of isoforms")


test_that("function plots anything", {
  dat <- testData()
  expect_silent(
    plotTranscripts(
      exons = dat$exons,
      counts = dat$counts ,
      primers = dat$primers,
      circs = dat$circs,
      minAspectRatio = .1
    )
  )
})

test_that("no error if only part is supplied", {
  dat <- testData()
  expect_silent(
    plotTranscripts(
      exons = dat$exons,
      #counts = dat$counts ,
      #primers = dat$primers,
      #circs = dat$circs,
      #minRatio = .1
    )
  )
  expect_silent(
    plotTranscripts(
      exons = dat$exons,
      counts = dat$counts ,
      #primers = dat$primers,
      #circs = dat$circs,
    )
  )
  expect_silent(
    plotTranscripts(
      exons = dat$exons,
      #counts = dat$counts ,
      #primers = dat$primers,
      circs = dat$circs,
    )
  )
  expect_silent(
    plotTranscripts(
      exons = dat$exons,
      #counts = dat$counts ,
      primers = dat$primers,
      circs = dat$circs,
    )
  )
  
})
