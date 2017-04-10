
unifyDiff <- function(x,y, ratio){
    points <- cbind(x,y)
    o <- order(points)
    deltas <- diff(points[o])
    nonZeroDeltas <- deltas[deltas>0]
    logs <- log(nonZeroDeltas/min(nonZeroDeltas))
    logs <- logs / max(logs) * log(ratio)
    # back
    nonZeroDeltas <- exp(logs)
    deltas[deltas > 0] <- nonZeroDeltas
    points[o] <- cumsum(c(1,deltas))
    list(points[seq_along(x)],
         points[seq_along(y) + length(x)])
}

normaliseData <- function(..., ratio=5){
    dat <- list(...)
    columns <- c("start_pos", "end_pos")
    positions <- do.call(rbind,
        lapply(names(dat), function(x){
            cbind(id=x,dat[[x]][, columns])
        })
    )
    result <- as.data.frame(
        do.call(cbind,
            unifyDiff(positions$start_pos, positions$end_pos, ratio=ratio))
    )
    names(result) <- columns
    split(result, positions$id)
}

testData <- function() {
  set.seed(239)
  res <- list()
  # create transcript-exon table
  exonsNum <- 10
  exons <- data.frame(start = 1000 * 1:exonsNum,
                      end = 1000 * 1:exonsNum + 600)
  transNum <- 5
  exons <-  lapply(1:transNum,
                   function(x) {
                     cbind(tx_id = paste0("ENST000000", x),
                           exons[sample(exonsNum, sample(exonsNum, 1) +
                                          1), ])
                   })
  exons <- do.call(rbind, exons)
  res$exons <- exons
  # create expression levels
  res$counts <- data.frame(id = unique(exons$tx_id),
                           count = round(2 ^ stats::runif(
                             transNum, min = -10, max = 20
                           )))
  
  # primers
  t1 <- exons$tx_id[1]
  e1 <- exons[exons$tx_id == t1, ]
  e1 <- e1[order(e1$start), ]
  res$primers <- data.frame(id = paste0("circ", c(1, 1, 2)),
                            rbind(
                              c(e1$end[1] - 10, e1$end[1]),
                              c(e1$start[2], e1$start[2] + 10),
                              c(e1$start[2] + 40, e1$start[2] + 60)
                            ))
  names(res$primers) <- c("id", "start", "end")
  # circ coord
  res$circs <-  data.frame(id = "circ1",
                           start = e1$start[4],
                           end = e1$end[6])
  #circ primers
  res$circPrimers <-
    with(res$circs,
         {
           data.frame(id = paste0("circ", c(1, 1, 2)),
                      rbind(
                        c(start[1], start[1] + 10),
                        c(end[1] - 10, end[1]),
                        c(start[1] + 60, start[1] + 80)
                      ))
         })
  names(res$circPrimers) <- names(res$primers)
  res
}
