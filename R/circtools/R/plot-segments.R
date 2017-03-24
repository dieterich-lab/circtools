plotSegments <- function(data, primers) {
    geneSize <- 5
    aspect   <- 4
    width    <- 12
    columns  <- c("start_pos", "end_pos")
    norms    <- normaliseData(data=data, primers=primers)
    primers[, columns] <- norms$primers
    data[, columns]    <- norms$data
    data$y_pos         <- as.factor(data$transcript_id)
# get lines for transcripts
    trans <- data                       %>%
        group_by(transcript_id)         %>%
        summarise(
            start_pos = min(start_pos),
            end_pos   = max(end_pos),
            y_pos     = unique(y_pos),
            count     = max(count))
    env=environment()
# define plotting themes
    rmLines <- theme(
            panel.grid.minor.x      = element_blank(),
            panel.border            = element_blank())
    rmXAxis <- theme(axis.text.x    = element_blank(),
            axis.ticks.x            = element_blank(),
            panel.grid.major.y      = element_blank(),
            panel.grid.minor.y      = element_blank())
    sc <- scale_x_continuous(limits = c(min(data$start_pos), max(data$end_pos)))
# plot segments and lines
    qGenes <- ggplot(environment=env) +
    geom_segment(data = trans,
            aes(x    = start_pos,
                xend = end_pos,
                y    = y_pos,
                yend = y_pos)) + 
    geom_segment(data=data,
            aes(x    = start_pos,
                xend = end_pos,
                y    = y_pos,
                yend = y_pos),
            size=geneSize) + 
    xlab(NULL) + ylab(NULL) + 
    sc + theme_bw() + rmLines + rmXAxis + 
    theme(panel.grid.major.y   = element_blank(),
            panel.grid.minor.y = element_blank())
    qExpr <- ggplot(environment=env) +
        geom_point(data = trans,
                aes(x=count, y=y_pos)) + xlab(NULL)  + 
        scale_x_log10() +
        theme_bw()      +
        rmLines         +
        theme(axis.text.y    = element_blank(),
                axis.title.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.line.y  = element_line(),
                axis.line.x  = element_line())
    qPrim <- ggplot(env=env)  +
        geom_segment(
                data    = primers,
                mapping = aes(
                    x    = start_pos,
                    xend = end_pos,
                    y    = name,
                    yend = name 
                    ),
                color   = "green",
                size = geneSize) +
        geom_line(
                data    = primers,
                mapping = aes(
                    x     = start_pos,
                    y     = name,
                    group = name),
                size = 1) +sc + xlab(NULL) + ylab(NULL) + 
        theme_bw() +
        rmLines + rmXAxis
    qPrim   <- ggplotGrob(qPrim)
    qGenes  <- ggplotGrob(qGenes)
    qExpr   <- ggplotGrob(qExpr)
    heights <- grid::unit.pmax(qExpr$heights, qGenes$heights)
    maxwidths <- grid::unit.pmax(qPrim$widths[2:5], qGenes$widths[2:5])
    qGenes$widths[2:5] <- as.list(maxwidths)
    qPrim$widths[2:5] <- as.list(maxwidths)
    qGenes$heights <- qExpr$heights <- heights
    # plot and estimate height
    hPrimers <- length(unique(primers$name)) * 11 * 2# * 2.56/72
    hGenes   <- length(unique(data$transcript_id)) *11 * 2# * 2.56/72
    width    <- hGenes + hPrimers
    wRight   <- .3 * width * aspect
    wLeft    <- .7 * width * aspect
    arrangeGrob(qPrim, qGenes, qExpr,
            layout_matrix=rbind(c(1,NA),c(2,3)),
            widths=unit(c(wLeft,wRight),"pt"),
            heights=unit(c(hPrimers, hGenes), "pt"))
}

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
                     cbind(transcript_id = paste0("ENS000000", x),
                           exons[sample(exonsNum, sample(exonsNum, 1) +
                                          1), ])
                   })
  exons <- do.call(rbind, exons)
  res$exons <- exons
  # create expression levels
  res$counts <- data.frame(transcript_id = unique(exons$transcript_id),
                           count = round(2 ^ runif(
                             transNum, min = -10, max = 20
                           )))
  
  # primers
  t1 <- exons$transcript_id[1]
  e1 <- exons[exons$transcript_id == t1, ]
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
