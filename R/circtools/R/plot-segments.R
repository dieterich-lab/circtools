
inches_per_line <- function(){
  par("csi") 
}
       
maxLabelWidth <- function(x){
  textLength <- max(strwidth(x, units = "inches", cex = 1))
  textLength / inches_per_line() * 1.3
}

margins <- function(left = 0,
                    right = 0,
                    top = 0,
                    bottom = 0,
                    x = c(left, right),
                    y = c(bottom, top)) {
  op <- par()
  par(mar = c(y[1], x[1], y[2], x[2]))
  op
}

xy_per_in <- function() par("cxy") / par("cin")

which_axis <- function(x = FALSE, y = FALSE) {
  op <- par(xaxt = ifelse(x, "s", "n"),
            yaxt = ifelse(y, "s", "n"))
  op
}

no_axis <- function(){
  which_axis() 
}

no_box <- function(){
  op <- par(bty = "n")
  op
}

getPanelHeight <- function(laneNumber){
  par("csi") * laneNumber 
}

getYLim <- function() par()$usr[3:4]

#' Plots isoforms structure, primer position and isoform counts if provided
#'
#' @param exons an interval data.frame
#' @param circs a counts data.frame
#' @param counts an interval data.frame
#' @param primers an interval data.frame
#' @param minRatio a minimal aspect ratio for plotted segments (height/width)
#' @param opts an option list
#' 
#' @return Used for its side effects. Plots intervals for exons,
#' primers and transcript counts if provided.
#' @export
#'
plotTranscripts <- function(exons,
                            counts,
                            primers = NULL,
                            circs = NULL,
                            minRatio = .2,
                            opts = list()) {
  # pre-defined
  numMarginLines <- 3
  widths <- c(2, 1)
  minSegmentAspect <- .2
  # calculate sizes of panels and segments
  segmentSize <- list(size = getPanelHeight(1),
                      minWidth = getPanelHeight(1) * minSegmentAspect)
  primersNum <- ifelse(missing(primers), 0, length(unique(primers$id)))
  upperPanelHeight <- getPanelHeight(primersNum)
  lowerPanelHeight <- getPanelHeight(numMarginLines + length(unique(exons$transcript_id)))
  # in relative units
  heights <-
    c(upperPanelHeight, lowerPanelHeight) / lowerPanelHeight
  layout(
    matrix(c(2, 1, 4, 3), ncol = 2),
    widths = widths,
    heights = heights,
    respect = TRUE
  )
  # plot exons -- bottom left
  labWidth <- maxLabelWidth(as.character(exons$transcript_id))
  op <- margins(left = labWidth, bottom = numMarginLines)
  # plot segments
  with(
    exons,
    plotRanges(
      ids       = transcript_id,
      starts    = start,
      ends      = end,
      seg_width = segmentSize$size,
      min_width = segmentSize$minWidth
    )
  )
  # add circ rectangles if defined
  isoformsYLim <- getYLim()
  segmentWidthInc <- segmentSize$size * xy_per_in()
  if (!missing(circs))
    with(
      circs,
      annotateCircs(
        ids = id,
        starts = start,
        ends = end,
        seg_width = segmentSize$size,
        alpha = .2
      )
    )
  exonsXLim <- with(exons, range(start, end))
  # plot primers -- upper left
  op <- margins(left = labWidth,
                top = 0,
                bottom = .5)
  with(
    primers,
    plotRanges(
      ids = id,
      starts = start,
      ends = end,
      seg_width = segmentSize$size,
      min_width = segmentSize$minWidth,
      xlim = exonsXLim
    )
  )
  # plot counts -- lower right
  par(bty = "o")
  margins(left   = 1,
          bottom = 3,
          right  = 1)
  with(counts,
       plotCounts(id = transcript_id,
                  count = count,
                  ylim = isoformsYLim))
}

plotRanges <- function(ids,
                       starts,
                       ends,
                       seg_width,
                       min_width = 0,
                       xlim = range(starts, ends)) {
  ylim <- c(.5, .5 + length(levels(ids)))
  no_axis()
  no_box()
  plot(
    0,
    type = "n",
    xlim = xlim,
    ylim = ylim,
    ylab = "",
    xlab = ""
  )
  y_pos <- as.numeric(ids)
  mtext(
    as.character(ids),
    side = 2,
    at   = y_pos,
    las  = 1,
    cex  = 1
  )
  seg_width_y <- seg_width * xy_per_in()[2]
  min_width_x <- xy_per_in()[1] * min_width
  o <- (ends - starts) < min_width_x
  ends[o] <- starts[o] + min_width_x
  rect(
    xleft   = starts,
    ybottom = y_pos - seg_width_y / 2,
    xright  = ends,
    ytop    = y_pos + seg_width_y / 2,
    col     = 1,
    border  = NA
  )
  print(starts)
  print(ends)
  print(y_pos)
}


annotateCircs <- function(ids, starts, ends, seg_width, alpha = .4) {
  stopifnot(length(start) == length(ends))
  stopifnot(length(start) == length(ids))
  stopifnot(seg_width > 0)
  stopifnot(alpha > 0 & alpha <= 1)
  seg_width <- seg_width * xy_per_in()[2]
  colors <- rainbow(length(ids), alpha = alpha)
  ylim <- par()$usr[3:4]
  rect(
    xleft = starts,
    xright = ends,
    ybottom = ylim[1] -  seg_width / 2,
    ytop = ylim[2] + seg_width / 2,
    col = colors,
    border = NA
  )
}

plotCounts <- function(id, count, ylim = c(.5, length(id) + .5)) {
  which_axis(x = TRUE)
  plot(
    x = count + .5,
    y = as.numeric(id),
    pch = 16,
    cex = 1.5,
    ylim = ylim,
    xlim = c(0.5, max(count) * 2.5),
    log = 'x',
    xaxs = "i"
  )
  segments(.5, y0 = as.numeric(id), count + .5, lwd = 2)
}

