
# get height of a text line in inches
inches_per_line <- function(){
  graphics::par("csi") 
}

# calculate text width in lines
maxLabelWidth <- function(x){
  textLength <- max(graphics::strwidth(x, units = "inches", cex = 1))
  textLength / graphics::par()$csi
}

# set plot outer margins
margins <- function(left = 0,
                    right = 0,
                    top = 0,
                    bottom = 0,
                    x = c(left, right),
                    y = c(bottom, top)) {
  op <- graphics::par()
  graphics::par(mar = c(y[1], x[1], y[2], x[2]))
  op
}

xy_per_in <- function() graphics::par("cxy") / graphics::par("cin")

# set which axis to plot
which_axis <- function(x = FALSE, y = FALSE) {
  op <- graphics::par(xaxt = ifelse(x, "s", "n"),
            yaxt = ifelse(y, "s", "n"))
  op
}

# drop axes
no_axis <- function(){
  which_axis() 
}

# no outer box around a plot
no_box <- function(){
  op <- graphics::par(bty = "n")
  op
}

getPanelHeight <- function(laneNumber){
  graphics::par("csi") * laneNumber 
}

# get vertical limits in user coordinates 
getYLim <- function() graphics::par()$usr[3:4]

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
                            counts = NULL,
                            primers = NULL,
                            circs = NULL,
                            minAspectRatio = .2,
                            opts = list()) {
  # TODO clean naming convention --> id, ids => seqnames ranges object?
  toDF <- function(x, idColumnName = "id") {
    if (!is.null(x) && is(x, "GRangesList")) {
      grList2df(x, idColumnName)
    } else if (!is.null(x)) {
      as.data.frame(x)
    } else
      NULL
  }
  exons <- toDF(exons, "tx_id")
  circs <- toDF(circs, "id")
  primers <- toDF(primers)
  # pre-defined
  numMarginLines <- 3
  widths <- c(2, 1)
  if (is.null(counts)) 
    widths[2] <- .1
  minAspectRatio <- .2
  # calculate sizes of panels and segments
  segmentSize <- .75
  primersNum <- ifelse(missing(primers), 0, length(unique(primers$id)))
  upperPanelHeight <- getPanelHeight(primersNum)
  lowerPanelHeight <- getPanelHeight(numMarginLines + length(unique(exons$id)))
  # in relative units
  heights <- c(upperPanelHeight, lowerPanelHeight) / lowerPanelHeight + .1
  graphics::layout(
    matrix(c(2, 1, 4, 3), ncol = 2),
    widths = widths,
    heights = heights
  )
  # plot exons -- bottom left
  # leave 0.5 + 0.5 = 1 margin lines around the labels
  labWidth <- maxLabelWidth(as.character(exons$tx_id)) + 1 
  op <- margins(left = labWidth, bottom = numMarginLines)
  # plot segments
  if (is.null(exons$tx_id))
    stop("No transcript id field named 'tx_id' in the `exons` argument")
  plotRanges(
    ids         = exons$tx_id,
    starts      = exons$start,
    ends        = exons$end,
    segmentSize = segmentSize,
    minWidth    = minAspectRatio,
    opts = opts
  )
  # add circ rectangles if defined
  isoformsYLim <- getYLim()
  if (!is.null(circs))
    annotateCircs(
      ids    = circs$CIRCID,
      starts = circs$start,
      ends   = circs$end,
      alpha = .1
    )
  exonsXLim <- range(exons$start, exons$end)
  # plot primers -- upper left
  if (!is.null(primers)) {
    graphics::abline(h = graphics::par()$usr[4] + .1,
           col = "deepskyblue1",
           lwd = inches_per_line()  * 96)
    op <- margins(left = labWidth,
                  top = 0,
                  bottom = .1)
    plotRanges(
      ids    = primers$id,
      starts = primers$start,
      ends   = primers$end,
      segmentSize = segmentSize,
      minWidth = minAspectRatio,
      xlim = exonsXLim,
      opts = list(col = "firebrick3")
    )
  } else {
    margins()
    graphics::plot.new()
  }
  # plot counts -- lower right
  if (!is.null(counts)) {
    graphics::par(bty = "u")
    margins(left   = 1,
            bottom = numMarginLines,
            right  = 0.2)
    plotCounts(id    = counts$id,
               count = counts$count)
  }
}

#' Plots segments for a list of intervals
#'
#' @param ids a character vector
#' @param starts a numeric vector
#' @param ends a numerica vector
#' @param segmentSize a segment size (height) in inches, so it can be the same
#' for several subplots
#' @param minWidth a minimal segment width in inches
#' @param xlim a range of interval coordinates on the plot. Used for alignment
#' of features at multiple plots
#'
#' @return Used for its side effect. Plots segments corresponding to given 
#' intervals. 
#' @export
#'
plotRanges <- function(ids,
                       starts,
                       ends,
                       segmentSize,
                       minWidth = 0,
                       xlim = range(starts, ends),
                       opts = list()) {
  if (is.null(opts$col)) opts$col <- "dodgerblue4"
  ids <- as.factor(ids)
  ylim <- c(0.5, .5 + length(levels(ids)))
  no_axis()
  no_box()
  graphics::plot(
    0,
    type = "n",
    xlim = xlim,
    ylim = ylim,
    ylab = "",
    xlab = ""
  )
  y_pos <- as.numeric(ids)
  graphics::mtext(
    as.character(ids),
    side = 2,
    line = .5,
    at   = y_pos,
    las  = 1,
    cex  = graphics::par()$cex
  )
  #seg_width_y <- segmentSize * xy_per_in()[2]
  seg_width_y <- segmentSize 
  #min_width_x <- xy_per_in()[1] * minWidth
  x_to_y <- xy_per_in()[1] / xy_per_in()[2]
  min_width_x <- segmentSize * x_to_y * minWidth
  o <- (ends - starts) < min_width_x
  ends[o] <- starts[o] + min_width_x
  graphics::rect(
    xleft   = starts,
    ybottom = y_pos - seg_width_y / 2,
    xright  = ends,
    ytop    = y_pos + seg_width_y / 2,
    col     = opts$col,
    border  = NA
  )
}


annotateCircs <- function(ids, starts, ends, alpha = .2) {
    stopifnot(length(starts) == length(ends))
    stopifnot(alpha > 0 & alpha <= 1)
    #colors <- rainbow(length(starts), s = .6, alpha = alpha)
    #colorsLine <- rainbow(length(starts), s = 1, alpha = 1)
    colors <- grDevices::adjustcolor("darkseagreen1", alpha=.1)
    colorsLine <- "darkolivegreen4"
    ylim <- graphics::par()$usr[3:4] + c(0, -.5) 
    step <- .5 / length(starts)
    graphics::rect(
      xleft = starts,
      xright = ends,
      ybottom = ylim[1] + step * (seq_along(starts) - 1),
      ytop = ylim[2] + step * (seq_along(starts) - 1), 
      col = colors,
      border = colorsLine,
      lwd = 2
    )
  }

plotCounts <- function(id, count, ylim = c(.5, length(id) + .5)) {
  which_axis(x = TRUE)
  graphics::plot(
    x = count + .5,
    y = as.numeric(id),
    pch = 16,
    cex = 1.5,
    ylim = ylim,
    xlim = c(0.5, max(count) * 2.5),
    log = 'x',
    xaxs = "i"
  )
  graphics::segments(.5, y0 = as.numeric(id), count + .5, lwd = 2)
}


# create a data.frame from a named list of ranges
# used to flatten exon-by-transcript lists to pass to a plotting function
grList2df <- function(grl, idColumn = "id") {
  nrows <- IRanges::elementNROWS(grl)
  id <- rep(names(nrows), times = nrows)
  rangesDF <- as.data.frame(
    GenomicsRanges::ranges(unlist(grl, use.names = FALSE)))
  rangesDF[[idColumn]] <- id
  rangesDF
}
