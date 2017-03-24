source("code/vis/plot-segments.R")

inches_per_line <- function(){
  par("csi") 
}
       
maxLabelWidth <- function(x){
  textLength <- max(strwidth(x, units = "inches", cex=1))
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

dat <- testData()

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
 op <- par(bty="n")
 op
}


plotTr <- function(exons,
                   counts,
                   primers = NULL,
                   circs = NULL,
                   seg_w = 1,
                   min_ratio = .2) {
  seg_width <-  seg_w * par("csi")
  nprimers <- ifelse(missing(primers), 0, length(unique(primers$id)))
  h1 <- nprimers * seg_width * 1.1# + 1 * par("csi")
  h2 <-
    length(unique(exons$transcript_id)) * seg_width * 1.0 + 3 * par("csi")
  h1 <- h1/h2
  layout(
    matrix(c(2, 1, 4, 3), ncol = 2),
    widths = c(2, 1),
    heights = c(h1,1),
    respect = TRUE
  )
  lab.width <- maxLabelWidth(as.character(exons$transcript_id))
  op <- margins(left = lab.width, bottom = 3)
  with(exons, plotRanges(transcript_id, start, end, seg_width))
  ylim_isoforms <- par()$usr[3:4]
  xy_seg <- seg_width * xy_per_in()
  if(!missing(circs))
    with(circs, annotateCircs(id, start, end, seg_width))
  xlims_exons <- with(exons,range(start, end))
  # primers
  op <- margins(left=lab.width, top = 0, bottom = .5)
  with(primers,
       plotRanges(id, start, end, seg_width, xlim = xlims_exons,
                 min_width = min_ratio)
  )
  # counts
  par(bty = "o")
  margins(left = 1,
          bottom = 3,
          right = 1)
  with(counts, plotCounts(transcript_id,count,ylim = ylim_isoforms))
}

plotRanges <- function(ids,
           starts,
           ends,
           seg_width,
           xlim = range(starts, ends),
           min_width=0) {
  ylim <- c(.5, .5 + length(levels(ids)))
  par()$fin
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
    at = y_pos,
    las = 1,
    cex = 1
  )
  seg_width_y <- seg_width * xy_per_in()[2]
  min_width_x <- seg_width * xy_per_in()[1] * min_width
  o <- (ends-starts) < min_width_x
  ends[o] <- starts[o] + min_width_x
  rect(
    xleft = starts,
    ybottom = y_pos - seg_width_y / 2,
    xright  = ends,
    ytop = y_pos + seg_width_y / 2,
    col = 1
    ,
    border = NA
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


plotTr(
  dat$exons,
  dat$counts ,
  primers = dat$primers,
  circs = dat$circs,
  seg_w = 1
)
