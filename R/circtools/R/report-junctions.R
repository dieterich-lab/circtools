

.htmlColumn <- function(vec) {
  lapply(vec, function(cell)
    htmltools::tags$td(htmltools::HTML(as.character(cell))))
}

.addColumnAttr <- function(col, attr, value, filter=seq_along(col)) {
  col[filter] <- lapply(col[filter], function(x) {
    if (!is.null(x$attribs[[attr]])) {
      x$attribs[[attr]] <- paste(x$attribs[[attr]], value)
    } else {
      x$attribs[[attr]] <- value
    }
    x
  })
  col
}

.toRow <- function(colList) {
  names(colList) <- NULL
  f <- function(...) {
   cols <- list(...)
   # remove NA due to spanned cells
   cols <- cols[!vapply(cols, anyNA, logical(1))]
   htmltools::tags$tr(cols)
  }
  do.call(Map, c(f=f, colList)) 
}

.getStream <- function(strand, side) {
  res <- rep("upstream", length(strand))
  res[ strand == '+' & side == 'left'] <- 'downstream'
  res[ strand == '-' & side == 'right'] <- 'downstream'
  res[ strand == '*' ] <- "n.a."
  res
}

getReportTemplate <- function(x) {
  tmpl <- system.file("templates", x, package = "circtools")
  readLines(tmpl)
}


reportCircs <- function(exSeq, file, template) {
  ex <- unlist(exSeq, use.names = FALSE)
  ex <- as.data.frame(ex)
  ex$stream <- .getStream(strand = ex$strand, side = ex$side)
  o <- order(ex$CIRCID, ex$stream=='downstream', ex$start, ex$end)
  ex <- ex[o,]
  if (missing(template))
    template <- system.file("templates", "report-circ.html",
                            package = "circtools")
  cols <- c(
    'CIRCID' = 'CIRCID',
    'exon_id' = 'exon id',
    'seqnames' = 'chr',
    'start' = 'start',
    'end' = 'end',
    'strand' = 'strand',
    'stream' = 'relative position',
    'seq' = 'sequence'
  )
  tdCols <- lapply(ex[names(cols)], .htmlColumn)
  tdCols$seq <- .addColumnAttr(tdCols$seq, 'class', 'seq')
  tdCols$seq <- .addColumnAttr(
    tdCols$seq, 'class', 'upstream', filter = ex$stream == 'upstream') 
  tdCols$seq <- .addColumnAttr(
    tdCols$seq, 'class', 'downstream', filter = ex$stream == 'downstream') 
  tdCols$CIRCID <- spanCells(tdCols$CIRCID, ex$CIRCID)
  header <- lapply(cols, htmltools::tags$th)
  names(header) <- NULL
  header <- htmltools::tags$tr(header)
  tbl <- htmltools::tags$table(header, .toRow(tdCols))
  reportBody <- htmltools::tags$body(tbl)
  html <- htmltools::htmlTemplate(file = template, body = reportBody)
  dep <- htmltools::htmlDependency(name = "circtools",
                        version = utils::packageVersion("circtools"),
                        src = c(file = system.file(package = "circtools")), 
                        stylesheet = "templates/css/seq-table.css")
  res <- renderTags(html)$html
  res <- attachDependencies(res, dep)
  htmltools::save_html(res, file)
}

#' @param x a list of tags
#' @param f a factor to span the cells
spanCells <- function(x, f) {
  r <- rle(f)
  indexesToKeep <- c(1, 1 + cumsum(r$lengths))
  indexesToKeep <- indexesToKeep[-length(indexesToKeep)]
  for (i in seq_along(indexesToKeep)) {
    x[[indexesToKeep[i]]]$attribs$rowspan <- r$lengths[i]
  }
  x[-indexesToKeep] <- NA
  x
}
