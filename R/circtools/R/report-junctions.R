

.htmlColumn <- function(vec) {
  lapply(vec, function(cell)
    htmltools::tags$td(htmltools::HTML(as.character(cell))))
}

.addColumnAttr <- function(col, attr, value) {
  lapply(col, function(x) {
    x$attribs[[attr]] <- value
    x})
}

.toRow <- function(colList) {
  names(colList) <- NULL
  do.call(Map, c(f=htmltools::tags$tr, colList)) 
}

.isUpstream <- function(strand, side) {
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
  ex$isUpstream <- .isUpstream(strand = ex$strand, side = ex$side)
  if (missing(template))
    template <- system.file("templates", "report-circ.html",
                            package = "circtools")
  cols <- c(
    'CIRCID' = 'CIRCID',
    'seqnames' = 'chr',
    'start' = 'start',
    'end' = 'end',
    'strand' = 'strand',
    'isUpstream' = 'relative position',
    'seq' = 'sequence'
  )
  tdCols <- lapply(ex[names(cols)], .htmlColumn)
  tdCols$seq <- .addColumnAttr(tdCols$seq, 'class', 'seq')
  header <- lapply(cols, htmltools::tags$th)
  names(header) <- NULL
  header <- htmltools::tags$tr(header)
  tbl <- htmltools::tags$table(header, .toRow(tdCols))
  reportBody <- htmltools::tags$body(tbl)
  html <- htmlTemplate(file = template, body = reportBody)
  dep <- htmlDependency(name = "circtools",
                        version = utils::packageVersion("circtools"),
                        src = c(href = system.file(package = "circtools")), 
                        stylesheet = "/templates/css/seq-table.css")
 html <- renderDocument(html, dep) 
 writeLines(html, con = file)
}


spanCells <- function(x) {
  
}
