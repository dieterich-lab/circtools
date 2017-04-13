
.template_env <- list(
        head="<!DOCTYPE html>
            <html>
            <head>
            <title>Sequences</title>
            <style>
            table {
                font-family: arial, sans-serif;
                font-size: 14px;
                border-collapse: collapse;
            }

            td, th {
                border: 1px solid #dddddd;
                text-align: left;
                padding: 8px;
            }

            </style>
            </head>
    ",body=NULL,
    end="</html>")

htmlTable <- function(...) {
 cols <- list(...)
 cols <- lapply(cols, function(x) {
   lapply(x, function(cell) htmltools::tags$td(HTML(cell)))
   })
 names(cols) <- NULL
 toRow <- function(...) Map(f = htmltools::tags$tr, ...)
 do.call(toRow, cols)
}

reportCircs <- function(exSeq, file) {
  ex <- unlist(exSeq, use.names = FALSE)
  ex <- as.data.frame(ex)
  cols <- c(
    'CIRCID' = 'CIRCID',
    'seqnames' = 'chr',
    'start' = 'start',
    'end' = 'end',
    'strand' = 'strand',
    'seq' = 'sequence'
  )
  tableBody <- do.call(htmlTable, ex[names(cols)])
  header <- lapply(cols, tags$th)
  names(header) <- NULL
  doc <- .template_env
  doc$body <- as.character(
    renderTags(tags$table(header, tableBody))$html)
  result <- paste(doc, collapse = "\n ")
  cat(result, file=file)
}

