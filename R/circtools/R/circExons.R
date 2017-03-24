
library(dplyr)
library(biomaRt)
library(Biostrings)


get_exons <- function(t_chr, t_start, t_end, mart){
    s <- getSequence(chr=t_chr, start=t_start, end=t_end, mart=human,
        type=c(
            "ensembl_gene_id",
            "ensembl_exon_id",
            "ensembl_transcript_id", 
            "chromosome_name",
            "exon_chrom_start", 
            "exon_chrom_end",
            "strand", 
            "external_gene_name"), 
        seqType="gene_exon") 
    left_exon<- s %>% filter(exon_chrom_start==t_start) 
    right_exon<- s %>% filter(exon_chrom_end==t_end) 
    list(left_exon=left_exon,right_exon=right_exon)
}

makeFields <- function(e){
    if(dim(e$right_exon)[1]==0||dim(e$left_exon)[1]==0) return(NULL)
    res <- list()
    res$ensembl_id <- e$right_exon$ensembl_gene_id[1]
    res$gene_name <- e$right_exon$external_gene_name[1]
    res$chromosome <- e$right_exon$chromosome_name
    res$strand <- unique(e$right_exon$strand)
    left_exon <- e$left_exon%>% 
        dplyr::select(gene_exon,exon_chrom_start, exon_chrom_end) %>%
                        top_n(n=1,exon_chrom_end) %>% unique()
    right_exon <- e$right_exon %>% 
        dplyr::select(gene_exon,exon_chrom_start, exon_chrom_end) %>%
                        top_n(n=-1,exon_chrom_start) %>% unique()
    if (res$strand==1){
        tmp <- left_exon
        left_exon <- right_exon
        right_exon <- tmp
    }
    res$left_coords<- paste(
        left_exon$exon_chrom_start, "-",left_exon$exon_chrom_end)
    res$right_coords <- paste(
        right_exon$exon_chrom_start, "-",right_exon$exon_chrom_end)
    res$left_seq <- left_exon$gene_exon
    res$right_seq<- right_exon$gene_exon
    res
}

getRow <- function(row){
    gene_name <- td(row$gene_name)
    ensembl_id <- td(row$ensembl_id)
    chromosome <- td(row$chromosome)
    coords <-  c(
        td(span(row$left_coords, style="color:red")),
        td(span(row$right_coords, style="color:blue")))
    strand <- td(c("+","-")[row$strand])
    s1span <- span(row$left_seq, style="color:red")
    s2span <- span(row$right_seq, style="color:blue")
    seq <- td(div(s1span, s2span, 
        style="word-wrap: break-word;width:1000px; font-family:monospace"),
    style="background-color:  #fff4c8 ")
    tr(gene_name, ensembl_id,strand, chromosome, coords, seq)
}
template_env <-list(
    doc=list(
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
    end="</html>"),
    reportTitle=
"  cDNA sequences of adjacent exons at circular splice junctions",
    description= 
"<span style=\"color:red\"> exons BEFORE the junction are in RED </span>;
<span style=\"color:blue\"> exons AFTER the junction are in BLUE </span>" 
)
writeReport <- function(file, mart, coords,
                reportTitle=template_env$reportTitle, 
                description=template_env$description){
    doc <- template_env$doc
    header <- tr(
            th("ensembl"),
            th("gene"),
            th("strand"),
            th("chr"),
            th("left exon coordinates"),
            th("right exon coordinates"),
            th("sequence"))
    description <- paste("<h1>", reportTitle, "</h1>",
                   "<p>", description, "</p>")
    exons  <- apply(X=coords, MARGIN=1,
        FUN=function(x){ print(x);get_exons(x[1], x[2], x[3], human)})
    results <- lapply(exons, makeFields)
    o <- order(unlist((lapply(results, function(x) x$gene_name))))
    doc$body <-body(table(header,description,lapply(results[o], getRow)))
    cat(file=file, paste(doc, collapse="\n"))
}

### from Hadley Wickham book

named <- function(x) {
  if (is.null(names(x))) return(NULL)
  x[names(x) != ""]
}
unnamed <- function(x) {
  if (is.null(names(x))) return(x)
  x[names(x) == ""]
}

tag <- function(tag) {
  force(tag)
  function(...) {
    args <- list(...)
    attribs <- html_attributes(named(args))
    children <- unlist((unnamed(args)))

    paste0(
      "<", tag, attribs, ">",
      paste(children, collapse = ""),
      "</", tag, ">"
    )
  }
}

html_attributes <- function(list) {
  if (length(list) == 0) return("")

  attr <- Map(html_attribute, names(list), list)
  paste0(" ", unlist(attr), collapse = "")
}

html_attribute <- function(name, value = NULL) {
  if (length(value) == 0) return(name) # for attributes with no value
  if (length(value) != 1) stop("value must be NULL or of length 1")

  if (is.logical(value)) {
    # Convert T and F to true and false
    value <- tolower(value)
  } else {
    value <- escape_attr(value)
  }
  paste0(name, " = '", value, "'")
}

escape_attr <- function(x) {
  x <- escape.character(x)
  x <- gsub("\'", '&#39;', x)
  x <- gsub("\"", '&quot;', x)
  x <- gsub("\r", '&#13;', x)
  x <- gsub("\n", '&#10;', x)
  x
}

escape.character <- function(x) {
  x <- gsub("&", "&amp;", x)
  x <- gsub("<", "&lt;", x)
  x <- gsub(">", "&gt;", x)
  x
}


body <-tag("body")
p <- tag("p")
pre <- tag("pre")
span <- tag("span")
table <- tag("table")
tr <- tag("tr")
td <- tag("td")
th <- tag("th")
style <- tag("style")
div <- tag("div")