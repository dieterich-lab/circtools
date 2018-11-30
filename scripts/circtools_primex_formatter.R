#!/usr/bin/env Rscript

# suppress loading messages
suppressMessages(library(formattable))
suppressMessages(library(kableExtra))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(colortools))

# switch to red warning color if more blast hits are found
high_count_number = 0

args <- commandArgs(trailingOnly = TRUE)

experiment_name <- args[2]

# set output to HTML
options(knitr.table.format = 'html')

# generic HTML header with bootstrap JS libraries
#############################################################################################################
html_header="
<html>
<head>

<meta charset=\"utf-8\">
<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">

  <link rel=\"stylesheet\" href=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css\">
  <script src=\"https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js\"></script>
  <script src=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js\"></script>

<script>
$(document).ready(function(){
    $('[data-toggle=\"popover\"]').popover();
});
</script>

<script>
$(document).ready(function(){
    $('[data-toggle=\"tooltip\"]').tooltip();
});
</script>

   <style type=\"text/css\">
        /* The max width is dependant on the container (more info below) */
        .popover{
            max-width: 50%; /* Max Width of the popover (depending on the container!) */
        }
    </style>

</head>
<body>"

html_header <- paste(html_header,"<h1>circtools primer design results for experiment ", experiment_name , "</h1>", sep="")

#############################################################################################################

# generate a divergent color scale with 11 shades
color_palette <- rev(brewer.pal(n = 5, name = 'RdBu'))

#default TM value
default_tm_value <- 60
default_gc_value <- 50
default_product_value <- 100

construct_color_column <- function(column, default_value, palette)
{
    top_val <- (max(column, na.rm=T) - default_value)
    bottom_val <- (default_value - min(column, na.rm=T))

    if (top_val > bottom_val){
        from <- default_value - top_val
        to <- default_value + top_val
    } else if (top_val == bottom_val) {
        from <- default_value - 1
        to <- default_value + 1
    } else {
        from <- default_value - bottom_val
        to <- default_value + bottom_val
    }

    return(as.character(cut(column,seq( from, to, length.out= length(palette)+1 ), labels=palette, include.lowest = T)))
}

generate_complementary_column <- function (input){
    return (complementary(input)[2])
}

# read data file name from args
data_file_name <- args[1]

# read whole file into data table
data_table <- read.csv(data_file_name, header = FALSE, sep = "\t")
data_table$circid <- paste(data_table$V1,data_table$V2,data_table$V3,data_table$V4,data_table$V5,data_table$V6,sep="_")

data_table$circid <- paste(sep="","<img src=",data_table$circid,".svg>")

# remove unused columns
data_table <- data_table[-c(6,9,10)]

# correctly name the columns
colnames(data_table) <- c(  "Annotation",
                            "Chr",
                            "Start",
                            "Stop",
                            "Strand",
                            "Left_",
                            "Right_",
                            "TM_left",
                            "TM_right",
                            "GC_left",
                            "GC_right",
                            "Product_size",
                            "BLAST_left",
                            "BLAST_right",
                            "ID"
                            )

data_table$right_tm_color  = construct_color_column(data_table$TM_right,default_tm_value,color_palette)
data_table$left_tm_color   = construct_color_column(data_table$TM_left,default_tm_value,color_palette)

data_table$left_gc_color   = construct_color_column(data_table$GC_left,default_gc_value,color_palette)
data_table$right_gc_color  = construct_color_column(data_table$GC_right,default_gc_value,color_palette)

data_table$product_color  = construct_color_column(data_table$Product_size,default_product_value,color_palette)


colnames_final <- c(        "Annotation",
                            "Chr",
                            "Start",
                            "Stop",
                            "Strand",
                            "TM forward",
                            "TM reverse",
                            "GC% forward",
                            "GC% reverse",
                            "Product size",
                            "Forward",
                            "BLAST",
                            "Reverse",
                            "BLAST"
                    )

# run_primer_design a column with BLAST hit counts
data_table$BLAST_left_count <- lengths(regmatches(data_table$BLAST_left, gregexpr(";", data_table$BLAST_left))) + 1
data_table$BLAST_right_count <- lengths(regmatches(data_table$BLAST_right, gregexpr(";", data_table$BLAST_right))) + 1

data_table$BLAST_left_count[data_table$BLAST_left_count == 1] = 0
data_table$BLAST_right_count[data_table$BLAST_right_count == 1] = 0

# replace ; with HTML linebreaks for hover popover text
data_table$BLAST_left <- gsub(";", "<br/><br/>", data_table$BLAST_left)
data_table$BLAST_right <- gsub(";", "<br/><br/>", data_table$BLAST_right)

# remove 0 entries from location columns for provided circRNA FASTA files
data_table$Chr <- gsub("\\b0\\b", "", data_table$Chr )
data_table$Start <- gsub("\\b0\\b", "", data_table$Start )

data_table$Stop  <- ifelse( data_table$Start == "" , "", data_table$Stop )

data_table$Strand<- gsub("\\b0\\b", "", data_table$Strand )

# clean up rownames to hide them lateron
rownames(data_table) <- c()

# main output table generation
output_table <- data_table %>%
    mutate(
    Product_size = color_bar(product_color)(Product_size),

    Forward = cell_spec(escape = F, Left_, popover = spec_popover( title = "Graphical represensation of designed primers and annotated circRNA structure\"data-html=\"True\"", position = "left", content =ID ), background = ifelse(BLAST_left_count > high_count_number, "red", "darkgreen"),
    color = ifelse(BLAST_left_count > high_count_number, "white", "white")),

    L = cell_spec(paste(BLAST_left_count),
    popover = spec_popover(content = BLAST_left, title = "Blast Hits\"data-html=\"True\"", position = "right"),
    background = ifelse(BLAST_left_count > high_count_number, "red", "darkgreen"),
    color = ifelse(BLAST_left_count > high_count_number, "white", "white"), bold = "true"),

    Reverse = cell_spec(escape = F, Right_, popover = spec_popover( title = "Graphical represensation of designed primers and annotated circRNA structure\"data-html=\"True\"", position = "left", content =ID ), background = ifelse(BLAST_right_count > high_count_number, "red", "darkgreen"),
    color = ifelse(BLAST_right_count > high_count_number, "white", "white")),

    R = cell_spec(paste(BLAST_right_count),
    popover = spec_popover(content = BLAST_right, title = "Blast Hits\"data-html=\"True\"", position = "left"),
    background = ifelse(BLAST_right_count > high_count_number, "red", "darkgreen"),
    color = ifelse(BLAST_right_count > high_count_number, "white", "white"), bold = "true"),

    TM_left = color_bar(left_tm_color)(TM_left),
    TM_right = color_bar(left_tm_color)(TM_right),


    GC_left = color_bar(left_gc_color)(GC_left),
    GC_right = color_bar(right_gc_color)(GC_right),


    Strand = formatter('span', style = style(font.weight = "bold"))(Strand)
    ) %>%
    select(- Left_) %>%
    select(- Right_) %>%
    select(- BLAST_left) %>%
    select(- BLAST_right) %>%
    select(- BLAST_left_count) %>%
    select(- BLAST_right_count) %>%
    select(- ID) %>%
    select(- right_tm_color) %>%
    select(- left_tm_color) %>%
    select(- right_gc_color) %>%
    select(- left_gc_color) %>%
    select(- product_color) %>%
    select(Annotation, everything()) %>%
    kable("html", escape = F, col.names=colnames_final) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "responsive"), full_width = T) %>%
    # column_spec(5, width = "3cm")
    add_header_above(c("Input circRNAs" = 5, "Designed Primers" = 9)) # %>%
    # group_rows("Group 1", 4, 7) %>%
    # group_rows("Group 1", 8, 10)
    # collapse_rows(columns = 1)

write(paste(html_header, output_table, sep=""), file = "")