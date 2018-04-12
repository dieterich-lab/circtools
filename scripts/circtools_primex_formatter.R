#!/usr/bin/env Rscript

# suppress loading messages
suppressMessages(library(formattable))
suppressMessages(library(kableExtra))
suppressMessages(library(dplyr))

# switch to red warning color if more blast hits are found
high_count_number = 5

args <- commandArgs(trailingOnly = TRUE)

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
            max-width: 30%; /* Max Width of the popover (depending on the container!) */
        }
    </style>

</head>
<body>
"
#############################################################################################################

# read data file name from args
data_file_name <- args[1]

# read whole file into data table
data_table <- read.csv(data_file_name, header = FALSE, sep = "\t")

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
                            "BLAST_right"
                            )

# generate a column with BLAST hit counts
data_table$BLAST_left_count <- lengths(regmatches(data_table$BLAST_left, gregexpr(";", data_table$BLAST_left)))
data_table$BLAST_right_count <- lengths(regmatches(data_table$BLAST_right, gregexpr(";", data_table$BLAST_right)))

# replace ; with HTML linebreaks for hover popover text
data_table$BLAST_left <- gsub(";", "<br/><br/>", data_table$BLAST_left)
data_table$BLAST_right <- gsub(";", "<br/><br/>", data_table$BLAST_right)

# clean up rownames to hide them lateron
rownames(data_table) <- c()

# main output table generation
output_table <- data_table %>%
    mutate(
    Product_size = color_bar('lightblue')(Product_size),

    Specificity_left = cell_spec(paste(BLAST_left_count, "off-site BLAST hits"),
    popover = spec_popover(content = BLAST_left, title = "Blast Hits\"data-html=\"True\"", position = "right"),
    background = ifelse(BLAST_left_count > high_count_number, "red", "darkgreen"),
    color = ifelse(BLAST_left_count > high_count_number, "white", "white")),

    Sequence_Left <- cell_spec(Left_, escape = F),
    Sequence_Right <- cell_spec(Right_, escape = F),

    Specificity_right = cell_spec(paste(BLAST_right_count, "off-site BLAST hits"),
    popover = spec_popover(content = BLAST_right, title = "Blast Hits\"data-html=\"True\"", position = "left"),
    background = ifelse(BLAST_right_count > high_count_number, "red", "darkgreen"),
    color = ifelse(BLAST_right_count > high_count_number, "white", "white")),

    TM_left = color_tile('white', 'red')(TM_left),
    TM_right = color_tile('white', 'red')(TM_right),
    GC_left = color_tile('white', 'orange')(GC_left),
    GC_right = color_tile('white', 'orange')(GC_right),
    Strand = formatter('span', style = style(font.weight = "bold"))(Strand)
    ) %>%
    select(- Left_) %>%
    select(- Right_) %>%
    select(- BLAST_left) %>%
    select(- BLAST_right) %>%
    select(- BLAST_left_count) %>%
    select(- BLAST_right_count) %>%
    select(Annotation, everything()) %>%
    kable("html", escape = F) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "responsive"), full_width = T) %>%
    # column_spec(5, width = "3cm")
    add_header_above(c("Input circRNA" = 5, "Designed Primers" = 9)) %>%
    # group_rows("Group 1", 4, 7) %>%
    # group_rows("Group 1", 8, 10)
    collapse_rows(columns = 1)

write(paste(html_header, output_table, sep=""), file = "")