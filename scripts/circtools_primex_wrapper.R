#!/usr/bin/env Rscript

suppressMessages(library(primex))
library(formattable)

args <- commandArgs(trailingOnly = TRUE)

html_header="
<head>
<meta charset=\"utf-8\">
<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">
<link rel=\"stylesheet\" href=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css\">

</head>
<body>
"

color_bar_NA <- function(color = "lightgray", fun = "proportion", ...) {
    fun <- match.fun(fun)
    replace_na <- function(fun, x, ...) {
        x[which(is.na(x))] = 0
        return(fun(as.numeric(x), ...))
    }
    formatter("span",
        style = function(x) style(display = "inline-block",
            direction = "rtl",
            "border-radius" = "4px",
            "padding-right" = "2px",
            "background-color" = csscolor(color),
            width = percent(replace_na(fun, x, ...))
        )
    )
}



data_file_name <- args[1]

con  <- file(data_file_name, open = "r")

data_table <- data.frame()

while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    myVector <- (strsplit(oneLine, "\t"))
        circ_data <- (strsplit(myVector[[1]][1], "_"))

    seqOpts <-  seqSettings(
                        seqId = myVector[[1]][1],
                        seq = c(myVector[[1]][3], myVector[[1]][2])
                        )

    seqOpts$SEQUENCE_OVERLAP_JUNCTION_LIST = NULL

    seqOpts$PRIMER_NUM_RETURN = 10

    sink("/dev/null")
    productSize(seqOpts, c(50, 160))
    primers <- design(seqOpts, returnStats = FALSE)
    sink()

    tmp_df <- primers$primers[-c(1,2,3,c(12:21))]


    tmp_df$Host <- circ_data[[1]][1]
    tmp_df$Chr <- circ_data[[1]][2]
    tmp_df$Start <- circ_data[[1]][3]
    tmp_df$Stop <- circ_data[[1]][4]
    tmp_df$Strand <- circ_data[[1]][5]

    tmp_df <- tmp_df[c(c(10:14),c(1:9))]


    colnames(tmp_df) <- c("Annotation", "Chr", "Start", "Stop", "Strand", "Left", "Right", "Pos_left", "Pos_right", "TM_left", "TM_right", "GC_left", "GC_right","Product")

    rownames(tmp_df) <- paste(myVector[[1]][1], rownames(tmp_df), sep="_")

    data_table <- rbind(data_table, tmp_df)
}

widget_formattable = formattable(data_table, list( Product  = color_bar('lightblue'),
                                                    Left = formatter('span', style=style(font.family = "monospace", font.size ="16")),
                                                        Right = formatter('span', style=style(font.family = "monospace", font.size ="16")),
                                                        TM_left = color_tile('white', 'red'),
                                                        TM_right = color_tile('white', 'red'),
                                                        GC_left = color_tile('white', 'orange'),
                                                        GC_right = color_tile('white', 'red'),
                                                        Strand = formatter('span', style=style(font.weight = "bold"))
                                                        ))



html_table = format_table(widget_formattable)

write(paste(html_header, html_table, sep=""), "./bla.html")


close(con)






