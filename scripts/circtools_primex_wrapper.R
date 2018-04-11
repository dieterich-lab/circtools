#!/usr/bin/env Rscript

suppressMessages(library(primex))
suppressMessages(library(formattable))
suppressMessages(library(kableExtra))
suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)

data_file_name <- args[1]

con  <- file(data_file_name, open = "r")

data_table <- data.frame()

while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    myVector <- (strsplit(oneLine, "\t"))
        circ_data <- (strsplit(myVector[[1]][1], "_"))
    #print(myVector[[1]][1])

    seqOpts <-  seqSettings(
                        seqId = myVector[[1]][1],
                        seq = c(myVector[[1]][3], myVector[[1]][2])
                        )

    seqOpts$SEQUENCE_OVERLAP_JUNCTION_LIST = NULL
    seqOpts$SEQUENCE_PRIMER_PAIR_OK_REGION_LIST <- paste(1,nchar(myVector[[1]][3])-10,1,nchar(myVector[[1]][2]),sep=",")
    seqOpts$PRIMER_NUM_RETURN = 10

    sink("/dev/null")
    productSize(seqOpts, c(50, 160))
    primers <- design(seqOpts, returnStats = TRUE)
    sink()

    tmp_df <- primers$primers[-c(1,2,3,c(12:21))]

    if (primers$options$PRIMER_PAIR_NUM_RETURNED > 0){
        rownames(tmp_df) <- paste(myVector[[1]][1], rownames(tmp_df), sep="_")
        data_table <- rbind(data_table, tmp_df)
    }
}


close(con)

write.table(data_table, file = "", quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = FALSE)
