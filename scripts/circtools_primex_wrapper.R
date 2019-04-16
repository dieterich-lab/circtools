#!/usr/bin/env Rscript

# suppress loading messages
suppressMessages(library(primex))
suppressMessages(library(formattable))
suppressMessages(library(kableExtra))
suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)

# temporary file from python script is the only argument
data_file_name <- args[1]
product_size <- unlist(lapply(strsplit(args[2],","), as.numeric))
junction_mode <- args[3]
num_pairs <- args[4]


# open file for reading
con  <- file(data_file_name, open = "r")

# create empty data table
data_table <- data.frame()

# read file line by line
while (length(current_line <- readLines(con, n = 1, warn = FALSE)) > 0) {

    # split line by tab
    line_column <- (strsplit(current_line, "\t"))

    # we split exon 2 in the middle and create 2 fake exons to continue with the process
    # because we need to have the BSJ somehow in the sequence

    if (line_column[[1]][3] == ""){
        line_column[[1]][3] <- substr(line_column[[1]][2], 1, nchar(line_column[[1]][2])/2)
        tmp1 <- substr(line_column[[1]][2], 1, nchar(line_column[[1]][2])/2)
        tmp2 <- substr(line_column[[1]][2], nchar(line_column[[1]][2])/2+1, nchar(line_column[[1]][2]))
        line_column[[1]][2] <- tmp1
        line_column[[1]][3] <- tmp2
    }

    # set minimal information for primer design
    seqOpts <-  seqSettings(
                    seqId = line_column[[1]][1],
                    seq = c(line_column[[1]][3], line_column[[1]][2])
                        )

    seqOpts <- productSize(seqOpts, c(product_size[1], product_size[2]))

    # empty overlap list
    seqOpts$SEQUENCE_OVERLAP_JUNCTION_LIST = NULL

    # make sure, that the primer actually COVER the BSJ

    # forward primer should cover BSJ
    if (junction_mode == "f"){
        seqOpts$SEQUENCE_PRIMER_PAIR_OK_REGION_LIST <- paste(nchar(line_column[[1]][3])-15,",30,,",sep="")
    # reverse primer should cover BSJ
    } else if (junction_mode == "r"){
        seqOpts$SEQUENCE_PRIMER_PAIR_OK_REGION_LIST <- paste(",,",nchar(line_column[[1]][3])-15,",30",sep="")

    # left primer only in exon2, right primer only in exon1
    } else {
        seqOpts$SEQUENCE_PRIMER_PAIR_OK_REGION_LIST <- paste(   1,
                                                                nchar(line_column[[1]][3]),
                                                                nchar(line_column[[1]][3])+1,
                                                                nchar(line_column[[1]][2])-1,
                                                                sep=","
                                                            )
    }

    # return 10 pairs at max
    seqOpts$PRIMER_NUM_RETURN = num_pairs

    # mute primer3 output, we only want the variable later
    sink("/dev/null")

    # restrict PCR product size to 50-160 bp
    #productSize(seqOpts, c(50, 180))
    primers <- design(seqOpts, returnStats = TRUE)
    # stop output redirect
    sink()

    # temporary data frame holds the primer results for this circRNA isoform
    tmp_df <- primers$primers[-c(1,2,3,c(12:21))]

    # make sure we found any primers at all
    if (!is.null(primers$options$PRIMER_PAIR_NUM_RETURNED) && primers$options$PRIMER_PAIR_NUM_RETURNED > 0){
        colnames(tmp_df) <- c(  "PRIMER_LEFT_SEQUENCE",
                                "PRIMER_RIGHT_SEQUENCE",
                                "PRIMER_LEFT",
                                "PRIMER_RIGHT",
                                "PRIMER_LEFT_TM",
                                "PRIMER_RIGHT_TM",
                                "PRIMER_LEFT_GC_PERCENT",
                                "PRIMER_RIGHT_GC_PERCENT",
                                "PRIMER_PAIR_PRODUCT_SIZE"
                                )

        # rename rows by circRNA-ID + running number
        rownames(tmp_df) <- paste(line_column[[1]][1], rownames(tmp_df), sep="_")

        # merge together new and processed results
        data_table <- rbind(data_table, tmp_df)

    } else {

        tmp_df <- data.frame("NA","NA","NA","NA","NA","NA","NA","NA","NA")
        colnames(tmp_df) <- c(  "PRIMER_LEFT_SEQUENCE",
                        "PRIMER_RIGHT_SEQUENCE",
                        "PRIMER_LEFT",
                        "PRIMER_RIGHT",
                        "PRIMER_LEFT_TM",
                        "PRIMER_RIGHT_TM",
                        "PRIMER_LEFT_GC_PERCENT",
                        "PRIMER_RIGHT_GC_PERCENT",
                        "PRIMER_PAIR_PRODUCT_SIZE"
                        )

        rownames(tmp_df) <- paste(line_column[[1]][1], rownames(tmp_df), sep="_")
        # merge together new and processed results
        data_table <- rbind(data_table, tmp_df)
    }


}

# close file connection
close(con)

# write tab-separated output to screen -> captured by parent python script
write.table(data_table,
            file = "",
            quote = FALSE,
            sep = "\t",
            eol = "\n",
            na = "NA",
            dec = ".",
            row.names = TRUE,
            col.names = FALSE
            )
