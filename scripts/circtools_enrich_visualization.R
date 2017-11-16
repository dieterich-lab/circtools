#!/usr/bin/env Rscript

# Copyright (C) 2017 Tobias Jakobi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# packages (silent load)
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(data.table))
suppressMessages(library(reshape2))
suppressMessages(library(plyr))
suppressMessages(library(gridExtra))


# from https://learnr.wordpress.com/2009/09/24/ggplot2-back-to-back-bar-charts/
# return the - now - negative y counts as positive
commapos <- function(x, ...) {
    format(abs(x), big.mark = ",", trim = TRUE,
    scientific = FALSE, ...)
}

# may not be optimal, but we do not want warning for now
options(warn = - 1)

# quiet mode
options(echo = FALSE)

# read arguments
args <- commandArgs(trailingOnly = TRUE)

arg_data_file_1 <- args[1]
arg_pval <- as.numeric(args[2])
arg_max_circRNAs <- as.integer(args[3])
arg_max_rbps <- as.integer(args[4])
arg_output_file_name <- args[5]
arg_label_sample_1 <- args[6]

# arguments for file 2 are ath the send so we can leave them empty
arg_data_file_2 <- args[7] # path is string
arg_label_sample_2 <- args[8] # path is string

# default output file name
if (is.na(arg_output_file_name)) {
    arg_output_file_name <- "./output.pdf"
}

# default label for file 1
if (is.na(arg_label_sample_1)) {
    arg_label_sample_1 <- "Sample 1"
}

# default label for file 2: should not be used however
if (is.na(arg_label_sample_2)) {
    arg_label_sample_2 <- "Sample 2"
}

message("Reading input file 1")
rbp_data_file_1 <- read.delim(arg_data_file_1, check.names = FALSE, header = F)

# colnames of processed circtools out files (adding 1 column for RBP)
colnames(rbp_data_file_1) <- c( "RBP",
                                "Annotation",
                                "chr",
                                "start",
                                "stop",
                                "strand",
                                "p_val_circular",
                                "raw_count_circ_rna",
                                "observed_input_peaks_circ_rna",
                                "length_circ_rna",
                                "length_normalized_count_circ_rna",
                                "number_of_features_intersecting_circ",
                                "circ_rna_confidence_interval_0.05",
                                "p_val_linear", "raw_count_host_gene",
                                "observed_input_peaks_host_gene",
                                "length_host_gene_without_circ_rna",
                                "length_normalized_count_host_gene",
                                "number_of_features_intersecting_linear",
                                "host_gene_confidence_interval_0.05",
                                "distance_normalized_counts"
)

rbp_data_file_1 <- rbp_data_file_1[rbp_data_file_1$p_val_circular < arg_pval & rbp_data_file_1$p_val_linear > arg_pval,]


if (!is.na(arg_data_file_2)) {
    message("Reading input file 2")
    rbp_data_file_2 <- read.delim(arg_data_file_2, check.names = FALSE, header = F)
    colnames(rbp_data_file_2) <- colnames(rbp_data_file_1)
    rbp_data_file_2 <- rbp_data_file_2[rbp_data_file_2$p_val_circular < arg_pval & rbp_data_file_2$p_val_linear > arg_pval,]
}


message("plotting data")

# we use PDF and standard A4 page size
pdf(arg_output_file_name, height = 8.2, width = 11.69 , title = paste("circtools RBP enrichment analysis - ", arg_label_sample_1))


tmp <- data.frame(rbp_data_file_1$RBP)
colnames(tmp) <- c("RBP")
tmp$annotation <- rbp_data_file_1$Annotation
tmp <- unique(tmp)
rbps <- table(tmp$RBP)

rbps <- sort(rbps, decreasing = TRUE)
rbp_df <- data.frame(rbps)
colnames(rbp_df) <- c("RBP", "Frequency")
rbp_df <- rbp_df[rbp_df$Frequency > arg_max_circRNAs,]

if (!is.na(arg_data_file_2)) {
    tmp <- data.frame(rbp_data_file_2$RBP)
    colnames(tmp) <- c("RBP")
    tmp$annotation <- rbp_data_file_2$Annotation
    tmp <- unique(tmp)
    rbps2 <- table(tmp$RBP)

    rbps2 <- sort(rbps2, decreasing = TRUE)
    rbp_df2 <- data.frame(rbps2)
    colnames(rbp_df2) <- c("RBP", "Frequency")
    rbp_df2 <- rbp_df2[rbp_df2$Frequency > arg_max_circRNAs,]
}


total <-rbp_df

if (!is.na(arg_data_file_2)) {
    total <- merge(rbp_df, rbp_df2, by = "RBP", all = TRUE)
    total[is.na(total)] <- 0
    colnames(total) <- c("RBP", "FrequencyA", "FrequencyB")
    total$sum <- total$FrequencyA + total$FrequencyB
    colnames(total) <- c("RBP", "FrequencyA", "FrequencyB", "total")
} else {
    total$sum <- total$Frequency
    colnames(total) <- c("RBP", "FrequencyA","total")
}

# limit to top X
total <- head(total,arg_max_rbps)

rbp_simple_plot <- ggplot(data=total) +
                        geom_bar(stat="identity", colour="black", size=0.1, aes(x=reorder(RBP, -total), y=FrequencyA, fill=RBP))

                        if (!is.na(arg_data_file_2)) {
                            rbp_simple_plot <- rbp_simple_plot + geom_bar(stat="identity", colour="black", size=0.1, aes(x=reorder(RBP, -total), y=-FrequencyB, fill=RBP))
                        }
                        rbp_simple_plot <- rbp_simple_plot +
                        labs(   title = paste(arg_label_sample_1, ": Number of circular RNAs per RBP", sep=""),
                                subtitle = paste("Counting circular RNAs (including different isoforms) with significantly enriched RBPs ( p <",
                                arg_pval, ")")) +
                        labs(y = "Number of circular RNAs") +
                        labs(x = "RNA binding protein") +
                        scale_y_continuous(labels = commapos) +
                        labs(caption = paste(   "based on data from ",
                                                length(unique(rbp_data_file_1$RBP)),
                                                " RBPs, showing RBPs with > ",
                                                arg_max_circRNAs,
                                                " circRNAs: ",
                                                date(),
                                                "",
                                                sep="")) +

                        theme(  plot.title = element_text(lineheight=0.8,
                                face=quote(bold)),
                                legend.justification = c(1, 1),
                                legend.position = "none",
                                axis.text.x = element_text(angle = 90, hjust = 1, size=14)

                        )


    print(rbp_simple_plot)

message(paste("Printing to",arg_output_file_name))

# disable "NULL" device message
invisible(capture.output(dev.off()))

# get rif of Rplots.pdf file
invisible(capture.output(file.remove(list.files(pattern = "Rplots.pdf"))))

message("Done")
