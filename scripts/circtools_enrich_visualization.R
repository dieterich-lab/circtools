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

########################################################################################################################
### This is the # circRNAs per RBP plot
### We count all isoforms as one circRNA here!
########################################################################################################################

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
    colnames(total) <- c("RBP", "FrequencyA", "FrequencyB", "sum")
    label_pos_1 <- (max(total$FrequencyA, na.rm = TRUE))
    label_pos_2 <- (max(total$FrequencyA, na.rm = TRUE))
} else {
    total$sum <- total$Frequency
    colnames(total) <- c("RBP", "FrequencyA","sum")
    label_pos_1 <- (max(total$FrequencyA, na.rm = TRUE))

}

# limit to top X
total <- head(total,arg_max_rbps)

# start actual plotting
rbp_simple_plot <- ggplot(data=total) +
                        geom_bar(stat="identity", colour="black", size=0.1, aes(x=reorder(RBP, -FrequencyA), y=FrequencyA, fill=RBP)) +
                        geom_label(aes(x = arg_max_rbps - 1 , y = label_pos_1, label = arg_label_sample_1), fill = "white")

                        if (!is.na(arg_data_file_2)) {
                            rbp_simple_plot <- rbp_simple_plot + geom_bar(stat="identity", colour="black", size=0.1, aes(x=reorder(RBP, -FrequencyA), y=-FrequencyB, fill=RBP))
                            rbp_simple_plot <- rbp_simple_plot + geom_label(aes(x = arg_max_rbps - 1, y = - label_pos_2, label = arg_label_sample_2), fill = "white")

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
                                                " RBPs, showing top ",
                                                arg_max_rbps,
                                                " RBPS: ",
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

########################################################################################################################
### This is the # RBPs per circRNA plot
### We count isoforms specifically and do not sum up!
########################################################################################################################

sample_list <- list(rbp_data_file_1, rbp_data_file_2)
sample_names <- list(arg_label_sample_1, arg_label_sample_2)

for (sample in seq(1 : 2)) {

    current_data <- sample_list[[sample]]

    sub_dataframe <- data.frame(table(current_data$Annotation))
    sub_dataframe$Var1 <- levels(droplevels(sub_dataframe$Var1))
    sub_dataframe <- sub_dataframe[order(- sub_dataframe$Freq),]
    to_run <- sub_dataframe$Var1[1 : arg_max_circRNAs]

    for (top_circ in to_run) {


        current_dataframe <- data.frame(unique(current_data[current_data$Annotation == top_circ, 4 : 5]))
        circrna_plot <- list()

        for (circ_isoform in 1 : nrow(current_dataframe))
        {
            tmp_frame <- data.frame(current_data[which(current_data$Annotation == top_circ &
                current_data$start == current_dataframe[circ_isoform, 1] &
                current_data$stop == current_dataframe[circ_isoform, 2] &
                current_data$observed_input_peaks_circ_rna > 0), c(1, 2, 3, 4, 5, 9)])
            colnames(tmp_frame) <- c("RBP", "Annotation", "chr", "start", "stop", "clip_peaks")
            tmp_frame <- tmp_frame[with(tmp_frame, order(- clip_peaks)),]
            tmp_frame <- head(tmp_frame, arg_max_circRNAs)

            miniframe <- data.frame(tmp_frame$RBP, tmp_frame$clip_peaks)
            colnames(miniframe) <- c("RBP", "clip_peaks")
                miniframe <- ddply(unique(miniframe), .(RBP), transform, border = rep(1, clip_peaks))

                theme_set(theme_grey(base_size = 8))
                circrna_plot[[circ_isoform]] <- ggplot(miniframe, aes(RBP)) +
                    geom_bar(aes(x = reorder(paste(RBP, ": ", clip_peaks, sep = ""), - clip_peaks), clip_peaks, fill = RBP), width = 1, size = 0.15, stat = "identity", color = "white") +
                    scale_y_continuous() +
                    coord_polar() +
                    theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks = element_blank()) +
                    labs(title = paste(sample_names[sample], ":\nComposition of RBP landscape for circRNA", tmp_frame[1, 2]),
                    subtitle = paste("Isoform ", circ_isoform, ": Chromsome ", tmp_frame[1, 3], ", ", commapos(as.integer(tmp_frame[1, 4])), "->", commapos(as.integer(tmp_frame[1, 5])), sep = "")) +
                    labs(y = "#eCLIP peaks within annotated circRNA") +
                    labs(x = "") +
                    labs(caption = paste("circRNAs enriched for RBP peaks compared to their host gene ( p <",
                    arg_pval, ")"))
        }

        if (nrow(current_dataframe) == 1){
            do.call(grid.arrange, c(circrna_plot, ncol = 1))
        } else if (nrow(current_dataframe) == 2){
            do.call(grid.arrange, c(circrna_plot, ncol = 2))
        } else {
            ml <- do.call(marrangeGrob, list(grobs=circrna_plot, nrow = 2, ncol = 2, top=NULL))
            print(ml)
        }
    }
}


########################################################################################################################
### Done with plotting, finishing up
########################################################################################################################

message(paste("Printing to",arg_output_file_name))

# disable "NULL" device message
invisible(capture.output(dev.off()))

# get rif of Rplots.pdf file
invisible(capture.output(file.remove(list.files(pattern = "Rplots.pdf"))))

message("Done")
