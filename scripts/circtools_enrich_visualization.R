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
suppressMessages(library(ggrepel)) # beautifies text labels
suppressMessages(library(data.table))
suppressMessages(library(reshape2))

# may not be optimal, but we do not want warning for now
options(warn = - 1)

# quiet mode
options(echo = FALSE)

# read arguments
args <- commandArgs(trailingOnly = TRUE)

arg_data_file <- args[1] # path is string
arg_pval <- as.numeric(args[2])
arg_min_circRNAs <- as.integer(args[3]) # path is string
arg_min_rbps <- as.integer(args[4]) # path is string
arg_min_top_colors <- as.integer(args[5]) # path is string
arg_output <- args[6] # path is string
arg_experiment <- args[7] # path is string


if (is.na(arg_output)){
    arg_output <- "./output.pdf"
}

if (is.na(arg_experiment)){
    arg_experiment <- "default"
}

message("Loading data file")
rbp_data <- read.delim(arg_data_file, check.names=FALSE, header = F)

colnames(rbp_data) <- c("RBP", "Annotation", "chr", "start", "stop", "strand", "p_val_circular", "raw_count_circ_rna", "observed_input_peaks_circ_rna", "length_circ_rna", "length_normalized_count_circ_rna", "number_of_features_intersecting_circ", "circ_rna_confidence_interval_0.05", "p_val_linear", "raw_count_host_gene", "observed_input_peaks_host_gene", "length_host_gene_without_circ_rna", "length_normalized_count_host_gene", "number_of_features_intersecting_linear", "host_gene_confidence_interval_0.05", "distance_normalized_counts")

# head(rbp_data)
message("plotting data")

rbp_data <- rbp_data[rbp_data$p_val_circular < arg_pval & rbp_data$p_val_linear > arg_pval, ]


# we use PDF and standard A4 page size
pdf(arg_output, height= 8.2, width=11.69 , title=paste("circtools RBP enrichment analysis - ",arg_experiment))

    # rbp plot

    tmp <- data.frame(rbp_data$RBP)
    colnames(tmp) <-c("RBP")
    tmp$annotation <- rbp_data$Annotation
    tmp <- unique(tmp)

    rbps <- table(tmp$RBP)
    rbps <- sort(rbps, decreasing = TRUE)
    rbp_df <- data.frame(rbps)
    colnames(rbp_df) <-c("RBP","Frequency")
    rbp_df <- rbp_df[rbp_df$Frequency>arg_min_circRNAs,]

    rbp_simple_plot <- ggplot(data=rbp_df, aes(x=reorder(RBP, -Frequency), y=Frequency, fill=RBP)) +
                        geom_bar(stat="identity", colour="black") +

                        labs(   title = "Number of circular RNAs per RBP",
                                subtitle = paste("Counting circular RNAs (including different isoforms) with significantly enriched RBPs ( p <",
                                arg_pval, ")")) +
                        labs(y = "Number of circular RNAs") +
                        labs(x = "RNA binding protein") +
                        labs(caption = paste(   "based on data from ",
                                                length(unique(rbp_data$RBP)),
                                                " RBPs, showing RBPs with > ",
                                                arg_min_circRNAs,
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

    # annotation plot

    annotation <- table(rbp_data$Annotation)
    annotation <- sort(annotation, decreasing = TRUE)
    annotation_df <- data.frame(annotation)
    colnames(annotation_df) <-c("Annotation","Frequency")
    annotation_df <- annotation_df[annotation_df$Frequency>arg_min_rbps,]


    circ_simple_plot <- ggplot(data=annotation_df, aes(x=reorder(Annotation, -Frequency), y=Frequency, fill=Annotation)) +
                        geom_bar(stat="identity", colour="black") +

                        labs(   title = "Number of RBP sites per circRNA",
                                subtitle = paste("Counting circular RNAs (including different isoforms) with significantly enriched RBPs ( p <",
                                arg_pval, ")")) +
                        labs(y = "Number of enriched binding sites") +
                        labs(x = "Circular RNA") +
                        labs(caption = paste(   "based on data from ",
                                                length(unique(rbp_data$Annotation)),
                                                " circRNAs, showing circRNAs with > ",
                                                arg_min_rbps,
                                                " RPBs: ",
                                                date(),
                                                "",
                                                sep="")) +
                        theme(  plot.title = element_text(lineheight=0.8,
                                face=quote(bold)),
                                legend.justification = c(1, 1),
                                legend.position = "none",
                                axis.text.x = element_text(angle = 90, hjust = 1, size=14)

                        )
    print(circ_simple_plot)


    # mixed plot

    heading <- unname(transpose(data.frame((table(rbp_data[rbp_data$Annotation,]$RBP))))[1,])
    heading <- unlist(heading)

    tmp_list <- list()
    counter <- 1

    for (circRNA in (rownames(annotation))){
        row <- as.numeric(unname(transpose(data.frame((table(rbp_data[rbp_data$Annotation==circRNA,]$RBP)))))[2,])
        tmp_list[[counter]] <- unlist(row)
        counter <- counter+1
    }
    message(paste(counter,"circRNAs processed"))
    circRNA_df <- as.data.frame(do.call("rbind", tmp_list))
    circRNA_df <- data.frame(rownames(annotation), circRNA_df)
    colnames(circRNA_df) <- c("CircRNA",heading)
    circRNA_df_m <- melt(head(circRNA_df,arg_min_top_colors),id.vars = "CircRNA")
    circRNA_df_m$CircRNA <- factor(circRNA_df_m$CircRNA, levels = circRNA_df_m$CircRNA)

    circ_simple_plot <- ggplot(data=circRNA_df_m, aes(x=CircRNA, y = value, fill = variable)) +
                        geom_bar(stat="identity") +
                        labs(   title = "Assignment of circRNAs to RPBs",
                                subtitle = paste("plotting colour-coded RBPs to allow quick visulation of different groups per circRNA ( p <",
                                            arg_pval, ")")) +
                        labs(y = "Number of enriched binding sites") +
                        labs(x = "RNA binding protein") +
                        labs(caption = paste(   "based on data from ",
                                                length(unique(rbp_data$Annotation)),
                                                " circRNAs and ",
                                                length(unique(rbp_data$RBP)),
                                                " RBPs, showing top ",
                                                arg_min_top_colors,
                                                " RPBs: ",
                                                date(),
                                                "",
                                                sep="")) +
                        theme(  plot.title = element_text(lineheight=0.8,
                                face=quote(bold)),
                                legend.justification = c(1, 1),
                                legend.position = "none",
                                axis.text.x = element_text(angle = 90, hjust = 1, size=14)

                        )
    print(circ_simple_plot)



# disable "NULL" device message
invisible(capture.output(dev.off()))

# get rif of Rplots.pdf file
invisible(capture.output(file.remove(list.files(pattern = "Rplots.pdf"))))

message("Done")
