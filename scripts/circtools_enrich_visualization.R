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
library(data.table)
library(reshape2)

# may not be optimal, but we do not want warning for now
options(warn = - 1)

# quiet mode
options(echo = FALSE)

# read arguments
args <- commandArgs(trailingOnly = TRUE)

arg_data_file <- args[1] # path is string

message("Loading data file")
rbp_data <- read.delim(arg_data_file, check.names=FALSE, header = F)
colnames(rbp_data) <- c("RBP","Annotation","Chr","Start","Stop","Strand","PVal","Simulated","Observed")


message("plotting data")

# we use PDF and standard A4 page size
pdf("test.pdf", height= 8.2, width=11.69 , title="circtools library analysis")


    # rbp plot

    rbps <- table(rbp_data$RBP)
    rbps <- sort(rbps, decreasing = TRUE)
    rbp_df <- data.frame(rbps)
    colnames(rbp_df) <-c("RBP","Frequency")
    rbp_df <- rbp_df[rbp_df$Frequency>40,]


    rbp_simple_plot <- ggplot(data=rbp_df, aes(x=reorder(RBP, Frequency), y=Frequency)) +
                        geom_bar(stat="identity", colour="black") +

                        # labs(   title = "Detected circular RNAs per million unique mapped reads",
                        #         subtitle = "plotting circRNA predictions from DCC and uniquely mapped reads from STAR") +
                        # labs(y = "Number of detected circular RNAs") +
                        # labs(x = "Sequencing library") +
                        # labs(fill = "Group") +
                        # labs(caption = paste(   "based on data from ",
                        #                         length(star_runs),
                        #                         " mapped libraries\ncreated: ",
                        #                         date(),
                        #                         "",
                        #                         sep="")) +
                        theme(  plot.title = element_text(lineheight=0.8,
                                face=quote(bold)),
                                legend.justification = c(1, 1),
                                legend.position = c(0.25, 0.98),
                                axis.text.x = element_text(angle = 90, hjust = 1, size=14)

                        )
                        # geom_label_repel(   data=circle_ratio,
                        #                     aes(x=name, y=num,
                        #                     label=num),
                        #                     box.padding = unit(0.0, "lines"),
                        #                     point.padding = unit(0.0, "lines")
                        #                 )
    print(rbp_simple_plot)

    # annotation plot

    annotation <- table(rbp_data$Annotation)
    annotation <- sort(annotation, decreasing = TRUE)
    annotation_df <- data.frame(annotation)
    colnames(annotation_df) <-c("Annotation","Frequency")
    annotation_df <- annotation_df[annotation_df$Frequency>10,]


    circ_simple_plot <- ggplot(data=annotation_df, aes(x=reorder(Annotation, Frequency), y=Frequency)) +
                        geom_bar(stat="identity", colour="black") +

                        # labs(   title = "Detected circular RNAs per million unique mapped reads",
                        #         subtitle = "plotting circRNA predictions from DCC and uniquely mapped reads from STAR") +
                        # labs(y = "Number of detected circular RNAs") +
                        # labs(x = "Sequencing library") +
                        # labs(fill = "Group") +
                        # labs(caption = paste(   "based on data from ",
                        #                         length(star_runs),
                        #                         " mapped libraries\ncreated: ",
                        #                         date(),
                        #                         "",
                        #                         sep="")) +
                        theme(  plot.title = element_text(lineheight=0.8,
                                face=quote(bold)),
                                legend.justification = c(1, 1),
                                legend.position = c(0.25, 0.98),
                                axis.text.x = element_text(angle = 90, hjust = 1, size=14)

                        )
                        # geom_label_repel(   data=circle_ratio,
                        #                     aes(x=name, y=num,
                        #                     label=num),
                        #                     box.padding = unit(0.0, "lines"),
                        #                     point.padding = unit(0.0, "lines")
                        #                 )
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
    # colnames(circRNA_df) <- heading
    # rownames(circRNA_df) <- (rownames(annotation))
    circRNA_df <- data.frame(rownames(annotation), circRNA_df)
    colnames(circRNA_df) <- c("CircRNA",heading)
    circRNA_df_m <- melt(head(circRNA_df,10),id.vars = "CircRNA")
    circRNA_df_m$CircRNA <- factor(circRNA_df_m$CircRNA, levels = circRNA_df_m$CircRNA)

    print(circRNA_df_m)

    # reorder(name, num)
    circ_simple_plot <- ggplot(data=circRNA_df_m, aes(x=CircRNA, y = value, fill = variable)) +
                        geom_bar(stat="identity") +

                        # labs(   title = "Detected circular RNAs per million unique mapped reads",
                        #         subtitle = "plotting circRNA predictions from DCC and uniquely mapped reads from STAR") +
                        labs(y = "") +
                        labs(x = "") +
                        # labs(fill = "Group") +
                        # labs(caption = paste(   "based on data from ",
                        #                         length(star_runs),
                        #                         " mapped libraries\ncreated: ",
                        #                         date(),
                        #                         "",
                        #                         sep="")) +
                        theme(  plot.title = element_text(lineheight=0.8,
                                face=quote(bold)),
                                legend.justification = c(1, 1),
                                legend.position = "none",
                                axis.text.x = element_text(angle = 45, hjust = 1, size=14)

                        )
                        # geom_label_repel(   data=circle_ratio,
                        #                     aes(x=name, y=num,
                        #                     label=num),
                        #                     box.padding = unit(0.0, "lines"),
                        #                     point.padding = unit(0.0, "lines")
                        #                 )
    print(circ_simple_plot)



# disable "NULL" device message
invisible(capture.output(dev.off()))

# get rif of Rplots.pdf file
invisible(capture.output(file.remove(list.files(pattern = "Rplots.pdf"))))

message("Done")







#     ######################### page one
#
#     # create data frame for ggplot2
#     raw_counts <- data.frame(circ_counts_summed, linear_counts_summed, colors)
#     colnames(raw_counts) <- c("circ","lin","group")
#
#     page_one <- ggplot( raw_counts, aes(x=circ, y=lin, color=as.factor(group), label=rownames(raw_counts))) +
#                         geom_point() +
#                         geom_label_repel(   data=raw_counts,
#                                             aes(x=circ, y=lin,
#                                             color=factor(group),
#                                             label=rownames(raw_counts)),
#                                             box.padding = unit(0.55, "lines"),
#                                             point.padding = unit(0.5, "lines")
#                                         ) +
#                         labs(   title = "Circular vs. linear read counts throughout selected libraries",
#                                 subtitle = "plotting non-normalized raw read counts from STAR") +
#                         labs(x = "Circular RNA read count (log scale)") +
#                         labs(y = "Linear RNA read count (log scale)") +
#                         labs(colour = "Group") +
#                         labs(caption = paste(   "based on data from ",
#                                                 length(star_runs),
#                                                 " mapped libraries\ncreated: ",
#                                                 date(),
#                                                 "",
#                                                 sep="")) +
#                         theme(  plot.title = element_text(lineheight=0.8,
#                                 face=quote(bold)),
#                                 # legend.justification = c(1, 1),
#                                 # legend.position = c(0.95, 0.95)
#                                 legend.justification = "top"
#
#                         ) +
#                         scale_x_log10() +
#                         scale_y_log10() +
#                         annotation_logticks(sides = "trbl")
#     print(page_one)
#
#     ######################### page two
#
#     # set number of circles
#     number_of_circles <- apply(CircRNACount[, - c(1 : 3)], 2, function(x){length(which(x > 1))})
#
#     # reset names
#     names(uniquely_mapped_reads) <- names(number_of_circles)
#
#     # create data frame for ggplot2
#     circle_counts <- data.frame(uniquely_mapped_reads, number_of_circles, colors)
#     colnames(circle_counts) <- c("unique","circles","group")
#
#     page_two <- ggplot( circle_counts, aes(x=unique, y=circles, color=as.factor(group), label=rownames(raw_counts))) +
#                         geom_point() +
#                         geom_label_repel(   data=circle_counts,
#                                             aes(x=unique, y=circles,
#                                             color=factor(group),
#                                             label=rownames(raw_counts)),
#                                             box.padding = unit(0.55, "lines"),
#                                             point.padding = unit(0.5, "lines")
#                                         ) +
#                         labs(   title = "Number of mapped reads vs. number of detected circles per library",
#                                 subtitle = "plotting uniquely mapped reads from STAR and circRNA predictions from DCC") +
#                         labs(x = "Number of mapped reads (log scale)") +
#                         labs(y = "Number of detected circles (log scale)") +
#                         labs(colour = "Group") +
#                         labs(caption = paste(   "based on data from ",
#                                                 length(star_runs),
#                                                 " mapped libraries\ncreated: ",
#                                                 date(),
#                                                 "",
#                                                 sep="")) +
#                         theme(  plot.title = element_text(lineheight=0.8,
#                                 face=quote(bold)),
#                                 # legend.justification = c(1, 1),
#                                 # legend.position = c(0.95, 0.95)
#                                 legend.justification = "top"
#                         ) +
#                         scale_x_log10() +
#                         scale_y_log10() +
#                         annotation_logticks(sides = "trbl")
#     print(page_two)
#
