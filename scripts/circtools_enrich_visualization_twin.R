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

arg_data_file <- args[1] # path is string
arg_data_file2 <- args[2] # path is string
arg_pval <- as.numeric(args[3])
arg_min_circRNAs <- as.integer(args[4]) # path is string
arg_min_rbps <- as.integer(args[5]) # path is string
arg_output <- args[6] # path is string
arg_label1 <- args[7] # path is string
arg_label2 <- args[8] # path is string



if (is.na(arg_output)){
    arg_output <- "./output.pdf"
}

if (is.na(arg_label1)){
    arg_label1 <- "default"
}

if (is.na(arg_label2)){
    arg_label2 <- "default"
}

message("Loading data file")
rbp_data <- read.delim(arg_data_file, check.names=FALSE, header = F)
rbp_data2 <- read.delim(arg_data_file2, check.names=FALSE, header = F)


colnames(rbp_data) <- c("RBP", "Annotation", "chr", "start", "stop", "strand", "p_val_circular", "raw_count_circ_rna", "observed_input_peaks_circ_rna", "length_circ_rna", "length_normalized_count_circ_rna", "number_of_features_intersecting_circ", "circ_rna_confidence_interval_0.05", "p_val_linear", "raw_count_host_gene", "observed_input_peaks_host_gene", "length_host_gene_without_circ_rna", "length_normalized_count_host_gene", "number_of_features_intersecting_linear", "host_gene_confidence_interval_0.05", "distance_normalized_counts")
colnames(rbp_data2) <- colnames(rbp_data)


# head(rbp_data)
message("plotting data")

rbp_data <- rbp_data[rbp_data$p_val_circular < arg_pval & rbp_data$p_val_linear > arg_pval, ]
rbp_data2 <- rbp_data2[rbp_data2$p_val_circular < arg_pval & rbp_data2$p_val_linear > arg_pval, ]

# we use PDF and standard A4 page size
pdf(arg_output, height= 8.2, width=11.69 , title=paste("circtools RBP enrichment analysis - ",arg_label1))

    # mixed plot

    annotation <- table(rbp_data$Annotation)
    annotation <- sort(annotation, decreasing = TRUE)
    annotation_df <- data.frame(annotation)
    colnames(annotation_df) <-c("Annotation","Frequency")
    annotation_df <- annotation_df[annotation_df$Frequency>arg_min_rbps,]

    annotation2 <- table(rbp_data2$Annotation)
    annotation2 <- sort(annotation2, decreasing = TRUE)
    annotation_df2 <- data.frame(annotation2)
    colnames(annotation_df2) <-c("Annotation","Frequency")
    annotation_df2 <- annotation_df2[annotation_df2$Frequency>arg_min_rbps,]


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
    circRNA_df_m <- melt(circRNA_df,id.vars = "CircRNA")
    circRNA_df_m$CircRNA <- factor(circRNA_df_m$CircRNA, levels = circRNA_df_m$CircRNA)


    heading <- unname(transpose(data.frame((table(rbp_data2[rbp_data2$Annotation,]$RBP))))[1,])
    heading <- unlist(heading)

    tmp_list <- list()
    counter <- 1

    for (circRNA in (rownames(annotation2))){
        row <- as.numeric(unname(transpose(data.frame((table(rbp_data2[rbp_data2$Annotation==circRNA,]$RBP)))))[2,])
        tmp_list[[counter]] <- unlist(row)
        counter <- counter+1
    }
    message(paste(counter,"circRNAs processed"))
    circRNA_df <- as.data.frame(do.call("rbind", tmp_list))
    circRNA_df <- data.frame(rownames(annotation2), circRNA_df)
    colnames(circRNA_df) <- c("CircRNA",heading)
    circRNA_df_m2 <- melt(circRNA_df,id.vars = "CircRNA")
    circRNA_df_m2$CircRNA <- factor(circRNA_df_m2$CircRNA, levels = circRNA_df_m2$CircRNA)


    total <- merge(circRNA_df_m,circRNA_df_m2,by=c("CircRNA","variable"))
    total$sum <- total$value.x + total$value.y

    for (circRNA in (unique(total$CircRNA))){
        total$all[total$CircRNA==circRNA]<-sum(total$sum[total$CircRNA==circRNA])
        total$a_max[total$CircRNA==circRNA]<-sum(total$value.x[total$CircRNA==circRNA])
        total$b_max[total$CircRNA==circRNA]<-sum(total$value.y[total$CircRNA==circRNA])
    }

    colnames(total) <- c("circRNA","RBP","A","B","rbp_sum","total_sum")

    num_rpbs <- length(unique(total$RBP))
    num_circs <- length(unique(total$circRNA))

    # only graph the top arg_min_circRNAs circRNAs
    total <- total[which(total$total_sum>arg_min_circRNAs),]

    circ_simple_plot <- ggplot(data=total) +
                        geom_bar(aes(x=reorder(circRNA, -rbp_sum, function(x) {sum(x)}), y = A, fill = RBP), stat="identity", size=0.1, colour="black") +
                        geom_bar(aes(x=reorder(circRNA, -rbp_sum, function(x) {sum(x)}), y = -B, fill = RBP), stat="identity", size=0.1, colour="black") +
                        labs(   title = "Assignment of circRNAs to RPBs",
                                subtitle = paste("plotting colour-coded RBPs per circRNA, ordered by number of total RPB sites ( p <",
                                            arg_pval, ")")) +
                        labs(y = "Number of enriched RBP binding sites") +
                        labs(x = "CircRNA") +
                        labs(caption = paste(   "based on data from ",
                                                num_circs,
                                                " circRNAs and ",
                                                num_rpbs,
                                                " RBPs, showing top ",
                                                arg_min_circRNAs,
                                                " enriched CircRNAs: ",
                                                date(),
                                                "",
                                                sep="")) +
                        scale_y_continuous(labels = commapos) +
                        theme(  plot.title = element_text(lineheight=0.8,
                                face=quote(bold)),
                                legend.justification = c(1, 1),
                                #legend.position = "none",
                                axis.text.x = element_text(angle = 45, hjust = 1, size=10)
                        ) +
                        geom_hline(yintercept=0, color = "black", size=0.5) +
                        geom_label(aes(x = 2, y = 10, label = arg_label1), fill = "white") +
                        geom_label(aes(x = 2, y = -10, label = arg_label2), fill = "white") +
                        guides(fill=guide_legend(ncol=2))

    print(circ_simple_plot)

######################### RBP -> circRNA


    # heading <- unname(transpose(data.frame((table(rbp_data$Annotation))))[1,])
    # heading <- unlist(heading)
    #
    # tmp_list <- list()
    # counter <- 1
    #
    # for (rbp in unique(rbp_data$RBP)){
    #     row <- as.numeric(unname(transpose(data.frame((table(rbp_data[rbp_data$RBP==rbp,]$Annotation)))))[2,])
    #     tmp_list[[counter]] <- unlist(row)
    #     counter <- counter+1
    # }
    # message(paste(counter,"RBPs processed"))
    # circRNA_df <- as.data.frame(do.call("rbind", tmp_list))
    #
    #
    #
    # circRNA_df <- data.frame(unique(rbp_data$RBP), circRNA_df)
    # colnames(circRNA_df) <- c("RBP",heading)
    #
    # circRNA_df_m <- melt(circRNA_df,id.vars = "RBP")
    # circRNA_df_m$RBP <- factor(circRNA_df_m$RBP, levels = circRNA_df_m$RBP)
    #
    #
    #
    # heading <- unname(transpose(data.frame((table(rbp_data2$Annotation))))[1,])
    # heading <- unlist(heading)
    #
    # tmp_list <- list()
    # counter <- 1
    #
    # for (rbp in unique(rbp_data2$RBP)){
    #     row <- as.numeric(unname(transpose(data.frame((table(rbp_data2[rbp_data2$RBP==rbp,]$Annotation)))))[2,])
    #     tmp_list[[counter]] <- unlist(row)
    #     counter <- counter+1
    # }
    #
    #
    # message(paste(counter,"RBPs processed"))
    # circRNA_df2 <- as.data.frame(do.call("rbind", tmp_list))
    # circRNA_df2 <- data.frame(unique(rbp_data2$RBP), circRNA_df2)
    # colnames(circRNA_df2) <- c("RBP",heading)
    # circRNA_df_m2 <- melt(circRNA_df2,id.vars = "RBP")
    # circRNA_df_m2$RBP <- factor(circRNA_df_m2$RBP, levels = circRNA_df_m2$RBP)
    #
    #
    #
    # total <- merge(circRNA_df_m,circRNA_df_m2,by=c("RBP","variable"))
    # total$sum <- total$value.x + total$value.y
    #
    #
    # for (rbp in (unique(total$RBP))){
    #     total$all[total$RBP==rbp]<-sum(total$sum[total$RBP==rbp])
    # }
    #
    # colnames(total) <- c("RBP","circRNA","A","B","rbp_sum","total_sum")
    #
    # num_rpbs <- length(unique(total$RBP))
    # num_circs <- length(unique(total$circRNA))
    #
    # # only graph the top arg_min_circRNAs circRNAs
    # total <- total[which(total$total_sum>arg_min_rbps),]
    #
    # nrow(total)
    #
    # circ_simple_plot <- ggplot(data=total) +
    #                     geom_bar(aes(x=reorder(RBP, -rbp_sum, function(x) {sum(x)}), y = A, fill = circRNA), stat="identity") +
    #                     geom_bar(aes(x=reorder(RBP, -rbp_sum, function(x) {sum(x)}), y = -B, fill = circRNA), stat="identity") +
    #                     labs(   title = "Assignment of RBPs to circRNAs",
    #                             subtitle = paste("plotting colour-coded circRNAs per RBP, ordered by number of total circRNAs ( p <",
    #                                         arg_pval, ")")) +
    #                     labs(y = "Number of enriched circRNAs") +
    #                     labs(x = "RNA binding protein") +
    #                     labs(caption = paste(   "based on data from ",
    #                                             num_circs,
    #                                             " circRNAs and ",
    #                                             num_rpbs,
    #                                             " RBPs, showing top ",
    #                                             arg_min_circRNAs,
    #                                             " enriched RBPs: ",
    #                                             date(),
    #                                             "",
    #                                             sep="")) +
    #                     scale_y_continuous(labels = commapos) +
    #                     theme(  plot.title = element_text(lineheight=0.8,
    #                             face=quote(bold)),
    #                             legend.justification = c(1, 1),
    #                             legend.position = "none",
    #                             axis.text.x = element_text(angle = 90, hjust = 1, size=14)
    #                     ) +
    #                     geom_hline(yintercept=0, color = "black", size=0.5) +
    #                     geom_label(aes(x = 10, y = 30, label = arg_label1), fill = "white") +
    #                     geom_label(aes(x = 10, y = -30, label = arg_label2), fill = "white")
    #
    # print(circ_simple_plot)
    #
    #



# disable "NULL" device message
invisible(capture.output(dev.off()))

# get rif of Rplots.pdf file
invisible(capture.output(file.remove(list.files(pattern = "Rplots.pdf"))))

message("Done")
