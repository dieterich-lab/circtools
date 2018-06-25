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
arg_colour_mode <- args[7]
arg_use_only_circ_data <- args[8]


# arguments for file 2 are ath the send so we can leave them empty
arg_label_sample_2 <- args[9]
arg_data_file_2 <- args[10]


# check circ-only mode
if (arg_use_only_circ_data != "True" & arg_use_only_circ_data != "False" ) {
    print(arg_use_only_circ_data)
    message("Please specify the data mode")
    quit()
}

# check colour mode
if (arg_colour_mode != "colour" & arg_colour_mode != "bw" ) {
    print(arg_colour_mode)
    message("Please specify the colour model: colour or bw")
    quit()
}

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


if (arg_use_only_circ_data == "True" ) {
    rbp_data_file_1 <- rbp_data_file_1[rbp_data_file_1$p_val_circular < arg_pval,]
} else {
    rbp_data_file_1 <- rbp_data_file_1[rbp_data_file_1$p_val_circular < arg_pval & rbp_data_file_1$p_val_linear > arg_pval,]
}

if (!is.na(arg_data_file_2)) {
    message("Reading input file 2")
    rbp_data_file_2 <- read.delim(arg_data_file_2, check.names = FALSE, header = F)
    colnames(rbp_data_file_2) <- colnames(rbp_data_file_1)

    if (arg_use_only_circ_data == "True" ) {
        rbp_data_file_2 <- rbp_data_file_2[rbp_data_file_2$p_val_circular < arg_pval,]
    } else {
        rbp_data_file_2 <- rbp_data_file_2[rbp_data_file_2$p_val_circular < arg_pval & rbp_data_file_2$p_val_linear > arg_pval,]
    }

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

print_isoforms <- 1

if (length(rbp_df) == 1){
    rbp_df$RBP <- rownames(rbp_df)
    colnames(rbp_df) <- c("A", "B")
    rbp_df <- rbp_df[,c("B","A")]
    rownames(rbp_df) <- NULL
    print_isoforms <- 0
}

colnames(rbp_df) <- c("RBP", "Frequency")
rbp_df <- head(rbp_df, arg_max_rbps)

if (!is.na(arg_data_file_2)) {
    tmp <- data.frame(rbp_data_file_2$RBP)
    colnames(tmp) <- c("RBP")
    tmp$annotation <- rbp_data_file_2$Annotation
    tmp <- unique(tmp)
    rbps2 <- table(tmp$RBP)

    rbps2 <- sort(rbps2, decreasing = TRUE)
    rbp_df2 <- data.frame(rbps2)

    if (length(rbp_df2) == 1){
        rbp_df2$RBP <- rownames(rbp_df2)
        colnames(rbp_df2) <- c("A", "B")
        rbp_df2 <- rbp_df2[,c("B","A")]
        rownames(rbp_df2) <- NULL
        print_isoforms <- 0
    }

    colnames(rbp_df2) <- c("RBP", "Frequency")
    rbp_df2 <-  head(rbp_df2, arg_max_rbps)
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


num_rbps <- length(unique(total$RBP))


# limit to top X
total <- head(total[order(-total$FrequencyA),] ,arg_max_rbps)


# start actual plotting
rbp_simple_plot <- ggplot(data=total) +
                        geom_bar(stat="identity", colour="black", size=0.1, aes(x=reorder(RBP, -FrequencyA), y=FrequencyA, fill=RBP)) +

                        geom_label(aes(x = arg_max_rbps - 1 , y = label_pos_1, label = arg_label_sample_1), fill = "white")

                        if (!is.na(arg_data_file_2)) {
                            rbp_simple_plot <- rbp_simple_plot + geom_bar(stat="identity", colour="black", size=0.1, aes(x=reorder(RBP, -FrequencyA), y=-FrequencyB, fill=RBP))
                            rbp_simple_plot <- rbp_simple_plot + geom_label(aes(x = arg_max_rbps - 1, y = - label_pos_2, label = arg_label_sample_2), fill = "white")

                        }
                        if (arg_colour_mode == "bw" ) {
                            rbp_simple_plot <- rbp_simple_plot + scale_fill_grey(start = 0, end = .9)
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
                                                " RBPs: ",
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
### We do not accumulate isoforms and only count unique hits
### This corresponds to the Top X single circRNAs plotted before
########################################################################################################################


if (!is.na(arg_data_file_2)) {
    tmp1 <- data.frame(table(rbp_data_file_1$Annotation))
    tmp2 <- data.frame(table(rbp_data_file_2$Annotation))
    circ_rna_selection <- merge(tmp1, tmp2, by = c("Var1"), all = TRUE)
    circ_rna_selection[is.na(circ_rna_selection)] <- 0
    circ_rna_selection$Freq <- circ_rna_selection$Freq.x + circ_rna_selection$Freq.y
    colnames(circ_rna_selection) <- c( "circRNA", "A", "B", "Total")

} else {
    circ_rna_selection <- data.frame(table(rbp_data_file_1$Annotation))
    colnames(circ_rna_selection) <- c( "circRNA", "A")

}

num_circs <- length(unique(circ_rna_selection$circRNA))

if (!is.na(arg_data_file_2)) {
    circ_rna_selection$Var1 <- levels(droplevels(circ_rna_selection$circRNA))
    circ_rna_selection <- circ_rna_selection[order(- circ_rna_selection$A),]
    circ_rna_selection <- circ_rna_selection[circ_rna_selection$Total > 0,]
    selected_circrnas <- head(circ_rna_selection$circRNA, arg_max_circRNAs)

} else {
    circ_rna_selection$Var1 <- levels(droplevels(circ_rna_selection$circRNA))
    circ_rna_selection <- circ_rna_selection[order(- circ_rna_selection$A),]
    circ_rna_selection <- circ_rna_selection[circ_rna_selection$A > 0,]
    selected_circrnas <- head(circ_rna_selection$Var1, arg_max_circRNAs)
}


circ_rna_selection <- circ_rna_selection[circ_rna_selection$Var1 %in% selected_circrnas ,]

label_pos <- (max(circ_rna_selection$A, na.rm = TRUE))

circ_simple_plot <- ggplot(data = circ_rna_selection) +
    geom_bar(aes(x = reorder(circRNA, - A), y = A, fill = circRNA), stat = "identity", size = 0.0, colour = "black") +

    geom_label(aes(x = arg_max_circRNAs - 4 , y = label_pos, label = arg_label_sample_1), fill = "white")

    if (!is.na(arg_data_file_2)) {
        circ_simple_plot <- circ_simple_plot + geom_bar(aes(x = reorder(circRNA, - Total), y = -B, fill = circRNA), stat = "identity", size = 0.0, colour = "black") +               geom_label(aes(x = arg_max_circRNAs - 4, y = - label_pos, label = arg_label_sample_2), fill = "white")
    }

    # geom_bar(aes(x = reorder(circRNA, - total_sum), y = - B, fill = RBP), , stat = "identity", size = 0.0, colour = "black") +
    circ_simple_plot <- circ_simple_plot +  labs(title = paste(arg_label_sample_1, ": Top CircRNAs (by RBP hits)", sep=""),
    subtitle = paste("plotting colour-coded RBPs per circRNA, ordered by number of (distinct) RBP hits ( p <",
    arg_pval, ")")) +
    labs(y = "Number of different RBPs found for circRNA") +
    labs(x = "CircRNA") +
    labs(caption = paste("based on data from ",
    num_circs,
    " circRNAs and ",
    num_rbps,
    " RBPs, showing top ",
    arg_max_circRNAs,
    " circRNAs: ",
    date(),
    "",
    sep = "")) +
    scale_y_continuous(labels = commapos)  #, limits = c(-50,50)

    if (arg_colour_mode == "bw" ) {
        circ_simple_plot <- circ_simple_plot + scale_fill_grey(start = 0, end = .9, name = "CircRNAs")
    } else {
        circ_simple_plot <- circ_simple_plot + scale_fill_discrete(name = "CircRNAs")
    }

    circ_simple_plot <- circ_simple_plot +
    theme(plot.title = element_text(lineheight = 0.8,
    face = quote(bold)),
    legend.justification = c(1, 1),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.3, "cm"),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10)
    ) +
    geom_hline(yintercept = 0, color = "black", size = 0.5) +
    guides(fill = guide_legend(ncol = 2))

print(circ_simple_plot)


########################################################################################################################
### This is the # circRNAs per RBP plot
### We accumulate isoforms!
########################################################################################################################

# only plot if we see more than one RBP
if (print_isoforms == 1){

    if (!is.na(arg_data_file_2)) {
        sample_list <- list(rbp_data_file_1, rbp_data_file_2)
        sample_names <- list(arg_label_sample_1, arg_label_sample_2)
    } else {
        sample_list <- list(rbp_data_file_1)
        sample_names <- list(arg_label_sample_1)
    }

    for (sample in seq(1 : length(sample_list))) {

        current_data <- sample_list[[sample]]
        sub_dataframe <- data.frame(table(current_data$Annotation))
        sub_dataframe$Var1 <- levels(droplevels(sub_dataframe$Var1))
        sub_dataframe <- sub_dataframe[order(- sub_dataframe$Freq),]

        if (nrow(sub_dataframe) < arg_max_circRNAs){
            arg_max_circRNAs <- nrow(sub_dataframe)
        }

        selected_circrnas <- sub_dataframe$Var1[1 : arg_max_circRNAs]

        for (top_circ in selected_circrnas) {


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
                tmp_frame <- head(tmp_frame, arg_max_rbps)
                miniframe <- data.frame(tmp_frame$RBP, tmp_frame$clip_peaks)
                colnames(miniframe) <- c("RBP", "clip_peaks")
                    miniframe <- ddply(unique(miniframe), .(RBP), transform, border = rep(1, clip_peaks))

                    theme_set(theme_grey(base_size = 8))
                    circrna_plot[[circ_isoform]] <- ggplot(miniframe, aes(RBP)) +
                        geom_bar(aes(x = reorder(paste(RBP, ": ", clip_peaks, sep = ""), - clip_peaks), clip_peaks, fill = RBP), width = 1, size = 0.15, stat = "identity", color = "white")

                        if (arg_colour_mode == "bw" ) {
                            circrna_plot[[circ_isoform]] <- circrna_plot[[circ_isoform]] + scale_fill_grey(start = 0, end = .9)
                        }

                        circrna_plot[[circ_isoform]] <- circrna_plot[[circ_isoform]] +
                        scale_y_continuous() +
                        coord_polar() +
                        theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks = element_blank()) +
                        labs(title = paste(sample_names[sample], ":\nComposition of RBP landscape for circRNA", tmp_frame[1, 2]),
                        subtitle = paste("Isoform ", circ_isoform, ": Chromsome ", tmp_frame[1, 3], ", ", commapos(as.integer(tmp_frame[1, 4])), "->", commapos(as.integer(tmp_frame[1, 5])), sep = "")) +
                        labs(y = "Accumulated eCLIP peak number observed per RBP and circRNA") +
                        labs(x = "") +
                        labs(caption = paste("Top circRNAs enriched for RBP peaks compared to their host gene ( p <",
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
}
########################################################################################################################
### These are the # RBPs peaks per circRNA - i.e. the maximum of the observed peaks
### for each circRNA (non isoform specific)
### Per column we see the circRNA X with the observed CLIP peaks summed up for each RBP
########################################################################################################################

if (!is.na(arg_data_file_2)) {
    tmp1 <- data.frame(table(rbp_data_file_1$Annotation))
    tmp2 <- data.frame(table(rbp_data_file_2$Annotation))
    circ_rna_selection <- merge(tmp1, tmp2, by = c("Var1"), all = TRUE)
    circ_rna_selection[is.na(circ_rna_selection)] <- 0
    circ_rna_selection$Freq <- circ_rna_selection$Freq.x + circ_rna_selection$Freq.y
} else {
    circ_rna_selection <- data.frame(table(rbp_data_file_1$Annotation))
}

circ_rna_selection$Var1 <- levels(droplevels(circ_rna_selection$Var1))
circ_rna_selection <- circ_rna_selection[order(- circ_rna_selection$Freq),]


circ_rna_selection <- circ_rna_selection[circ_rna_selection$Freq > 0,]

selected_circrnas <- head(circ_rna_selection$Var1, arg_max_circRNAs)

data1 <- as.data.table(rbp_data_file_1[, c(1, 2, 9)])
data1 <- (unique(data1[data1[, .I[observed_input_peaks_circ_rna == max(observed_input_peaks_circ_rna)], by = list(RBP, Annotation)]$V1]))

if (!is.na(arg_data_file_2)) {
    data2 <- as.data.table(rbp_data_file_2[, c(1, 2, 9)])
    data2 <- (unique(data2[data2[, .I[observed_input_peaks_circ_rna == max(observed_input_peaks_circ_rna)], by = list(RBP, Annotation)]$V1]))
    total <- merge(data1, data2, by = c("Annotation", "RBP"), all = TRUE)
    total[is.na(total)] <- 0
    colnames(total) <- c("circRNA", "RBP", "A", "B")
    total$rbp_sum <- total$A + total$B

} else {
    total <- data1
    colnames(total) <- c( "RBP", "circRNA", "A")
    total[is.na(total)] <- 0
}

num_circs <- length(unique(total$circRNA))
num_rbps <- length(unique(total$RBP))

for (top_circ in unique(total$circRNA)) {
    if (!is.na(arg_data_file_2)) {
        total$total_sum[total$circRNA == top_circ] <- sum(total[total$circRNA == top_circ , rbp_sum])
    }
    total$total_sum_A[total$circRNA == top_circ] <- sum(total[total$circRNA == top_circ , A])

    if (!is.na(arg_data_file_2)) {
        total$total_sum_B[total$circRNA == top_circ] <- sum(total[total$circRNA == top_circ , B])
    }
}

selected_circrnas <- as.character((head(unique(total[order(- total_sum_A),circRNA]),arg_max_circRNAs)))
total <- total[total$circRNA %in% selected_circrnas ,]

if (!is.na(arg_data_file_2)) {
    label_pos <- (max(total$total_sum, na.rm = TRUE) / 2) - 50
} else {
    label_pos <- (max(circ_rna_selection$A, na.rm = TRUE))
}
circ_simple_plot <- ggplot(data = total) +
    geom_bar(aes(x = reorder(circRNA, - total_sum_A), y = A, fill = RBP), stat = "identity", size = 0.0, colour = "black") +
    geom_label(aes(x = arg_max_circRNAs - 4 , y = label_pos, label = arg_label_sample_1), fill = "white")

    if (!is.na(arg_data_file_2)) {
        circ_simple_plot <- circ_simple_plot + geom_bar(aes(x = reorder(circRNA, - total_sum_A), y = - B, fill = RBP), stat = "identity", size = 0.0, colour = "black") +         geom_label(aes(x = arg_max_circRNAs - 4, y = - label_pos, label = arg_label_sample_2), fill = "white")
    }

    # geom_bar(aes(x = reorder(circRNA, - total_sum), y = - B, fill = RBP), , stat = "identity", size = 0.0, colour = "black") +
    circ_simple_plot <- circ_simple_plot +  labs(title = paste(arg_label_sample_1, ": Assignment of RBP CLIP peaks to circRNAs", sep=""),

    subtitle = paste("plotting colour-coded RBPs per circRNA, ordered by number of accumulated RBP eCLIP peaks ( p <",
    arg_pval, ")")) +
    labs(y = "Number of enriched RBP binding sites detected in the CLIP data set") +
    labs(x = "CircRNA") +
    labs(caption = paste("based on data from ",
    num_circs,
    " circRNAs and ",
    num_rbps,
    " RBPs, showing top ",
    arg_max_circRNAs,
    " circRNAs: ",
    date(),
    "",
    sep = "")) +
    scale_y_continuous(labels = commapos)  #, limits = c(-50,50)

    if (arg_colour_mode == "bw" ) {
        circ_simple_plot <- circ_simple_plot + scale_fill_grey(start = 0, end = .9, name = "RBPs")
    } else {
        circ_simple_plot <- circ_simple_plot + scale_fill_discrete(name = "RBPs")
    }

    circ_simple_plot <- circ_simple_plot +
    theme(plot.title = element_text(lineheight = 0.8,
    face = quote(bold)),
    legend.justification = c(1, 1),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.3, "cm"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10)
    ) +
    geom_hline(yintercept = 0, color = "black", size = 0.5) +
    guides(fill = guide_legend(ncol = 2))

print(circ_simple_plot)


########################################################################################################################
### Done with plotting, finishing up
########################################################################################################################

message(paste("Printing to",arg_output_file_name))

# disable "NULL" device message
invisible(capture.output(dev.off()))

# get rif of Rplots.pdf file
invisible(capture.output(file.remove(list.files(pattern = "Rplots.pdf"))))

message("Done")
