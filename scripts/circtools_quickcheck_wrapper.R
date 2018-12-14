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

# utiliy functions

getUniqueMappings <- function(STARfolder)
{
    # Read STAR log file
    logFile <- paste(STARfolder, "Log.final.out", sep = "/")

    # extract unique mappings
    unique_mapped <- read.delim(logFile, as.is = T, sep = "|")
    unique_mapped <- as.numeric(gsub("\t", "", unique_mapped[7, 2]))
    return(unique_mapped)
}

# may not be optimal, but we do not want warning for now
options(warn = - 1)

# quiet mode
options(echo = FALSE)

# read arguments
args <- commandArgs(trailingOnly = TRUE)

arg_dcc_data <- args[1] # path is string
arg_star_folder <- args[2] # path is string
arg_output_directory <- args[3] # string
arg_condition_list <- strsplit(args[4], ",")[[1]] # list of strings
arg_grouping <- unlist(lapply(strsplit(args[5], ","), as.numeric)) # list of strings
arg_colour_mode <- args[6]
arg_cleanup_string <- args[7]
arg_starfolder_suffix <- args[8]
arg_remove_last <- as.numeric(args[9])
arg_remove_first <- as.numeric(args[10])
arg_remove_columns <- unlist(lapply(strsplit(args[11], ","), as.numeric))

# check colour mode
if (arg_colour_mode != "colour" & arg_colour_mode != "bw" ) {
    print(arg_colour_mode)
    message("Please specify the colour model: colour or bw")
    quit()
}

group_length <- length(arg_condition_list)

## load complete data set
message("Loading CircRNACount")
CircRNACount <- read.delim(paste(arg_dcc_data, "CircRNACount", sep="/"), check.names=FALSE, header = T)

message("Loading LinearRNACount")
LinearCount <- read.delim(paste(arg_dcc_data, "LinearCount", sep="/"), check.names=FALSE, header = T)

# check columns to remove
if (arg_remove_columns != "0" && length(arg_remove_columns) > 0) {

    LinearCount <- LinearCount[, -arg_remove_columns]
    CircRNACount <- CircRNACount[, -arg_remove_columns]
}

message("Parsing data")
# remove DCC artifacts from columns names
names(LinearCount)<-gsub(arg_cleanup_string,"",names(LinearCount))
names(CircRNACount)<-gsub(arg_cleanup_string,"",names(CircRNACount))



min_length = nchar(names(LinearCount[which.min(unlist(lapply(names(LinearCount[, - c(1 : 3)]), nchar)))+3]))

if (min_length - arg_remove_first - arg_remove_last < 0){
    message("Values specified via -L and -F would result in a least one too short sample name.
    Please reduce the cuttoff values.")
    q()
}


names(LinearCount)<-gsub(arg_cleanup_string,"",names(LinearCount))
names(CircRNACount)<-gsub(arg_cleanup_string,"",names(CircRNACount))

names(LinearCount) <- c(names(LinearCount[, c(1 : 3)]),
strtrim(names(LinearCount[, - c(1 : 3)]), nchar(names(LinearCount[, - c(1 : 3)]))-arg_remove_last))

names(CircRNACount) <- c(names(LinearCount[, c(1 : 3)]),
strtrim(names(CircRNACount[, - c(1 : 3)]), nchar(names(CircRNACount[, - c(1 : 3)]))-arg_remove_last))

names(LinearCount) <- c(names(LinearCount[, c(1 : 3)]),
substring(names(LinearCount[, - c(1 : 3)]), arg_remove_first))

names(CircRNACount) <- c(names(LinearCount[, c(1 : 3)]),
substring(names(CircRNACount[, - c(1 : 3)]), arg_remove_first))

# summing up counts per column (i.e. library)
circ_counts_summed <- apply(CircRNACount[, - c(1 : 3)], 2, sum)
linear_counts_summed <- apply(LinearCount[, - c(1 : 3)], 2, sum)

num_samples <- ncol(LinearCount[, - c(1 : 3)])

message(paste("Found ", num_samples, " data columns in provided DCC data", sep=""))

message(paste(group_length, " different groups provided", sep=""))

if (length(arg_grouping) < num_samples) {
  message("Assuming (1,2),(1,2),(1,2),... sample grouping")
  dummy_list <- rep(arg_grouping, (num_samples/group_length))
  colors <- unlist(lapply(seq(1, num_samples), function(x) {
    return(arg_condition_list[dummy_list[x]])
  }))
} else {
  message("Setting sample groups manually")
  colors <- unlist(lapply(seq(1, num_samples), function(x) {
    return(arg_condition_list[arg_grouping[x]])
  }))
    counts <- c(rep(0,num_samples))
  generated_sample_names <- unlist(lapply(seq(1, num_samples), function(x) {
      counts[arg_grouping[x]] <<- counts[arg_grouping[x]]+1
    return(paste(arg_condition_list[arg_grouping[x]], counts[arg_grouping[x]], sep = " "))
  }))
}
# names(circ_counts_summed) <- generated_sample_names
# names(linear_counts_summed) <- generated_sample_names

# get unique mapping reads
## which star runs are in the DCC output?
star_columns <- colnames(CircRNACount[, - c(1 : 3)])

# get paths for those directories
star_runs <- list.dirs(arg_star_folder, full.names = TRUE)
star_runs <- star_runs[-1]

# only use the folders ending with _STARmapping (Dieterich lab default)
if (arg_starfolder_suffix != 0){
    star_runs <- star_runs[endsWith(star_runs, arg_starfolder_suffix)]
}
# only use the folders of full mappings (not the per-mate ones)
star_runs <- star_runs[!grepl("*mate*", star_runs)]

# check columns to remove
if (arg_remove_columns != "0" && length(arg_remove_columns) > 0) {
    tmp <- unlist(lapply(arg_remove_columns, function(x){x-3}))
    star_runs <- star_runs[-tmp]
}

# new empty list
uniquely_mapped_reads <- numeric();

for (run in star_runs)
{
    uniquely_mapped_reads <- c(uniquely_mapped_reads, getUniqueMappings(run));
}

# set names for unique reads
names(uniquely_mapped_reads) <- names(circ_counts_summed)

message("plotting data")

# we use PDF and standard A4 page size
pdf(paste(arg_output_directory, ".pdf", sep = "") , height= 8.2, width=11.69 , title="circtools library analysis")

    ######################### page one

    # create data frame for ggplot2
    raw_counts <- data.frame(circ_counts_summed, linear_counts_summed, colors)
    colnames(raw_counts) <- c("circ","lin","group")

    # rownames(raw_counts) <- generated_sample_names
    # rownames(raw_counts) <- generated_sample_names

    if (arg_colour_mode == "bw" ) {
        page_one <- ggplot( raw_counts, aes(x=circ, y=lin, shape = as.factor(group) , label=rownames(raw_counts))) +
                        geom_point(size = 2) + scale_color_grey(start = 0, end = 0) + labs(shape = "Group")
    } else {
        page_one <- ggplot( raw_counts, aes(x=circ, y=lin, color=as.factor(group), label=rownames(raw_counts))) +
                        geom_point() + labs(color = "Group")
    }

                        page_one <- page_one +
                        geom_label_repel(   data=raw_counts,
                                            aes(x=circ, y=lin,
                                            #color=factor(group),
                                            label=rownames(raw_counts)),
                                            box.padding = unit(0.55, "lines"),
                                            point.padding = unit(0.5, "lines")
                                        ) +
                        labs(   title = "Circular vs. linear read counts throughout selected libraries",
                                subtitle = "plotting non-normalized raw read counts from STAR") +
                        labs(x = "Circular RNA read count (log scale)") +
                        labs(y = "Linear RNA read count (log scale)") +
                        labs(caption = paste(   "based on data from ",
                                                length(star_runs),
                                                " mapped libraries\ncreated: ",
                                                date(),
                                                "",
                                                sep="")) +
                        theme(  plot.title = element_text(lineheight=0.8,
                                face=quote(bold)),
                                # legend.justification = c(1, 1),
                                # legend.position = c(0.95, 0.95)
                                legend.justification = "top"
                        ) +
                        scale_x_log10(
                            breaks = scales::trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x))
                        ) +
                        scale_y_log10(
                            breaks = scales::trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x))
                        ) +
                        annotation_logticks(sides = "trbl")

    if (arg_colour_mode == "bw" ) {
                        page_one <- page_one +  labs(colour = "Group") + labs(shape = "Group")
    } else {
                        page_one <- page_one + labs(colours = "Group")
    }

    print(page_one)

    ######################### page two

    # set number of circles
    number_of_circles <- apply(CircRNACount[, - c(1 : 3)], 2, function(x){length(which(x > 1))})

    # reset names
    names(uniquely_mapped_reads) <- names(number_of_circles)

    # create data frame for ggplot2
    circle_counts <- data.frame(uniquely_mapped_reads, number_of_circles, colors)
    colnames(circle_counts) <- c("unique","circles","group")

    if (arg_colour_mode == "bw" ) {
        page_two <- ggplot( circle_counts, aes(x=unique, y=circles, shape = as.factor(group) , label=rownames(raw_counts))) +
                        geom_point(size = 2) + scale_color_grey(start = 0, end = 0) + labs(shape = "Group")
    } else {
        page_two <- ggplot( circle_counts, aes(x=unique, y=circles, color=as.factor(group), label=rownames(raw_counts))) +
                        geom_point() + labs(colours = "Group")
    }

                        page_two <- page_two +
                        geom_label_repel(   data=circle_counts,
                                            aes(x=unique, y=circles,
                                            #color=factor(group),
                                            label=rownames(raw_counts)),
                                            box.padding = unit(0.55, "lines"),
                                            point.padding = unit(0.5, "lines")
                                        ) +
                        labs(   title = "Number of mapped reads vs. number of detected circles per library",
                                subtitle = "plotting uniquely mapped reads from STAR and circRNA predictions from DCC") +
                        labs(x = "Number of mapped reads (log scale)") +
                        labs(y = "Number of detected circles (log scale)") +
                        labs(caption = paste(   "based on data from ",
                                                length(star_runs),
                                                " mapped libraries\ncreated: ",
                                                date(),
                                                "",
                                                sep="")) +
                        theme(  plot.title = element_text(lineheight=0.8,
                                face=quote(bold)),
                                # legend.justification = c(1, 1),
                                # legend.position = c(0.95, 0.95)
                                legend.justification = "top"
                        ) +
                        scale_x_log10(
                            breaks = scales::trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x))
                        ) +
                        scale_y_log10(
                            breaks = scales::trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x))
                        ) +
                        annotation_logticks(sides = "trbl")

    if (arg_colour_mode == "bw" ) {
                        page_two <- page_two + labs(colour = "Group") + labs(shape = "Group")
    } else {
                        page_two <- page_two + labs(colour = "Group")
    }

    print(page_two)

    ######################### page three

    ref = order(number_of_circles / (uniquely_mapped_reads / 1000000))

    # create data frame for ggplot2
    circle_ratio <- data.frame(rownames(raw_counts), as.integer((number_of_circles / (uniquely_mapped_reads / 1000000))), colors)
    colnames(circle_ratio) <- c("name","num","group")
    circle_ratio <- circle_ratio[with(circle_ratio, order(num)), ]

    page_three <- ggplot(data=circle_ratio, aes(x=reorder(name, num), y=num, fill=as.factor(group))) +
                        geom_bar(stat="identity", colour="black")

                        if (arg_colour_mode == "bw" ) {
                            page_three <- page_three + scale_fill_grey(start = 0.4, end = 1)
                        }

                        page_three <- page_three +
                        labs(   title = "Detected circular RNAs per million unique mapped reads",
                                subtitle = "plotting circRNA predictions from DCC and uniquely mapped reads from STAR") +
                        labs(y = "Number of detected circular RNAs") +
                        labs(x = "Sequencing library") +
                        labs(fill = "Group") +
                        labs(caption = paste(   "based on data from ",
                                                length(star_runs),
                                                " mapped libraries\ncreated: ",
                                                date(),
                                                "",
                                                sep="")) +
                        theme(  plot.title = element_text(lineheight=0.8,
                                face=quote(bold)),
                                legend.justification = c(1, 1),
                                legend.position = c(0.25, 0.98),
                                axis.text.x = element_text(angle = 45, hjust = 1, size=14)

                        ) +
                        geom_label_repel(   data=circle_ratio,
                                            aes(x=name, y=num,
                                            label=num),
                                            box.padding = unit(0.0, "lines"),
                                            point.padding = unit(0.0, "lines")
                                        )
    print(page_three)

invisible(capture.output(dev.off()))

# get rif of Rplots.pdf file
invisible(capture.output(file.remove(list.files(pattern = "Rplots.pdf"))))

message("Done")