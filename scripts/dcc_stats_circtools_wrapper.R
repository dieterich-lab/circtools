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


arg_replictes <- as.integer(args[2]) # integer
arg_condition_list <- strsplit(args[3], ",")[[1]] # list of strings
arg_condition_columns <- lapply(strsplit(args[4], ","), as.numeric) # list of integers
arg_condition_columns <- unlist(arg_condition_columns)
arg_groups <- unlist(lapply(strsplit(args[5], ","), as.numeric)) # list of strings


## load complete data set
message("Loading CircRNACount")
CircRNACount <- read.delim(paste(arg_dcc_data, "CircRNACount", sep="/"), check.names=FALSE, header = T)

message("Loading LinearRNACount")
LinearCount <- read.delim(paste(arg_dcc_data, "LinearCount", sep="/"), check.names=FALSE, header = T)

message("Parsing data")
# remove DCC artifacts from columns names
names(LinearCount)<-gsub("_STARmappingChimeric.out.junction","",names(LinearCount))
names(CircRNACount)<-gsub("_STARmappingChimeric.out.junction","",names(CircRNACount))

# summing up counts per column (i.e. library)
circ_counts_summed <- apply(CircRNACount[, - c(1 : 3)], 2, sum)
linear_counts_summed <- apply(LinearCount[, - c(1 : 3)], 2, sum)

# get unique mapping reads
## which star runs are in the DCC output?
star_columns <- colnames(CircRNACount[, - c(1 : 3)])

# get paths for those directories
star_runs <- list.files(arg_star_folder, full.names = TRUE)

# only use the folders ending with _STARmapping (Dieterich lab default)
star_runs <- star_runs[endsWith(star_runs,"_STARmapping")]

# only use the folders of full mappings (not the per-mate ones)
star_runs <- star_runs[-grep(c("mate"), star_runs)]

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
pdf(paste(arg_output_directory, ".pdf", sep = "") , paper="a4", title="circtools library analysis")

# "page" one: circular RNA read count plottet against linear read count

    plot(   x = circ_counts_summed,
            y = linear_counts_summed,
            log = "xy",
            xlab = "CircReadCount",
            ylab = "LinearReadCount",
            col = rep(c("blue", "red"), 10),
            main = "Circular vs. linear read counts throughout selected libraries"
    )

    text(   x = circ_counts_summed,
            y = linear_counts_summed,
            names(circ_counts_summed),
            col = "black",
            cex = 0.8,
            pos = 1,
            srt = 45
    )

    legend("bottomright", c("RNASe+", "RNASe-"), fill = c("red", "blue"))

# "page" two: Number of mapped reads vs. number of detected circles per library

    # set number of circles
    number_of_circles <- apply(CircRNACount[, - c(1 : 3)], 2, function(x){length(which(x > 1))})

    # reset names
    names(uniquely_mapped_reads) <- names(number_of_circles)

    plot(   x = uniquely_mapped_reads,
            y = number_of_circles,
            log = "xy",
            ylab = "Number of detected circles",
            xlab = "Number of mapped reads",
            col = rep(c("blue", "red"), 9),
            main = "Number of mapped reads vs. number of detected circles per library"
    )

    text(   x = uniquely_mapped_reads,
            y = number_of_circles,
            names(uniquely_mapped_reads),
            cex = 0.5,
            col = "black",
            srt = 45
    )

    legend("topright", c("RNASe+", "RNASe-"), fill = c("red", "blue"))


# "page" two: Number of mapped reads vs. number of detected circles per library
ref = order(number_of_circles / (uniquely_mapped_reads / 1000000))
barplot(sort(number_of_circles / (uniquely_mapped_reads / 1000000)), las = 2, col = rep(c("blue", "red"), 9)[ref], main = "Circles_per_million_unique_mapper")
legend("topleft", c("RNASe+", "RNASe-"), fill = c("red", "blue"))

dev.off();

message("Done")