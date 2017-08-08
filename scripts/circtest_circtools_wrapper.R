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

library(CircTest)

# pre load libraries so we don't get messages later:
suppressMessages(library(aod))
suppressMessages(library(ggplot2))
suppressMessages(library(plyr))

# may not be optimal, but we do not want warning for now
options(warn=-1)

# fix 'could not find zip' error from openxlsx
Sys.setenv(R_ZIPCMD = "/usr/bin/zip")
library(openxlsx)

# fixed params:
MAX_LINES <- 50000

# let's now get the command line arguments

# quiet mode
options(echo=FALSE)

# read arguments
args <- commandArgs(trailingOnly = TRUE)

# assign input to script variables
arg_dcc_data <- args[1] # path is string
arg_replictes <- as.integer(args[2]) # integer
arg_condition_list <- strsplit(args[3],",")[[1]] # list of strings
arg_condition_columns <- lapply(strsplit(args[4],","), as.numeric) # list of integers
arg_condition_columns <- unlist(arg_condition_columns)
arg_output_name <- args[5] # string
arg_max_fdr <-as.numeric(args[6]) # float
arg_max_plots <- as.integer(args[7]) # string
arg_filter_sample <- as.integer(args[8]) # integer
arg_filter_count <- as.integer(args[9]) # integer
arg_groups <-   unlist(lapply(strsplit(args[10],","), as.numeric)) # list of strings
arg_output_label <- args[11] # string
arg_percent_filter <-as.numeric(args[12]) # float


run_CircTest = function(CircRNACount, LinearCount, CircCoordinates, groups, indicators, label, filename, filer.sample,
                        filter.count, max_fdr, max.plots, replicates, percent_filter) {

    message("Filtering circRNA counts")
    CircRNACount_filtered <- Circ.filter(circ = CircRNACount, linear = LinearCount,
    Nreplicates = replicates, filter.sample = filer.sample, filter.count =  filter.count, percentage = percent_filter)

    message("Filtering circRNA coordinates")
    CircCoordinates_filtered <- CircCoordinates[rownames(CircRNACount_filtered),]

    message("Filtering linear RNA counts")
    LinearCount_filtered <- LinearCount[rownames(CircRNACount_filtered),]

    message("running circTest")

    data = Circ.test(CircRNACount_filtered, LinearCount_filtered, CircCoordinates_filtered, group = groups, alpha=max_fdr)

    message("Generating plots")

    pdf(file = paste(filename,".pdf",sep=""), width= 8.2, height=11.69, title="circtools circtest analysis")

    max <- nrow(data$summary_table)

    # if we havbe too many results, cut them down to the user threshold
    if (nrow(data$summary_table) > max.plots) {
        max <- max.plots
    }

    for (i in rownames(data$summary_table[1 : max,])) {

        # invisible(capture.output(p[[i]] <- Circ.ratioplot( CircRNACount_filtered, LinearCount_filtered,
        # CircCoordinates_filtered, plotrow=i, size=16, gene_column=4, groupindicator1 = indicators,
        # x = "", y = "", lab_legend = label)))

        invisible(capture.output(Circ.ratioplot(CircRNACount_filtered, LinearCount_filtered, CircCoordinates_filtered,
        plotrow = i, size = 16, gene_column = 4, groupindicator1 = indicators,
        x = "", y = "", lab_legend = label)))
    }
    dev.off()

    message("creating Excel sheets")

    ## Significant result show in a summary table

    write.table(na.omit(data$summary_table[1 : MAX_LINES,]),
                file = paste(filename, ".csv", sep = "") ,
                quote = FALSE,
                sep = "\t",
                eol = "\n",
                row.names = FALSE,
                col.names = FALSE)

    wb <- createWorkbook()

    addWorksheet(wb, sheetName = "Significant circles")
    writeDataTable(wb, sheet = 1, x = data$summary_table[1 : MAX_LINES,])

    idx <- rownames(data$summary_table[1 : MAX_LINES,])

    addWorksheet(wb, sheetName = "Circle Counts")
    writeDataTable(wb, sheet = 2, x = CircRNACount[idx,])

    addWorksheet(wb, sheetName = "Linear Counts")
    writeDataTable(wb, sheet = 3, x = LinearCount[idx,])

    saveWorkbook(wb, paste(filename, ".xlsx", sep = "") , overwrite = TRUE)
}

###########################################################

## load complete data set
message("Loading CircRNACount")
CircRNACount <- read.delim(paste(arg_dcc_data, "CircRNACount", sep="/"), header = T)

message("Loading LinearRNACount")
LinearCount <- read.delim(paste(arg_dcc_data, "LinearCount", sep="/"), header = T)

message("Loading CircCoordinates")
CircCoordinates <- read.delim(paste(arg_dcc_data, "CircCoordinates", sep="/"), header = T)

# call the main function to run circTest

run_CircTest(
    CircRNACount[, c(1 : 3, arg_condition_columns)], # we always need the first 3 columns
    LinearCount[, c(1 : 3, arg_condition_columns)], # we always need the first 3 columns
    CircCoordinates,
    rep(arg_groups, arg_replictes),
    rep(arg_condition_list, arg_replictes),
    arg_output_label,
    arg_output_name,
    arg_filter_sample,
    arg_filter_count,
    arg_max_fdr,
    arg_max_plots,
    arg_replictes,
    arg_percent_filter
)
