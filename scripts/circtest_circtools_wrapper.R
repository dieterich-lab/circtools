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
library(gridExtra)

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
arg_condition_columns <- as.integer(strsplit(args[4],",")[[1]]) # list of integers
arg_output_name <- args[5] # string
arg_max_fdr <- as.numeric(args[6]) # float
arg_output_size <- args[7] # string
arg_filter_sample <- as.integer(args[8]) # integer
arg_filter_count <- as.integer(args[9]) # integer


run_CircTest = function(CircRNACount, LinearCount, CircCoordinates, groups, indicators, label, xlsName, pdfName) {

    message("Filtering circRNA counts")
    CircRNACount_filtered <- Circ.filter(circ = CircRNACount, linear = LinearCount,
    Nreplicates = 3, filter.sample = 3, filter.count = 5, percentage = 0.01)

    message("Filtering circRNA coordinates")
    CircCoordinates_filtered <- CircCoordinates[rownames(CircRNACount_filtered),]

    message("Filtering linear RNA counts")
    LinearCount_filtered <- LinearCount[rownames(CircRNACount_filtered),]

    message("running circTest")

    data = Circ.test(CircRNACount_filtered, LinearCount_filtered, CircCoordinates_filtered, group = groups)

    message("plotting")

    pdf(file = pdfName, width = 8.27, height = 11.69)  # DIN A4 in inches

    if (nrow(data$summary_table) > 100) {
        max <- 100
    } else {
        max <- nrow(data$summary_table)
    }

    # print(head(data$summary_table))

    for (i in rownames(data$summary_table[1 : max,])) {

        # p[[i]] <- Circ.ratioplot( CircRNACount_filtered, LinearCount_filtered,
        # CircCoordinates_filtered, plotrow=i, size=16, gene_column=4,
        print(i)

        Circ.ratioplot(CircRNACount_filtered, LinearCount_filtered, CircCoordinates_filtered,
        plotrow = i, size = 16, gene_column = 4, groupindicator1 = indicators,
        x = "", y = "", lab_legend = label)
    }
    # print(p) do.call(grid.arrange,(c(p,ncol=2, nrow=4))) dev.off();

    dev.off()

    message("creating Excel sheets")

    ## Significant result show in a summary table


    write.table(na.omit(data$summary_table[1 : MAX_LINES,]),
                file = paste(xlsName, ".csv", sep = "") ,
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

    addWorksheet(wb, sheetName = "CircleSkips")
    writeDataTable(wb, sheet = 4, x = CircSkip[idx,])

    saveWorkbook(wb, paste(xlsName, ".xlsx", sep = "") , overwrite = TRUE)
}

###########################################################

## load complete data set
message("Loading CircRNACount")
CircRNACount <- read.delim(paste(arg_dcc_data, "CircRNACount", sep="/"), header = T)

message("Loading LinearRNACount")
LinearCount <- read.delim(paste(arg_dcc_data, "LinearCount", sep="/"), header = T)

message("Loading CircCoordinates")
CircCoordinates <- read.delim(paste(arg_dcc_data, "CircCoordinates", sep="/"), header = T)

print(arg_dcc_data)
print(arg_replictes)
print(arg_condition_list)
print(arg_condition_columns)
print(arg_output_name)
print(arg_max_fdr)
print(arg_output_size)
print(arg_filter_sample)
print(arg_filter_count)


# run_CircTest(
# CircRNACount[, c(1 : 3, 10 : 15)],
# LinearCount[, c(1 : 3, 10 : 15)],
# CircCoordinates,
# c(rep(c(2, 1), 3)),
# c(rep(c("cond1", "cond2"), 3)),
# "Run",
# "./circTest/out",
# "./circTest/out.pdf"
# )

warnings()
