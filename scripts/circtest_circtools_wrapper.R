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

###########################################################
runCirc = function(CircRNACount, LinearCount, CircCoordinates, groups, indicators,
    label, xlsName, pdfName) {

    message("Filtering circRNA counts")
    CircRNACount_filtered <- Circ.filter(circ = CircRNACount, linear = LinearCount,
        Nreplicates = 3, filter.sample = 3, filter.count = 5, percentage = 0.01)

    message("Filtering circRNA coordinates")
    CircCoordinates_filtered <- CircCoordinates[rownames(CircRNACount_filtered),
        ]

    message("Filtering linear RNA counts")
    LinearCount_filtered <- LinearCount[rownames(CircRNACount_filtered), ]

    message("running circTest")

    data = Circ.test(CircRNACount_filtered, LinearCount_filtered, CircCoordinates_filtered,
        group = groups)

    message("plotting")

    pdf(file = pdfName, width = 8.27, height = 11.69)  # DIN A4 in inches

    if (nrow(data$summary_table) > 100) {
        max <- 100
    } else {
        max <- nrow(data$summary_table)
    }

    # print(head(data$summary_table))

    for (i in rownames(data$summary_table[1:max, ])) {

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

    ## fix 'not find zip' error from openxlsx
    Sys.setenv(R_ZIPCMD = "/usr/bin/zip")
    require(openxlsx)

		write.table(na.omit(data$summary_table[1:10000, ]), file = paste(xlsName,".csv",sep="") , quote = FALSE, sep = "\t",
            eol = "\n", row.names = FALSE, col.names = FALSE)

    wb <- createWorkbook()
    addWorksheet(wb, sheetName = "Significant circles")
    writeDataTable(wb, sheet = 1, x = data$summary_table[1:10000, ])

    idx <- rownames(data$summary_table[1:10000, ])

    addWorksheet(wb, sheetName = "Circle Counts")
    writeDataTable(wb, sheet = 2, x = CircRNACount[idx, ])

    addWorksheet(wb, sheetName = "Linear Counts")
    writeDataTable(wb, sheet = 3, x = LinearCount[idx, ])

    addWorksheet(wb, sheetName = "CircleSkips")
    writeDataTable(wb, sheet = 4, x = CircSkip[idx, ])

    saveWorkbook(wb, paste(xlsName,".xlsx",sep="") , overwrite = TRUE)

}

###########################################################

## load complete data set
message("Loading CircRNACount")
CircRNACount <- read.delim("CircRNACount", header = T)

message("Loading LinearRNACount")
LinearCount <- read.delim("LinearCount", header = T)

message("Loading CircCoordinates")
CircCoordinates <- read.delim("CircCoordinates", header = T)

message("Loading CircSkipJunctions")
CircSkip <- read.delim("CircSkipJunctions", header = T)

runCirc(
	CircRNACount[, c(1:3,10:15)],
	LinearCount[, c(1:3,10:15)],
	CircCoordinates,
	c(rep(c(2,1),3)),
	c(rep(c("cond1","cond2"),3)),
	"Run",
	"./circTest/out",
	"./circTest/out.pdf"
)


warnings()
