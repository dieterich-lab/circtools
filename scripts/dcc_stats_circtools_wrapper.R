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

getStats <- function(STARfolder)
{
    # Read STAR log file
    logFile <- paste(STARfolder, "Log.final.out", sep = "/");

    #check if file exists
    log <- read.delim(logFile, as.is = T, sep = "|");
    log <- as.numeric(gsub("\t", "", log[7, 2]))
    return(log);
}

# may not be optimal, but we do not want warning for now
options(warn = - 1)

# quiet mode
options(echo = FALSE)

# read arguments
args <- commandArgs(trailingOnly = TRUE)

arg_dcc_data <- args[1] # path is string
arg_replictes <- as.integer(args[2]) # integer
arg_condition_list <- strsplit(args[3], ",")[[1]] # list of strings
arg_condition_columns <- lapply(strsplit(args[4], ","), as.numeric) # list of integers
arg_condition_columns <- unlist(arg_condition_columns)
arg_groups <- unlist(lapply(strsplit(args[5], ","), as.numeric)) # list of strings
arg_output_directory <- args[6] # string
arg_star_folder <- args[7] # path is string

q()

## load complete data set
message("Loading CircRNACount")
CircRNACount <- read.delim(paste(arg_dcc_data, "CircRNACount", sep="/"), header = T)

message("Loading LinearRNACount")
LinearCount <- read.delim(paste(arg_dcc_data, "LinearCount", sep="/"), header = T)


circ_counts_summed <- apply(CircRNACount[, - c(1 : 3)], 2, sum)
linear_counts_summed <- apply(LinearCount[, - c(1 : 3)], 2, sum)


todo <- scan("/Volumes/new_home/manuscripts//Define_Circles_Step1/workflow/DCC_new/star_folder.txt", character())
todoNam <- scan("/Volumes/new_home/manuscripts//Define_Circles_Step1/workflow/DCC_new/star_name.txt", character())

UniqueMappers <- numeric();

for (k in todo)
{
    UniqueMappers <- c(UniqueMappers, getStats(k));
}

names(UniqueMappers) <- names(circ_counts_summed)

#tissue.group<-factor(rep(c("liver","heart","cerebellum","hippocampus"),6));
#age.group<-factor(c(rep("young",12),rep("old",12)));

pdf(paste(arg_output_directory, "Step1_circSeq_junction_read_counts.pdf", sep = "/"));
plot(x = circ_counts_summed, y = linear_counts_summed, log = "xy", xlab = "CircReadCount", ylab = "LinearReadCount", col = rep(c("blue", "red"), 9))
text(x = circ_counts_summed, y = linear_counts_summed, names(circ_counts_summed), col = "black", cex = 0.5, pos = 1, srt = 45)
#plot(x=eins,y=zwei,log="xy",xlab="CircReadCount",ylab="LinRead",col=tissue.group,pch=as.numeric(age.group))
#legend("bottomleft", c("liver","heart","cerebellum","hippocampus"),fill=tissue.group)
legend("bottomright", c("RNASe+", "RNASe-"), fill = c("red", "blue"))
dev.off();



NoOfCircles <- apply(CircRNACount[, - c(1 : 3)], 2, function(x){length(which(x > 1))})
names(UniqueMappers) <- names(NoOfCircles)

pdf(paste(arg_output_directory, "Step1_circles_per_unique_mapped_reads.pdf", sep = "/"));
plot(x = UniqueMappers, y = NoOfCircles, log = "xy", ylab = "Number of detected circles", xlab = "Number of mapped reads", col = rep(c("blue", "red"), 9));
text(x = UniqueMappers, y = NoOfCircles, names(UniqueMappers), cex = 0.5, col = "black", srt = 45);
legend("topright", c("RNASe+", "RNASe-"), fill = c("red", "blue"))

ref = order(NoOfCircles / (UniqueMappers / 1000000))
barplot(sort(NoOfCircles / (UniqueMappers / 1000000)), las = 2, col = rep(c("blue", "red"), 9)[ref], main = "Circles_per_million_unique_mapper")
legend("topleft", c("RNASe+", "RNASe-"), fill = c("red", "blue"))

dev.off();
