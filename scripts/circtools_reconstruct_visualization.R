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

### imports
message("Loading required packages")
suppressMessages(library(ggplot2))
suppressMessages(require(data.table))
suppressMessages(require(ggsignif))
message("Done loading packages")

# may not be optimal, but we do not want warning for now
options(warn = - 1)

create_quantile_plot <- function(table) {

  libs <- unique(table$group)
  maxval <- max(table$length, na.rm=T)

  for (lib in libs) {
    subset <- table[table$group==lib,]
    dt <- data.table(x=subset$group,y=subset$length)
    dens <- density(dt$y)
    df <- data.frame(x=dens$x, y=dens$y)
    probs <- c(0.90, 0.95, 0.99)
    quantiles <- quantile(dt$y, prob=probs)
    probs <- c(probs,paste("> ",probs[length(probs-1)]))
    df$quant <- factor(findInterval(df$x,quantiles))
    plot <- ggplot(df, aes(x,y)) + geom_line(colour="black", size=0.1) + geom_ribbon(aes(ymin=0, ymax=y, fill=quant))
    if (arg_colour_mode == "bw" ) {
        plot <- plot + scale_fill_grey(start = 0, end = .9, guide="legend", labels=probs, name="Quantiles")
    } else {
        plot <- plot + scale_fill_hue(guide="legend", labels=probs, name="Quantiles")
    }
    plot <- plot +
    labs(title = "Circular RNA reconstruction results",
    subtitle = paste("CircRNA length - quantile plot for sample \"", lib, "\"", sep="")) +
    labs(y = "Density") +
    labs(x = "CircRNA length") +
    xlim(0,maxval) +
    theme(legend.position = c(.95, .95),
    legend.justification = c("right", "top"))
    print(plot)
  }
}

args <- commandArgs(trailingOnly = TRUE)

arg_reconstruct_data <- args[1] # path is string
arg_groups <-   unlist(lapply(strsplit(args[2],","), as.numeric)) # list of ints
arg_condition_list  <- strsplit(args[3],",")[[1]] # list of strings
arg_colour_mode <- args[4]

geom_signif_list <- (combn(arg_condition_list,2, simplify=F))

# check colour mode
if (arg_colour_mode != "colour" & arg_colour_mode != "bw" ) {
    print(arg_colour_mode)
    message("Please specify the colour model: colour or bw")
    quit()
}

group_length <- length(arg_condition_list)
message(paste(group_length, " different groups provided", sep=""))
names <- unlist(lapply(seq(1, length(arg_groups)), function(x) {return(arg_condition_list[arg_groups[x]])}))

exon_counts <- system(paste("awk \'{split($11,array,\",\"); sum = 0; for( i = 1; i <= length(array); i++ ) {sum += array[i]}; print FILENAME\"\t\"$3-$2\"\t\"$10\"\t\"sum}\' *.exon_counts.bed | sed \'s/\\.exon_counts\\.bed//g\' | egrep -v \\W0\\W",sep=""), intern = TRUE)
exon_counts <- read.delim(textConnection(exon_counts), header = F, sep = "\t")


mate_status <- system(paste("awk \'{print FILENAME\"\\t\"$6\"\\t\"$7\"\\t\"$8}\' *.mate_status.txt | sed \'s/\\.mate_status\\.txt//g\' | grep -v single  | awk \'{print $0\"\\t\"$3/($2+$3+$4)}\'",sep=""), intern = TRUE)
mate_status <- read.delim(textConnection(mate_status), header = F, sep = "\t")

isoforms <- system(paste("grep -H -o -n \',\' *alternative_splicing.txt | cut -d : -f 1,2 |  uniq -c |  sed \'s/\\.alternative.*//g\' | awk \'{print $2\"\\t\"$1}\' ",sep=""), intern = TRUE)
isoforms <- read.delim(textConnection(isoforms), header = F, sep = "\t")



circ_length <- system(paste("awk \'{if($2>0){print FILENAME\"\\t\"$2}}\' *.coverage_profiles/cluster_association.all_circles.tsv | sed \'s/\\.coverage_profiles\\/cluster_association\\.all_circles\\.tsv//g\' | grep -v length",sep=""), intern = TRUE)
circ_length <- read.delim(textConnection(circ_length), header = F, sep = "\t")


colnames(exon_counts) <- c("Library2", "Length_total", "Exons", "Length_exons")
colnames(isoforms) <- c("Library", "num")
colnames(circ_length) <- c("Library", "length")
colnames(mate_status) <- c("Library", "Single", "Double", "Unknown", "Ratio")

sequencing_lib_name <- levels(unique(isoforms$Library))
detailed_names <- names

data <- data.frame((unique(isoforms$Library)))
colnames(data) <- c("Library")

for (lib in unique(isoforms$Library)) {
  data$count[data$Library==lib] <- nrow(circ_length[circ_length$Library==lib,])
}


for (i in 1 : length(sequencing_lib_name)) {
    exon_counts$group[exon_counts$Library == sequencing_lib_name[i]] <- detailed_names[i]
    mate_status$group[mate_status$Library == sequencing_lib_name[i]] <- detailed_names[i]
    isoforms$group[isoforms$Library == sequencing_lib_name[i]] <- detailed_names[i]
    circ_length$group[circ_length$Library == sequencing_lib_name[i]] <- detailed_names[i]
    data$group[data$Library == sequencing_lib_name[i]] <- detailed_names[i]
}

############### bare fuchs data, circRNA length
pdf("results.pdf", width= 8, height=5.5 , title="circtools reconstruction results")

# # total raw number of FUCHS circRNAs
plot <- ggplot(data = data, aes(x = group, y = count, fill = group)) +
    geom_boxplot(colour="black", size=0.1)

    if (arg_colour_mode == "bw" ) {
        plot <- plot + scale_fill_grey(start = 0, end = .9)
    }
    plot <- plot + geom_signif(comparisons = geom_signif_list, step_increase = 0.1, map_signif_level = T) +
    labs(title = "Circular RNA reconstruction results",
    subtitle = "Absolute count - normal scale") +
    labs(y = "Total number of circular RNAs") +
    labs(x = "Sample") +
    labs(fill = "Sample") +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
print(plot)

# # total raw number of FUCHS circRNAs
plot <- ggplot(data = data, aes(x = group, y = count, fill = group)) +
  geom_boxplot(colour="black", size=0.1)

    if (arg_colour_mode == "bw" ) {
        plot <- plot + scale_fill_grey(start = 0, end = .9)
    }
    plot <- plot + geom_signif(comparisons = geom_signif_list, step_increase = 0.1, map_signif_level = T) +
    labs(title = "Circular RNA reconstruction results",
		subtitle = "Absolute count - log10 scale") +
    labs(y = "Total number of circular RNAs (log10)") +
    labs(x = "Sample") +
    labs(fill = "Sample") +
    scale_y_log10() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
print(plot)

# Quantile length plots

create_quantile_plot(circ_length)

# more stuff: single/double breakpoint fragments, exon length
plot <- ggplot(data = isoforms, aes(x = group, y = num, fill = group)) +
    geom_boxplot() + theme_classic()

    if (arg_colour_mode == "bw" ) {
        plot <- plot + scale_fill_grey(start = 0, end = .9)
    }
    plot <- plot + geom_signif(comparisons = geom_signif_list, step_increase = 0.1, map_signif_level = T) +
    labs(title = "Circular RNA reconstruction results",
    subtitle = "# circRNA isoforms per host gene") +
    labs(y = "Number of isoforms") +
    labs(x = "Sample") +
    labs(fill = "Sample") +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
print(plot)

plot <- ggplot(data = exon_counts, aes(x = group, y = Length_total, fill = group)) +
    geom_boxplot() + theme_classic()

    if (arg_colour_mode == "bw" ) {
        plot <- plot + scale_fill_grey(start = 0, end = .9)
    }
    plot <- plot + geom_signif(comparisons = geom_signif_list, step_increase = 0.1, map_signif_level = T) +
    scale_y_log10() +
    labs(title = "Circular RNA reconstruction results",
    subtitle = "Total length of circRNAs") +
    labs(y = "Total length (exons + introns)") +
    labs(x = "Sample") +
    labs(fill = "Sample") +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
print(plot)

plot <- ggplot(data = exon_counts, aes(x = group, y = Length_exons, fill = group)) +
    geom_boxplot() + theme_classic()

    if (arg_colour_mode == "bw" ) {
        plot <- plot + scale_fill_grey(start = 0, end = .9)
    }
    plot <- plot + geom_signif(comparisons = geom_signif_list, step_increase = 0.1, map_signif_level = T) +
    scale_y_log10() +
    labs(title = "Circular RNA reconstruction results",
    subtitle = "Exon-based length of circRNAs") +
    labs(y = "Length (exons only)") +
    labs(x = "Sample") +
    labs(fill = "Sample") +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
print(plot)

plot <- ggplot(data = exon_counts, aes(x = group, y = Exons, fill = group)) +
    geom_boxplot() + theme_classic()
    if (arg_colour_mode == "bw" ) {
        plot <- plot + scale_fill_grey(start = 0, end = .9)
    }
    plot <- plot + geom_signif(comparisons = geom_signif_list, step_increase = 0.1, map_signif_level = T) +
    labs(title = "Circular RNA reconstruction results",
    subtitle = "Number of exons per circRNA") +
    labs(y = "# exons per circRNA") +
    labs(x = "Sample") +
    labs(fill = "Sample") +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
print(plot)

plot <- ggplot(data = mate_status, aes(x = group, y = Single, fill = group)) +
    geom_boxplot() + theme_classic()

    if (arg_colour_mode == "bw" ) {
        plot <- plot + scale_fill_grey(start = 0, end = .9)
    }
    plot <- plot + geom_signif(comparisons = geom_signif_list, step_increase = 0.1, map_signif_level = T) +
    scale_y_log10() +
    labs(title = "Circular RNA reconstruction results",
    subtitle = "Single breakpoint circRNAs (absolute)") +
    labs(y = "# single breakpoints") +
    labs(x = "Sample") +
    labs(fill = "Sample") +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
print(plot)

plot <- ggplot(data = mate_status, aes(x = group, y = Double, fill = group)) +
    geom_boxplot() + theme_classic()
    if (arg_colour_mode == "bw" ) {
        plot <- plot + scale_fill_grey(start = 0, end = .9)
    }
    plot <- plot + geom_signif(comparisons = geom_signif_list, step_increase = 0.1, map_signif_level = T) +
    scale_y_log10() +
		labs(title = "Circular RNA reconstruction results",
		subtitle = "Double breakpoint circRNAs (absolute)") +
    labs(y = "# double breakpoints") +
    labs(y = "# exons per circRNA") +
    labs(x = "Sample") +
    labs(fill = "Sample") +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
print(plot)

# plot <- ggplot(data = mate_status, aes(x = group, y = Unknown, fill = group)) +
#     geom_boxplot() + theme_classic()
#
#     if (arg_colour_mode == "bw" ) {
#         plot <- plot + scale_fill_grey(start = 0, end = .9)
#     }
#     plot <- plot +
#     scale_y_log10() +
# 		labs(title = "Circular RNA reconstruction results",
# 		subtitle = "Unknown breakpoint circRNAs (absolute)") +
#     labs(y = "# unknown breakpoints") +
#     labs(y = "# exons per circRNA") +
#     labs(x = "Sample") +
#     labs(fill = "Sample") +
#     theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
# print(plot)

plot <- ggplot(mate_status, aes(x = group, y = Ratio, fill = group)) +
		geom_boxplot() + theme_classic()
    if (arg_colour_mode == "bw" ) {
        plot <- plot + scale_fill_grey(start = 0, end = .9)
    }
    plot <- plot + geom_signif(comparisons = geom_signif_list, step_increase = 0.1, map_signif_level = T) +
		labs(title = "Circular RNA reconstruction results",
		subtitle = "Ratio of double breakpoints") +
    labs(y = "Ratio double breakpoints") +
    labs(x = "Sample") +
    labs(fill = "Sample") +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
print(plot)

invisible(dev.off())
