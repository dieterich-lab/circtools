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

suppressMessages(library(ballgown))
suppressMessages(library(edgeR))
suppressMessages(library(ggbio))
suppressMessages(library(ggfortify))
suppressMessages(library(openxlsx))
suppressMessages(library(GenomicRanges))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(biomaRt))

message("Done loading packages")

# may not be optimal, but we do not want warning for now
options(warn=-1)

# fix 'could not find zip' error from openxlsx
Sys.setenv(R_ZIPCMD = "/usr/bin/zip")
library(openxlsx)


# let's now get the command line arguments

# quiet mode
suppressMessages(options(echo=FALSE))

# read arguments
args <- commandArgs(trailingOnly = TRUE)

# assign input to script variables

arg_dcc_data <- args[1] # path is string
arg_replicates <- unlist(lapply(strsplit(args[2],","), as.numeric)) # list of integers
arg_condition_list <- strsplit(args[3],",")[[1]] # list of strings
arg_condition_list <- unlist(arg_condition_list)
arg_condition_columns <- lapply(strsplit(args[4],","), as.numeric) # list of integers
arg_condition_columns <- unlist(arg_condition_columns)
arg_groups <- lapply(strsplit(args[5],","), as.numeric) # list of integers
arg_groups <- unlist(arg_groups)
arg_output_directory <- args[6] # string
arg_ballgown_directory <- args[7] # string
arg_gtf_file <- args[8] # string
arg_circTest_file <- args[9] # string
arg_head_header <- as.logical(args[10])
arg_species <- args[11]


# load a pkg from a string
arg_ens_db <- switch(
    arg_species,
    "hs" = "hsapiens_gene_ensembl",
    "rn" = "rnorvegicus_gene_ensembl",
    "mm" = "mmuscules_gene_ensembl",
    "ss" = "sscrofa_gene_ensembl"
)

## load complete data set
message("Loading CircRNACount")
CircRNACount <- read.delim(paste(arg_dcc_data, "CircRNACount", sep="/"), header = T)

message("Loading CircCoordinates")
CircCoordinates <- read.delim(paste(arg_dcc_data, "CircCoordinates", sep="/"), header = T)

# read sub directories containing the ballgown runs and return list
ballgownRuns <- as.list(list.files(arg_ballgown_directory, full.names = TRUE))

# output directory
baseDir <- arg_output_directory

# group mapping
group <- unlist(lapply(arg_groups, function(x) {return(arg_condition_list[x])}))

# sample<>replicate mapping
id <- unlist(lapply(seq(1, length(arg_replicates)), function(x) {return((paste(arg_condition_list[x], arg_replicates[x], sep="_R")))}))

bg.dccDF <- data.frame( id=id, group=group)

bg_dirs_to_work <- unlist(lapply(arg_condition_columns, function(x) {return(ballgownRuns[x-3])}))

message("Starting ballgown processing")

bg <- ballgown(bg_dirs_to_work, verbose=TRUE)

message("Preparing necessary data structures")

whole_exon_table <- eexpr(bg, 'all')# eexpr -> exon level
whole_intron_table <- iexpr(bg, 'all')# iexpr -> intron level

t2g<-indexes(bg)$t2g # transcript / gene table
e2t<-indexes(bg)$e2t # exon / transcript table

e2g<-unique(merge(e2t, t2g, by.x=2, by.y=1)[, 2:3])# exon / gene table merging
idx<-names(which(table(e2g[, 1])==1))# idx to sort out non-gene exons
e2g<-subset(e2g, e2g[, 1]%in% idx)# exon / gene transcript / merging

# redo everything based on counts
# use mrcount -> multi-map-corrected number of reads overlapping the exon/intron
nonMMexonCount<-whole_exon_table[, c(1, grep("mrcount", colnames(whole_exon_table)))]

# only exons > 10 reads throught all mappings
idx<-which(apply(nonMMexonCount[, -1], 1, sum)>10)
nonMMexonCount<-nonMMexonCount[idx,]

# merge exon couns and ( exon / gene / transcript) table
e2g.counts<-merge(nonMMexonCount, e2g, by.x=1, by.y=1)

# extract first and last columns
geneBaseTable<-e2g.counts[, c(1, ncol(e2g.counts))]

# combine tables, add chr, strand, start and stop
# we now have a large table with exon to gene mappings and base annotation
geneBaseTable<-merge(geneBaseTable, whole_exon_table[, 1:5], by.x=1, by.y=1)

# add exon number to gene name. e.g.: ENSG00000223972 -> ENSG00000223972.8
geneBaseTable[, 1]<-gsub(" ", "", apply(geneBaseTable[, c(2, 1)],
1, paste, sep="", collapse="."))

#print(head(geneBaseTable))
# bundle together gene names, gene/exon combos in a data frame
geneIDs <- data.frame(GeneID=geneBaseTable[, "g_id"],
Gene.Exon=geneBaseTable[, "e_id"])

message("Setting treatment and conditions")

treatment <- group
condition <- id
batch <- unlist(lapply(arg_replicates,function(x){paste("R",x,sep="")}))

# remove first and last column
e2g.minimal=e2g.counts[,-c(1,ncol(e2g.counts))]

# indices of exons in e2g.minimal table with > 40 counts throughout all samples
idx<-which(apply(e2g.minimal,1,sum)>4*10)

# select > 40 count exons from the slimmed exon / gene table
e2g.minimal=e2g.minimal[idx,]

# set appropriate gene names
rownames(e2g.minimal)<-geneIDs[idx,2]

# only keep gene + exon names that also have > 40 count
geneIDs=geneIDs[idx,]

# extract multi-exon genes
genesWithMultipleExons=names(which(table(as.character(geneIDs[,1]))>1))

# extract single exon genes
genesWithSingleExon=names(which(table(as.character(geneIDs[,1]))==1))

message(paste("Found ",length(genesWithMultipleExons), " multi exon genes", sep = ""))

message(paste("Found ", length(genesWithSingleExon), " single exon genes", sep = ""))

# subset the geneID list to only contain multi-exon genes
geneIDs=subset(geneIDs,geneIDs[,1] %in% genesWithMultipleExons)

# intersect the minimal e2g list with with the geneIDs to get exon (gene)
# -> count relation
e2g.minimal=e2g.minimal[intersect(rownames(e2g.minimal),geneIDs[,2]),]

# DE analysis with edgeR starts here
z <- DGEList(counts=e2g.minimal, genes=geneIDs)
z <- calcNormFactors(z)

# build the design for edgeR
design <- model.matrix(~ treatment)

message("Starting dispersion estimation")
z <- estimateGLMCommonDisp(z,design)
z <- estimateGLMTrendedDisp(z,design)
z <- estimateGLMTagwiseDisp(z,design)

message("Fitting model...")
fit <- glmFit(z,design)

message("Testing for differential exon usage")
# Test For Differential Exon Usage
# Given a negative binomial generalized log-linear model fit at the exon level,
# test for differential exon usage between experimental conditions.
diffSplicedGenes<-diffSpliceDGE(  fit,
                                geneid="GeneID",
                                exonid="Gene.Exon",
                                prior.count=0.125,
                                verbose=TRUE
                              )

# Top Table Of Differentially Spliced Genes Or Exons
# Top table ranking the most differentially spliced genes or exons
topSplicedGenes=topSpliceDGE( diffSplicedGenes,
                            FDR=0.01,
                            test="Simes",
                            number=nrow(diffSplicedGenes$gene.geneIDs)
                          )

# merge information from gene2 table with the DS exons from above
splicedExonDF <- merge(
                      geneBaseTable,
                      data.frame(
                                    ID = names(diffSplicedGenes$exon.p.value),
                                    Pval=diffSplicedGenes$exon.p.value,
                                    log2FC=diffSplicedGenes$coefficients
                                  ),
                      by.x=1, by.y=1
                    )

#####################################################

# build new table with selected columns
splicedExonDFfixed=splicedExonDF[,c("chr","start","end","log2FC","Pval")];
splicedExonDFfixed[,2]<-splicedExonDFfixed[,2]-1 # fix start pos (0-based now)
splicedExonDFfixed[,3]<-splicedExonDFfixed[,3] # fix stop pos (added +1 as it was too short)
splicedExonDFfixed[,5]<-log2(splicedExonDFfixed[,5]) # get log2 for Pval column

message("Writing bed files...")
############# Writing BED files ##############

# write BED track header for FC track
write(  paste("track type=bedGraph name=\"edgeR_exon_foldChange\"",
       "description=\"Exon fold change over all total RNA data sets\"",
       "visibility=full color=200,100,0 altColor=0,100,200 priority=20",
       sep=""),file=paste(baseDir,"exon_fc_track.bedgraph", sep=""));

# write actual content (FCs)
write.table(
          splicedExonDFfixed[,1:4],
          file=paste(baseDir,"exon_fc_track.bedgraph",sep=""),
          quote=F,
          sep="\t",
          append=T,
          row.names=F,
          col.names=F
        )

############

# write BED track header for pValue track
write(paste("track type=bedGraph name=\"edgeR_exon_Pvalue\"",
         "description=\"Pvalue over all data sets\"",
         "visibility=full color=200,100,0 altColor=0,100,200 priority=20",
         sep=""), file=paste(baseDir,"exon_pval_track.bedgraph", sep=""));

# write actual content (PValues)
write.table(  splicedExonDFfixed[,c(1:3,5)],
            file=paste(baseDir,"exon_pval_track.bedgraph",sep=""),
            quote=F,
            sep="\t",
            append=T,
            row.names=F,
            col.names=F
          )

############# Done writing BED files ##############

#Read back-splice enrichment

# head(splicedExonDFfixed[,5])
#colnames(dccDF)<-c("chr","start","end","strand")
splicedExonDFfixed <- subset(
                            splicedExonDFfixed,
                            splicedExonDFfixed[,5]<=-5 & splicedExonDFfixed[,4]>0
                          )

# create GRanges object from diff. spliced genes table and assign PValue and FC
multiExonRanges <- makeGRangesFromDataFrame(splicedExonDFfixed[,1:3]);
mcols(multiExonRanges)$log2FC <- splicedExonDFfixed[,4]
mcols(multiExonRanges)$Pval <- splicedExonDFfixed[,5]

# Extract genes with single exons based on the table created before
singleExonDF=subset(geneBaseTable[,c("chr","start","end","strand")],
                  geneBaseTable[,2] %in% genesWithSingleExon
                  )

# create appropriate GRanges objects
singleExonRanges <- makeGRangesFromDataFrame(singleExonDF);

ensembl <- useMart("ensembl", dataset = arg_ens_db, host = "ensembl.org")

genemap <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "description"), values = as.character(topSplicedGenes$GeneID), mart = ensembl)

colnames(genemap) <-  c("SYMBOL", "GENEID", "DESCRIPTION")

topSplicedGenesMartData <- merge(
  topSplicedGenes,
  genemap,
  by.x = "GeneID",
  by.y = "GENEID",
  all.x = TRUE
)

# correct again for to 0-based positions
dccDF<-CircCoordinates
dccDF[,2]<-dccDF[,2]-1
dccDF[,3]<-dccDF[,3] #was too short had to extend again

# write BED file
colnames(dccDF)<-c("chr","start","end","GeneName","JType","strand")

message("Writing DCC prediction BED file")

#Print bed file, read counts
write(paste("track name=\"DCC_predictions\" description=\"Predictions ",
          "over all data sets\" itemRgb=\"on\"", sep=""),
           file=paste(baseDir,"dcc_predictions_track.bed",sep=""));

# create new data frame for improved visualisation
# only take useful columns from DCC DF

new_df <- data.frame(dccDF[,c(1,2,3,4)])
new_df$score <- 0
new_df$strand <- dccDF[,6]
new_df$tstart <- dccDF[,2]+as.integer(((dccDF[,3]-dccDF[,2])*0.10))
new_df$tstop <- dccDF[,3]-as.integer(((dccDF[,3]-dccDF[,2])*0.10))
new_df$rgb <- paste(0,sep=",")
new_df$blocks <- 2
new_df$bsize <- paste(as.integer(((dccDF[,3]-dccDF[,2])*0.05)),as.integer(((dccDF[,3]-dccDF[,2])*0.05)),",", sep=",")
new_df$bstart <- paste(as.integer(dccDF[,2]+as.integer(((dccDF[,3]-dccDF[,2])*0.10))),(dccDF[,3]-as.integer(((dccDF[,3]-dccDF[,2])*0.10))-as.integer(((dccDF[,3]-dccDF[,2])*0.05))) ,sep=",")

# new_df = as.matrix(new_df)



write.table(new_df,
            file=paste(baseDir,"dcc_predictions_track.bed",sep=""),
            quote=F,
            sep="\t",
            append=T,
            row.names=F,
            col.names=F
          )

#genomicRanges <- makeGRangesFromDataFrame(dccDF);
#do we need to test if enrichment happens inside of BSJ ?!?

message("Reading and integrating CircTest results")

# pull together circle predictions from DCC
CircPred <- data.frame(
                      Gene=unique(CircCoordinates[,"Gene"]),
                      DCC_predicted=rep(1,length(unique(CircCoordinates[,"Gene"])))
                    )

# enrich top spliced Genes with DCC Circle prediction
topSplicedGenesMartDCC=merge(
                            topSplicedGenesMartData,
                            CircPred,
                            by.x=5,
                            by.y=1,
                            all.x=T
                          )

circTestSummary<- read.delim(
                        arg_circTest_file,
                        header=arg_head_header,
                        as.is=T
                      )[,c(1:9)]

cols.ct <- c("Chr", "Start", "End", "Gene", "JunctionType", "Strand", "Start.End.Region", "OverallRegion", "sig_p")
colnames(circTestSummary) <- cols.ct

#write bed file
# assign circTest summary results to DF
dccDF<-circTestSummary[,c(1,2,3,6)]

# reorder columns
colnames(dccDF)<-c("chr","start","end","strand")

# create a GRanges object from DCC data frame
genomicRanges <- makeGRangesFromDataFrame(dccDF)

# compute overlap with exon enrichment for multi and single exon genes
multiExonOverlap <- as.matrix(findOverlaps(genomicRanges,multiExonRanges))
singleExonOverlap <- as.matrix(findOverlaps(genomicRanges,singleExonRanges))

# overwrite DCC data frame with more columns
dccDF<-circTestSummary[,1:7]
dccDF[,2]<-dccDF[,2]-1
# reorder columns
colnames(dccDF)<-c("chr","start","end","GeneName","JType","strand")

message("Writing back splice junction enriched BED file")

# Print BED file, read counts
write(paste("track name=\"DCC_BSJ_enriched\"",
          "description=\"Back-Splice junction enriched over all data ",
          "sets 1% FDR\" color=\"green\"", sep=""),
           file=paste(baseDir,"dcc_bsj_enriched_track.bed",sep=""));

write.table(  dccDF,
            file=paste(baseDir,"dcc_bsj_enriched_track.bed",sep=""),
            quote=F,
            sep="\t",
            append=T,
            row.names=F,
            col.names=F
          )

# merge circTest results @ unique exon overlaps with topSplicedGenes from
# ballgown; also remove second column
RNAse_RenrichedCircTest <- merge(
                          circTestSummary[unique(multiExonOverlap[,1]),],
                          topSplicedGenesMartData,by.x=4,by.y=5
                        )

# sort by PValue
RNAse_RenrichedCircTest <- RNAse_RenrichedCircTest[order(RNAse_RenrichedCircTest[,"sig_p"]),]


colnames(RNAse_RenrichedCircTest) <- c( "Gene",
                                        "Chr",
                                        "Start",
                                        "End",
                                        "JunctionType",
                                        "Strand",
                                        "Start.End.Region",
                                        "OverallRegion",
                                        "sig_p",
                                        "GeneID",
                                        "NExons",
                                        "P.Value",
                                        "FDR",
                                        "description"
                                        )


# create a data frame with all gene names from circtest and mark them
CircbackSpliceEnrich <- data.frame(
                                  Gene=unique(circTestSummary[,"Gene"]),
                                  RNaseR_enriched=rep(1,length(unique(circTestSummary[,"Gene"])))
                                )

# sort the top spliced genes from DCC with mart annotation by FDR
topSplicedGenesMartDCC <- topSplicedGenesMartDCC[order(topSplicedGenesMartDCC[,"FDR"]),]

#print out single exon subset in circTest summary
# circTestSummary[unique(singleExonOverlap[,1]),]


# combine DCC top spliced genes with the BS enriched data frame
mainTable <- merge( topSplicedGenesMartDCC,
                  CircbackSpliceEnrich,
                  by.x=1,
                  by.y=1,
                  all.x=T
                  )
message("Writing Excel file")

# create an excel workbook
wb <- createWorkbook()

# exons with 1% FDR from the main table
addWorksheet(wb, sheetName = "Exon FDR 1% (ballgown)")
writeDataTable(wb, sheet = 1, x=mainTable[order(mainTable[,"FDR"]),])

# Enriched Back-Splice junctions / FDR 1%
addWorksheet(wb, sheetName = "enriched BSJ FDR 1% (CircTest)")

print(head((RNAse_RenrichedCircTest)))
writeDataTable(wb, sheet = 2, x=RNAse_RenrichedCircTest)



# Back-Splice junctions from circTest / FDR 1%

colnames(circTestSummary) <- c( "Chr",
                                "Start",
                                "End",
                                "Gene",
                                "JunctionType",
                                "Strand",
                                "Start.End.Region",
                                "OverallRegion",
                                "sig_p"
                                )

addWorksheet(wb, sheetName = "Other BSJ FDR 1%")
writeDataTable( wb,
              sheet = 3,
              x=circTestSummary[grep("not_annotated",circTestSummary[,"Gene"]),]
              )

new_sheet <- splicedExonDF[order(splicedExonDF[,"Pval"]),]
new_sheet <- merge(new_sheet,topSplicedGenesMartDCC,by.x=2,by.y=2)
new_sheet <- new_sheet[order(new_sheet[,"Pval"]),]

addWorksheet(wb, sheetName = "Exon events")
writeDataTable( wb,
              sheet = 4,
              x=new_sheet
              )

# close workbook
saveWorkbook(wb, paste(baseDir,"diff_exon_enrichment.xlsx",sep=""), overwrite = TRUE)

message("Writing additional CSV output")

# write simple csv files
write.table( mainTable[order(mainTable[,"FDR"]),],
          file=paste(baseDir,"exon_enrichment.csv",sep=""), sep ="|")
write.table( RNAse_RenrichedCircTest,
          file=paste(baseDir,"bsj_enrichment.csv",sep=""), sep ="|")

message("Exon analysis finished")

## Plotting section
#
# plotSubset <- mainTable[order(mainTable[,"FDR"]),]
# plotSubset <- subset(plotSubset[,"GeneID"],plotSubset[,"RNaseR_enriched"]==1);
#
# message("Plotting top enriched RNaseR enriched circles")
#
#
# txdb <- makeTxDbFromGFF(arg_gtf_file, format="gtf")
# g<-genes(txdb)
#
# model <- exonsBy(txdb, by = "tx")
# exons <- exons(txdb)
#
# message("Generating PDF plots")
# sink(file("/dev/null", open = "wt"), type = c("output", "message"))
#
# # for the top30 enriched genes
# for (current_gene in topEnriched)
# {
# # overlap current enriched gene with exons from txdb
# currentExon <- subsetByOverlaps(exons, g[current_gene])
#
# reducedExon <- reduce(currentExon) #reduce could work here
# print(reducedExon)
# exonIDs <- grep(current_gene, names(diffSplicedGenes$exon.p.value), value = T)
#
# exonValues <- subset(geneBaseTable, geneBaseTable[, "e_id"] %in% exonIDs)
# exonValues <- merge(
#   exonValues,
#   data.frame(
#     exonID = grep(current_gene, names(diffSplicedGenes$exon.p.value), value = T),
#     fc = diffSplicedGenes$coefficients[grep(current_gene, names(diffSplicedGenes$exon.p.value))]
#   ),
#   by.x = 1,
#   by.y = 1
# )
#
# values(reducedExon)$RNAseR_foldchange <- rep(NA, length(reducedExon))
#
# for (z in 1:nrow(exonValues))
# {
#   targ <- which((start(reducedExon) == exonValues[z, "start"]) & (end(reducedExon) == exonValues[z, "end"]))
#   values(reducedExon)$RNAseR_foldchange[targ] <- exonValues[z, "fc"]
#
#   if (!isTRUE(targ > 0)) {
#     for (start_tmp in -10:10){
#       for (stop_tmp in -10:10){
#           targ <- which((start(reducedExon) == exonValues[z, "start"]+start_tmp) & (end(reducedExon) == exonValues[z, "end"]+stop_tmp))
#             if (isTRUE(targ > 0)) {
#               values(reducedExon)$RNAseR_foldchange[targ] <- exonValues[z, "fc"]
#             }
#         }
#       }
#     }
#   }
#
# values(reducedExon)$significant <- rep("FALSE", length(reducedExon))
#
# multiExonOverlap <- as.matrix(findOverlaps(genomicRanges, reducedExon))
#
# if (nrow(multiExonOverlap) > 0)
# {
#   values(reducedExon[multiExonOverlap[, 2]])$significant <- "TRUE"
# }
#
# # lower part with transcript landscape around the selected gene
# plot <- autoplot(
#   txdb,
#   g[current_gene],
#   main = current_gene,
#   label = FALSE,
#   legend = FALSE,
#   xlab = "Position in genome",
#   ylab = "Transcripts",
#   axis.text.x = "axis x",
#   axis.text.y = "axis y",
#   base_size = 20,
#    label.size = 3
# )+ theme(legend.position = "none")
#
# plotRangesLinkedToData(
#   reducedExon,
#   main = current_gene,
#   stat.y = c("RNAseR_foldchange"),
#   stat.ylab = "RNaseR fold change",
#   annotation = list(plot),
#   sig = "significant",
#   sig.col = c("#d7191c", "#2c7bb6"),
#   width.ratio = 0.6,
#   theme.stat = theme_bw(),
#   theme.align = theme_bw(),
#   linetype = 6
# )
#
# # an a4 page landscape
# ggsave(
#   width = 10,
#   height = 8,
#   paste(baseDir, current_gene, ".pdf", sep = ""),
#   title = paste("Enrichment plot: ", current_gene , sep = "")
# )
# }
# sink(file=NULL)
#
# message("Finished plotting, exiting")
