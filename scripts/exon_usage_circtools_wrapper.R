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

# suppressMessages(library(ballgown))
# suppressMessages(library(edgeR))
# suppressMessages(library(biomaRt))
# suppressMessages(library(ggbio))
# suppressMessages(library(ggfortify))
# suppressMessages(library(openxlsx))
# suppressMessages(library(GenomicRanges))
# suppressMessages(library(GenomicFeatures))

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
arg_ballgown_directory <- args[12] # string
arg_gtf_file <- args[13] # string


# read sub directories containing the ballgown runs
ballgownRuns <- list.files(arg_ballgown_directory)


print(ballgownRuns)

# set ballgown data input dir
ballgownData <- "../ballgown_data/"


q()
# Read in circRNA predictions by DCC from here
DCCdataDir <- ".";

runAnalysis = function(sampleName, mockColumns, treatedColumns) {

  # output directory
  baseDir <- paste("./ballgown/",sampleName,"/", sep="")

  mock <- unique (grep(paste(mockColumns,collapse="|"), ballgownFile$V1, value=TRUE))
  treated <- unique (grep(paste(treatedColumns,collapse="|"), ballgownFile$V1, value=TRUE))

  sampleNames<- c(
                    paste(sampleName,"- 1", sep=""),
                    paste(sampleName,"- 2", sep=""),
                    paste(sampleName,"- 3", sep=""),
                    paste(sampleName,"+ 1", sep=""),
                    paste(sampleName,"+ 2", sep=""),
                    paste(sampleName,"+ 3", sep="")
                  )

  bg.dccDF <- data.frame( id=sampleNames,
                          group=c(
                                     c(
                                       rep(paste(sampleName,"-", sep=""),3),
                                       rep(paste(sampleName,"+", sep=""),3)
                                       )
                                  )
                        )

  bg <- ballgown(c(treated,mock),verbose=TRUE)

  pData(bg) <- bg.dccDF
  #save(bg, file='bg.rda')

  # load('bg.rda')
  bg
  print("preparing and rearranging some tables")
  # 'all' -> meas='all', all measurements

  whole_exon_table = eexpr(bg, 'all') # eexpr -> exon level
  whole_intron_table = iexpr(bg, 'all') # iexpr -> intron level

  t2g<-indexes(bg)$t2g # transcript / gene table
  e2t<-indexes(bg)$e2t # exon / transcript table

  e2g<-unique(merge(e2t,t2g,by.x=2,by.y=1)[,2:3]) # exon / gene table merging
  idx<-names(which(table(e2g[,1])==1)) # idx to sort out non-gene exons
  e2g<-subset(e2g,e2g[,1] %in% idx) # exon / gene transcript / merging

  # redo everything based on counts
  # use mrcount -> multi-map-corrected number of reads overlapping the exon/intron
  nonMMexonCount<-whole_exon_table[,c(1,grep("mrcount",colnames(whole_exon_table)))]

  # only exons > 10 reads throught all mappings
  idx<-which(apply(nonMMexonCount[,-1],1,sum)>10)
  nonMMexonCount<-nonMMexonCount[idx,]

  # merge exon couns and ( exon / gene / transcript) table
  e2g.counts<-merge(nonMMexonCount,e2g,by.x=1,by.y=1)

  # extract first and last columns
  geneBaseTable<-e2g.counts[,c(1,ncol(e2g.counts))]

  # combine tables, add chr, strand, start and stop
  # we now have a large table with exon to gene mappings and base annotation
  geneBaseTable<-merge(geneBaseTable,whole_exon_table[,1:5],by.x=1,by.y=1)

  # add exon number to gene name. e.g.: ENSG00000223972 -> ENSG00000223972.8
  geneBaseTable[,1]<-gsub(" ","",apply(geneBaseTable[,c(2,1)],
                          1,paste,sep="",collapse="."))

                          print(head(geneBaseTable))

  # bundle together gene names, gene/exon combos in a data frame
  geneIDs <- data.frame(  GeneID=geneBaseTable[,"g_id"],
                          Gene.Exon=geneBaseTable[,"e_id"])

  print("setting treatment and conditions")

  # 3 RNase+, 3 RNase-
  treatment <- c(
                  rep(paste(sampleName,"-", sep=""),3),
                  rep(paste(sampleName,"+", sep=""),3)
                )

  # 6 samples
  condition <- c(
                    paste(sampleName,"- 1", sep=""),
                    paste(sampleName,"- 2", sep=""),
                    paste(sampleName,"- 3", sep=""),
                    paste(sampleName,"+ 1", sep=""),
                    paste(sampleName,"+ 2", sep=""),
                    paste(sampleName,"+ 3", sep="")
                  )

  batch <- c("R1","R2","R3","R1","R2","R3")


  # remove first and last column
  e2g.minimal=e2g.counts[,-c(1,ncol(e2g.counts))]
   # indices of exons in e2g.minimal table with > 40 counts throughout all samples
  idx<-which(apply(e2g.minimal,1,sum)>4*10)

  # select > 40 count exons from exon / gene table
  e2g.counts<-e2g.counts[idx,];

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

  print("multi exon genes:")
  length(genesWithMultipleExons)

  print("single exon genes:")
  length(genesWithSingleExon)

  # subset the geneID list to only contain multi-exon genes
  geneIDs=subset(geneIDs,geneIDs[,1] %in% genesWithMultipleExons)

  # intersect the minimal e2g list with with the geneIDs to get exon (gene)
  # -> count relation
  e2g.minimal=e2g.minimal[intersect(rownames(e2g.minimal),geneIDs[,2]),]

  # DE analysis with edgeR starts here
  z <- DGEList(counts=e2g.minimal, genes=geneIDs)
  z <- calcNormFactors(z)

  # build the design for edgeR
  design <- model.matrix(~ batch + treatment)

  print(design)

  print("starting dispersion estimation...")
  z <- estimateGLMCommonDisp(z,design)
  z <- estimateGLMTrendedDisp(z,design)
  z <- estimateGLMTagwiseDisp(z,design)

  print("fitting model...")
  fit <- glmFit(z,design)


  print("testing for differential exon usage...")
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
  # save.image()
  # load(".RData")

  # build new table with selected columns
  splicedExonDFfixed=splicedExonDF[,c("chr","start","end","log2FC","Pval")];
  splicedExonDFfixed[,2]<-splicedExonDFfixed[,2]-1 # fix start pos (0-based now)
  splicedExonDFfixed[,3]<-splicedExonDFfixed[,3]-1 # fix stop pos (0-based now)
  splicedExonDFfixed[,5]<-log2(splicedExonDFfixed[,5]) # get log2 for Pval column

  print("writing bed files...")
  ############# Writing BED files ##############

  # write BED track header for FC track
  write(  paste("track type=bedGraph name=\"edgeR_exon_foldChange\"",
           "description=\"Exon fold change over all total RNA data sets\"",
           "visibility=full color=200,100,0 altColor=0,100,200 priority=20",
           sep=""),file=paste(baseDir,"exon_fc_track.bedgraph", sep="/"));

  # write actual content (FCs)
  write.table(
              splicedExonDFfixed[,1:4],
              file=paste(baseDir,"exon_fc_track.bedgraph",sep="/"),
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
             sep=""), file=paste(baseDir,"exon_pval_track.bedgraph", sep="/"));

  # write actual content (PValues)
  write.table(  splicedExonDFfixed[,c(1:3,5)],
                file=paste(baseDir,"exon_pval_track.bedgraph",sep="/"),
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


  # disabled, biomart data cached to speed up script
  print("interfacing biomart... this may take some time...")
  # create biomart object with HS data set
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

  # get gene ID, description and hgcn symbol
  martData=getBM(
                  attributes=c("hgnc_symbol","ensembl_gene_id", "entrezgene", "description"),
                  filter="ensembl_gene_id",
                  values=topSplicedGenes$GeneID,
                  mart=mart,
                  verbose = FALSE
                )

  # print("caching biomart data...")
  # save(martData, mart, file='biomart.rda')

  ## load the cached version of the biomart query
  ## has to be updated when the gene list input changes!
  # load(file='biomart.rda')
  print("done with biomart...")


  # enrich the top spliced gene swith the mart data
  topSplicedGenesMartData=merge(topSplicedGenes, martData, by.x=1, by.y=2)

  # read in circ RNA annotation from DCC
  CircCoordinates <- read.delim(
                                paste(DCCdataDir,"CircCoordinates",sep="/"),
                                header=T,
                                as.is=T
                               )
  # correct again for to 0-based positions
  dccDF<-CircCoordinates[,c(1,2,3,4,5,6)]
  dccDF[,2]<-dccDF[,2]-1
  dccDF[,3]<-dccDF[,3]-1

  # write BED file
  colnames(dccDF)<-c("chr","start","end","GeneName","JType","strand")

  #Print bed file, read counts
  write(paste("track name=\"DCC_predictions\" description=\"Predictions ",
              "over all data sets\" color=\"red\"", sep=""),
               file=paste(baseDir,"dcc_predictions_track.bed",sep="/"));

  write.table(  dccDF,
                file=paste(baseDir,"dcc_predictions_track.bed",sep="/"),
                quote=F,
                sep="\t",
                append=T,
                row.names=F,
                col.names=F
              )

  #genomicRanges <- makeGRangesFromDataFrame(dccDF);
  #do we need to test if enrichment happens inside of BSJ ?!?

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
                            paste(DCCdataDir,"/circTest/circRNA_",
                            sampleName, "_RNaseR_P_signif_1percFDR.csv",sep=""),
                            header=T,
                            as.is=T
                          )
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
  dccDF<-circTestSummary[,2:7]

  # reorder columns
  colnames(dccDF)<-c("chr","start","end","GeneName","JType","strand")

  # Print BED file, read counts
  write(paste("track name=\"DCC_BSJ_enriched\"",
              "description=\"Back-Splice junction enriched over all data ",
              "sets 1% FDR\" color=\"green\"", sep=""),
               file=paste(baseDir,"dcc_bsj_enriched_track.bed",sep="/"));

  write.table(  dccDF,
                file=paste(baseDir,"dcc_bsj_enriched_track.bed",sep="/"),
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
                            )[,-2]
  # sort by PValue
  RNAse_RenrichedCircTest <- RNAse_RenrichedCircTest[order(RNAse_RenrichedCircTest[,"sig_p"]),]

  # create a data frame with all gene names from circtest and mark them
  CircbackSpliceEnrich <- data.frame(
                                      Gene=unique(circTestSummary[,"Gene"]),
                                      RNaseR_enriched=rep(1,length(unique(circTestSummary[,"Gene"])))
                                    )

  # sort the top spliced genes from DCC with mart annotation by FDR
  topSplicedGenesMartDCC <- topSplicedGenesMartDCC[order(topSplicedGenesMartDCC[,"FDR"]),]

  #print out single exon subset in circTest summary
  circTestSummary[unique(singleExonOverlap[,1]),]


  # combine DCC top spliced genes with the BS enriched data frame
  mainTable <- merge( topSplicedGenesMartDCC,
                      CircbackSpliceEnrich,
                      by.x=1,
                      by.y=1,
                      all.x=T
                      )

  # create an excel workbook
  wb <- createWorkbook()

  # exons with 1% FDR from the main table
  addWorksheet(wb, sheetName = "Exon FDR 1% (ballgown)")
  writeDataTable(wb, sheet = 1, x=mainTable[order(mainTable[,"FDR"]),])

  # Enriched Back-Splice junctions / FDR 1%
  addWorksheet(wb, sheetName = "enriched BSJ FDR 1% (CircTest)")
  writeDataTable(wb, sheet = 2, x=RNAse_RenrichedCircTest)

  # Back-Splice junctions from circTest / FDR 1%
  addWorksheet(wb, sheetName = "Other BSJ FDR 1%")
  writeDataTable( wb,
                  sheet = 3,
                  x=circTestSummary[grep("not_annotated",circTestSummary[,"Gene"]),]
                  )

  # close workbook
  saveWorkbook(wb, paste(baseDir,"diff_exon_enrichment.xlsx",sep="/"), overwrite = TRUE)

  # write simple csv files
  write.table( mainTable[order(mainTable[,"FDR"]),],
              file=paste(baseDir,"exon_enrichment.csv",sep="/"), sep ="|")
  write.table( RNAse_RenrichedCircTest,
              file=paste(baseDir,"bsj_enrichment.csv",sep="/"), sep ="|")

  ## Plotting section

  plotSubset <- mainTable[order(mainTable[,"FDR"]),]
  plotSubset <- subset(plotSubset[,"GeneID"],plotSubset[,"RNaseR_enriched"]==1);

  # select the top10 enriched genes
  top30Enriched=RNAse_RenrichedCircTest[1:5,"GeneID"]

  txdb <- makeTxDbFromGFF("/home/tjakobi/work/data/GRCh38_85/GRCh38.85.gtf", format="gtf")
  g<-genes(txdb)

  model <- exonsBy(txdb, by = "tx")
  exons <- exons(txdb)

  # for the top10 enriched genes
  for (n in top30Enriched)
  {
    # overlap current enriched gene with exons from txdb
    currentExon <- subsetByOverlaps(exons, g[n])

    reducedExon = reduce(currentExon) #reduce could work here
    exonIDs <-
      grep(n, names(diffSplicedGenes$exon.p.value), value = T)


    exonValues <-
      subset(geneBaseTable, geneBaseTable[, "e_id"] %in% exonIDs)
    exonValues <- merge(
      exonValues,
      data.frame(
        exonID = grep(n, names(diffSplicedGenes$exon.p.value), value = T),
        fc = diffSplicedGenes$coefficients
        [grep(n, names(diffSplicedGenes$exon.p.value))]
      ),
      by.x = 1,
      by.y = 1
    )

    values(reducedExon)$RNAseR_foldchange <-
      rep(NA, length(reducedExon))

    for (z in 1:nrow(exonValues))
    {
      targ = which((start(reducedExon) == exonValues[z, "start"]) &
                     (end(reducedExon) == exonValues[z, "end"]))
      values(reducedExon)$RNAseR_foldchange[targ] <-
        exonValues[z, "fc"]
    }

    values(reducedExon)$significant <-
      rep("FALSE", length(reducedExon))

    multiExonOverlap <-
      as.matrix(findOverlaps(genomicRanges, reducedExon))

    if (nrow(multiExonOverlap) > 0)
    {
      values(reducedExon[multiExonOverlap[, 2]])$significant <- "TRUE"
    }

    print(reducedExon)

    # lower part with transcript landscape around the selected gene
    plot <- autoplot(
      txdb,
      g[n],
      main = n,
      label = FALSE,
      legend = FALSE,
      xlab = "",
      ylab = ""
    )

    plotRangesLinkedToData(
      reducedExon,
      stat.y = c("RNAseR_foldchange"),
      stat.ylab = "",
      annotation = list(plot),
      sig = "significant",
      sig.col = c("grey", "red"),
      width.ratio = 0.8,
      theme.stat = theme_gray(),
      theme.align = theme_gray(),
      linetype = 3
    )
    # an a4 page landscape
    ggsave(
      width = 11.69,
      height = 8.27,
      paste(baseDir, n, ".pdf", sep = ""),
      title = paste("Enrichment plot: ", n , sep = "")
    )
  }

}

################## SETUP

runAnalysis(
              sampleName = "sample",
              treatedColumns = vector,
              mockColumns = mvector
            )
warnings()
