#!/usr/bin/env Rscript

suppressMessages(library(primex))

args <- commandArgs(trailingOnly = TRUE)

exon1 <- args[1]
exon2 <- args[2]
name <- args[3]

seqOpts <- seqSettings(
  seqId = name,
  seq = c(exon2,exon1))

seqOpts$SEQUENCE_OVERLAP_JUNCTION_LIST = NULL

seqOpts$PRIMER_NUM_RETURN = 10

sink("/dev/null")
productSize(seqOpts,c(50,160))
primers <- design(seqOpts, returnStats = FALSE)
sink()

primers$primers
