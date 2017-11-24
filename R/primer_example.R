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

library("circtools")
library("EnsDb.Hsapiens.v86")

db <- EnsDb.Hsapiens.v86

library(BSgenome.Hsapiens.NCBI.GRCh38)
bsg <- BSgenome.Hsapiens.NCBI.GRCh38


gr <- GRanges(
    seqnames = Rle(c("1")),
    ranges = IRanges(1223244, end = 1223968, names = "SDF4"))
gr
circData <- CircData(db, gr)

circData$sjGeneIds

exSeqAll <- getExonSeqs(circData = circData, bsg = bsg, type = "all")
exShortesSeq <- getExonSeqs(circData = circData, bsg = bsg, type = "shortest")
reportCircs(exSeq = exShortesSeq, file = "report.html")
primers <- designPrimers(exSeq = exShortesSeq, db = db, bsg = bsg)

bcl6EnsId <- circData$sjGeneIds
bcl6EnsId
pdf(height= 8.2, width=11.69)
plotCirc(circGenes = bcl6EnsId,
         circData = circData,
    primers = primers$primers[[gr$sjId]],
         opts = list(normalise = TRUE))

dev.off()