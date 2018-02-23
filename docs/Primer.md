circtools: primer design module
================
#### 2017-06-21 /  [Alexey Uvarovskii](https://github.com/alexey0308)

Introduction
------------

The package aims to simplify analysis of the RNA-seq experiments for circular RNA detection and quantification. Once a candidate list of possible circulirised RNA is provided (e.g. after analysis with splice-aware alignment tools like STAR and circRNA detection software like DCC or CIRI), it may be necessary to validate such transcripts by PCR.

We implemented a set of functions which could assist to design such validation experiments. A biologist might be interested to know

-   which linear trabscripts are expressed and at which level
-   which transcripts have common exons with the predicted circRNA
-   how to design primers for the chosen linear and circular transcripts.

A proposed workflow
-------------------

1.  Prepare the input data set:
    -   the gene model: a GTF/GFF file or EnsDb objects
    -   transcripts counts
    -   splice junction coordinates for the circular candidates.

2.  Generate sequences of the exons around the splice junctions.
3.  Plot the gene model to see relation of linear and circular transcripts.
4.  Design optimal primers for the circular and linear transcripts of interest.
5.  Plot the gene model and the primers to validate the design.

Transcript plot
---------------

One would like to know, how the predicted circRNA relate to the linear transcripts. In order to see, which transcripts are expressed and what is their gene model (i.e. exon content), we implemented a plotting function. It shows a block structure for every annotated transcript, a coordinate range which is covered by the circRNA candidate and plots read counts for the transcripts if they are provided.

Annotation source
-----------------

Coordinates and types of genomic features must be known in advance, and it can be provided in the form of `GTF/GFF` file or `EnsDb` objects. In the frame of the package we encourage one to use the Ensemble annotation files or R packages, which can be downloaded from the [Ensembl site](http://www.ensembl.org) or installed via the Bioconductor ecosystem:

``` r
source("https://bioconductor.org/biocLite.R")
biocLite("EnsDb.Hsapiens.v86")
```

In this example we are using the `ensembldb` package for annotation:

``` r
suppressPackageStartupMessages(
  library("EnsDb.Hsapiens.v86")
)
db <- EnsDb.Hsapiens.v86
```

In the first step, the helper function create an object, which keep relations between the annotated transcripts and coordinates for the predicted circRNAs. We assume that a user have a `data.frame` or `GRanges` object with splice junctions coordinates from which circRNAs are derived. The input must contain the chromosome, start, end and the strand of the splice junction.

Generate mock data
------------------

Here we create a splice junction GRanges:

``` r
geneName <- "BCL6"
circs <- createCirc(geneName, db)
```

    ## Warning: 'TxidFilter' is deprecated.
    ## Use 'TxIdFilter' instead.
    ## See help("Deprecated")

    ## Warning: 'TxidFilter' is deprecated.
    ## Use 'TxIdFilter' instead.
    ## See help("Deprecated")

``` r
circs
```

    ## GRanges object with 2 ranges and 1 metadata column:
    ##       seqnames                 ranges strand |                  sjId
    ##          <Rle>              <IRanges>  <Rle> |           <character>
    ##   [1]        3 [187734869, 187745727]      - | 3:187734869-187745727
    ##   [2]        3 [187734869, 187737088]      - | 3:187734869-187737088
    ##   -------
    ##   seqinfo: 1 sequence from GRCh38 genome

Let us simulate some numbers for read counts for transcripts of the BCL6 gene:

``` r
counts <- makeCounts(geneName, db)
tail(counts)
```

    ##                 id count
    ## 6  ENST00000450123     0
    ## 7  ENST00000470319     0
    ## 8  ENST00000479110   878
    ## 9  ENST00000480458   910
    ## 10 ENST00000496823     0
    ## 11 ENST00000621333     0

The workflow entry point
------------------------

Using the splice junction table and annotation object, prepare the `CircData` object:

``` r
suppressPackageStartupMessages(
  library(circtools)
)
circData <- CircData(db, circs)
```

Plot the gene model, circular transcripts and read counts for the BCL6 gene:

``` r
bcl6EnsId <- circData$sjGeneIds
bcl6EnsId
```

    ## [1] "ENSG00000113916"

``` r
plotCirc(circGenes = bcl6EnsId,
         circData = circData,
         counts = counts, 
         opts = list(normalise = FALSE))
```

![](img/unnamed-chunk-9-1.png)

Report sequences of the splice junction exons
---------------------------------------------

Most probably, an experimentalist is interested to obtain the sequencies of the splice junction exons. The sequencies will be used for the following primer design to validate the discovered transcripts using the PCR.

To achieve it, besides the gene model annotation, one needs a fasta file or an R package with the corresponding genome sequence. a BSgenome Bioconductor package can be used:

``` r
suppressPackageStartupMessages(
  library(BSgenome.Hsapiens.NCBI.GRCh38))
bsg <- BSgenome.Hsapiens.NCBI.GRCh38
```

Several exons with the same start but different lengths can be included in the annotation. By default, all described exons, which start or end at the position of the circular splice junction will be reported. It is possible to include the shortes or the longest sequence by setting the `type` argument:

``` r
# for all exons use
exSeqAll <- getExonSeqs(circData = circData, bsg = bsg, type = "all")
exShortesSeq <- getExonSeqs(circData = circData, bsg = bsg, type = "shortest")
exShortesSeq[['3:187734869-187737088']]
```

    ## GRanges object with 2 ranges and 5 metadata columns:
    ##       seqnames                 ranges strand |         exon_id
    ##          <Rle>              <IRanges>  <Rle> |     <character>
    ##   [1]        3 [187734869, 187734882]      - | ENSE00002535122
    ##   [2]        3 [187736090, 187737088]      - | ENSE00001666929
    ##               gene_id                  sjId        side
    ##           <character>           <character> <character>
    ##   [1] ENSG00000113916 3:187734869-187737088        left
    ##   [2] ENSG00000113916 3:187734869-187737088       right
    ##                           seq
    ##                <DNAStringSet>
    ##   [1]          AAGCAAGGCATTGG
    ##   [2] TCTCATTGAC...TGCTCATTTG
    ##   -------
    ##   seqinfo: 1 sequence from GRCh38 genome

Experimentalists might be interested in obtaining a list with exon sequences and their coordinates for every circular splice junction.

Having the list of sequences from the `getExonSeqs` function and the `CircData` object, let us create an HTML report to present it to the biologists:

``` r
reportCircs(exSeq = exShortesSeq, file = "report.html")
```

The resulting html file includes the information on the exon coordinates and sequences for the following primer design.

![The HTML report with splice junction exon sequences.](img/report.png)

Design and validate primers
---------------------------

To get *in silico* optimized primer sequences, one needs simply to invoke `designPrimers` function on the splice junction exons object:

``` r
primers <- designPrimers(exSeq = exShortesSeq, db = db, bsg = bsg)
```

The result is a list with an item for every splice junction. There are two records: `primers` and `products`. Every item consists of a list of primers for possible circular transcripts: if there several exons, which correspons to the same splice junction, all possible combinations of their pairs will be used for primer design.

The priducts are

``` r
str(primers$products)
```

    ## List of 2
    ##  $ 3:187734869-187737088:'data.frame':   1 obs. of  5 variables:
    ##   ..$ sjId      : chr "3:187734869-187737088"
    ##   ..$ upExonId  : chr "ENSE00002535122"
    ##   ..$ downExonId: chr "ENSE00001666929"
    ##   ..$ circSeq   : chr "AAGCAAGGCATTGGTCTCATTGACAGCCCTGCTCCTTGGAGATTGTTTTTGTGGGTAGTCTTGTGTGTGGCATTGGTGGAATGGCTGAATCTAGGAGACGCGGCGTGTCCA"| __truncated__
    ##   ..$ seqId     : chr "3:187734869-187737088"
    ##  $ 3:187734869-187745727:'data.frame':   1 obs. of  5 variables:
    ##   ..$ sjId      : chr "3:187734869-187745727"
    ##   ..$ upExonId  : chr "ENSE00002535122"
    ##   ..$ downExonId: chr "ENSE00001564372"
    ##   ..$ circSeq   : chr "AAGCAAGGCATTGGATACCATCGTCTTGGGCCCGGGGAGGGAGAGCCACCTTCAGGCCCCTCGAGCCTCGAACCGGAACCTCCAAATCCGAGACGCTCTGCTTATGAGGAC"| __truncated__
    ##   ..$ seqId     : chr "3:187734869-187745727"

and the primers

``` r
primers$primers$`3:187734869-187737088`
```

    ## GRangesList object of length 1:
    ## $3:187734869-187737088 
    ## GRanges object with 3 ranges and 6 metadata columns:
    ##       seqnames                 ranges strand |                    seq
    ##          <Rle>              <IRanges>  <Rle> |            <character>
    ##   [1]        3 [187734869, 187734876]      - | GGCATTGGTCTCATTGACAGCC
    ##   [2]        3 [187737075, 187737088]      - | GGCATTGGTCTCATTGACAGCC
    ##   [3]        3 [187736974, 187736990]      - |      TCTGGACACGCCGCGTC
    ##              type efficiency                 seqId          upExon
    ##       <character>  <numeric>           <character>     <character>
    ##   [1]     forward  0.9434878 3:187734869-187737088 ENSE00002535122
    ##   [2]     forward  0.9434878 3:187734869-187737088 ENSE00002535122
    ##   [3]     reverse  0.9853872 3:187734869-187737088 ENSE00002535122
    ##              downExon
    ##           <character>
    ##   [1] ENSE00001666929
    ##   [2] ENSE00001666929
    ##   [3] ENSE00001666929
    ## 
    ## -------
    ## seqinfo: 1 sequence from an unspecified genome; no seqlengths

``` r
circ <- "3:187734869-187737088"
plotCirc(sjIds = circ,
         #circGenes = bcl6EnsId,
         circData = circData,
         counts = counts, 
         primers = primers$primers[[circ]],
         opts = list(normalise = FALSE))
circ <-"3:187734869-187745727"
plotCirc(sjIds = circ,
         #circGenes = bcl6EnsId,
         circData = circData,
         counts = counts, 
         primers = primers$primers[[circ]],
         opts = list(normalise = FALSE))
```

![](img/unnamed-chunk-16-1.png)![](img/unnamed-chunk-16-2.png)

### Filter by counts and easy view

Sometimes it is cleaner to keep only expressed transcripts. One can specify a threshold for read count in `countThres` argument. In addition, by default, all the coordinated used for plotting are transformed for easier interpretation of relative position. It can be turned off in `opts$normalise = FALSE`.

``` r
circ <- "3:187734869-187737088"
plotCirc(sjIds = circ,
         #circGenes = bcl6EnsId,
         circData = circData,
         counts = counts, 
         primers = primers$primers[[circ]],
         countThres = 1,
         opts = list(normalise =TRUE))
circ <- "3:187734869-187745727"
plotCirc(sjIds = circ,
         #circGenes = bcl6EnsId,
         circData = circData,
         counts = counts, 
         primers = primers$primers[[circ]],
         countThres = 1,
         opts = list(normalise =TRUE))
```

![](img/unnamed-chunk-17-1.png)![](img/unnamed-chunk-17-2.png)

Session
-------

``` r
sessionInfo()
```

    ## R version 3.4.0 (2017-04-21)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 16.04.2 LTS
    ## 
    ## Matrix products: default
    ## BLAS: /usr/lib/openblas-base/libblas.so.3
    ## LAPACK: /usr/lib/libopenblasp-r0.2.18.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=de_DE.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats4    parallel  methods   stats     graphics  grDevices utils    
    ## [8] datasets  base     
    ## 
    ## other attached packages:
    ##  [1] BSgenome.Hsapiens.NCBI.GRCh38_1.3.1000
    ##  [2] BSgenome_1.44.0                       
    ##  [3] rtracklayer_1.36.0                    
    ##  [4] Biostrings_2.44.0                     
    ##  [5] XVector_0.16.0                        
    ##  [6] circtools_0.0.0.9000                  
    ##  [7] EnsDb.Hsapiens.v86_2.1.0              
    ##  [8] ensembldb_2.0.1                       
    ##  [9] AnnotationFilter_1.0.0                
    ## [10] GenomicFeatures_1.28.0                
    ## [11] AnnotationDbi_1.38.0                  
    ## [12] Biobase_2.36.2                        
    ## [13] GenomicRanges_1.28.1                  
    ## [14] GenomeInfoDb_1.12.0                   
    ## [15] IRanges_2.10.1                        
    ## [16] S4Vectors_0.14.1                      
    ## [17] BiocGenerics_0.22.0                   
    ## [18] rmarkdown_1.5                         
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] SummarizedExperiment_1.6.1    lattice_0.20-35              
    ##  [3] DECIPHER_2.4.0                htmltools_0.3.6              
    ##  [5] yaml_2.1.14                   interactiveDisplayBase_1.14.0
    ##  [7] XML_3.98-1.7                  DBI_0.6-1                    
    ##  [9] BiocParallel_1.10.1           matrixStats_0.52.2           
    ## [11] GenomeInfoDbData_0.99.0       stringr_1.2.0                
    ## [13] zlibbioc_1.22.0               ProtGenerics_1.8.0           
    ## [15] memoise_1.1.0                 evaluate_0.10                
    ## [17] knitr_1.16                    biomaRt_2.32.0               
    ## [19] httpuv_1.3.3                  BiocInstaller_1.26.0         
    ## [21] Rcpp_0.12.11                  xtable_1.8-2                 
    ## [23] backports_1.1.0               DelayedArray_0.2.2           
    ## [25] mime_0.5                      Rsamtools_1.28.0             
    ## [27] AnnotationHub_2.8.1           digest_0.6.12                
    ## [29] stringi_1.1.5                 shiny_1.0.3                  
    ## [31] rprojroot_1.2                 grid_3.4.0                   
    ## [33] tools_3.4.0                   bitops_1.0-6                 
    ## [35] magrittr_1.5                  RCurl_1.95-4.8               
    ## [37] lazyeval_0.2.0                RSQLite_1.1-2                
    ## [39] Matrix_1.2-10                 httr_1.2.1                   
    ## [41] R6_2.2.1                      GenomicAlignments_1.12.1     
    ## [43] compiler_3.4.0
