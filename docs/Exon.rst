Alternative exon module
********************************************************

The circtools exon usage module was implemented to detect and analyze differential exon usage in circRNA data sets. As an example, the module may list exons as significant differentially spliced in RNaseR treated sample compared to untreated samples, thus pointing out exons that may be part of a circRNA. 

``circtools exon`` requires mapped sequencing reads that are passed to `StringTie <https://ccb.jhu.edu/software/stringtie/>`_ in order to generate data readable by the `ballgown <https://bioconductor.org/packages/release/bioc/html/ballgown.html>`_ R package. 


Required tools and packages
--------------------------------

``exon`` depends on R and a few additional R packages, namely

* ballgown
* edgeR
* ggbio
* ggfortify
* openxlsx
* GenomicRanges
* GenomicFeatures

The ``exon`` circtools module as well as all R dependencies are automatically installed during the circtools installation procedure.


General usage
--------------

A call to ``circtools exon --help`` shows all available command line flags:

.. code-block:: bash

    usage: circtools [-h] -d DCC_DIR -l CONDITION_LIST -c CONDITION_COLUMNS -g
                     GROUPING -r REPLICATES -b BALLGOWN_DATA -G GTF_FILE -C
                     CIRCTEST_FILE [-H HAS_HEADER] [-o OUTPUT_DIRECTORY]
                     [-n OUTPUT_PREFIX]
    
    circular RNA exon usage analysis
    
    optional arguments:
      -h, --help            show this help message and exit
    
    Required:
      -d DCC_DIR, --DCC DCC_DIR
                            Path to the detect/DCC data directory
      -l CONDITION_LIST, --condition-list CONDITION_LIST
                            Comma-separated list of conditions which should be
                            compared E.g. "RNaseR +","RNaseR -"
      -c CONDITION_COLUMNS, --condition-columns CONDITION_COLUMNS
                            Comma-separated list of 1-based column numbers in the
                            detect/DCC output which should be compared; e.g.
                            10,11,12,13,14,15
      -g GROUPING, --grouping GROUPING
                            Comma-separated list describing the relation of the
                            columns specified via -c to the sample names specified
                            via -l; e.g. -g 1,2 and -r 3 would assign sample1 to
                            each even column and sample 2 to each odd column
      -r REPLICATES, --replicates REPLICATES
                            Comma-separated list describing the relation of the
                            samples specified via -g to the sample names specified
                            via -l; e.g. -g 1,2 and -r 3 would assign sample1 to
                            each even column and sample 2 to each odd column
      -b BALLGOWN_DATA, --ballgown-data BALLGOWN_DATA
                            Path to the ballgown data directory
      -G GTF_FILE, --gtf-file GTF_FILE
                            Path to the GTF file containing the employed genome
                            annotation
      -C CIRCTEST_FILE, --circtest-output CIRCTEST_FILE
                            Path to the CircTest CSV file containing the CircTest
                            results
    
    Additional options:
      -H HAS_HEADER, --has-header HAS_HEADER
                            Do the CircTest result files have a header? [Default:
                            No]
    
    Output options:
      -o OUTPUT_DIRECTORY, --output-directory OUTPUT_DIRECTORY
                            The output directory for files created by circtest
                            [Default: .]
      -n OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                            The output name (prefix) for files created by circtest
                   

Sample call
^^^^^^^^^^^^
.. code-block:: bash

    circtools exon -d /home/tjakobi/work/projects/circRNA/encode_paper/DCC_CD_07_2017/ -r 1,1,2,2,3,3 -l M,P,M,P,M,P -c 10,11,12,13,14,15 -g1,2,1,2,1,2 -C /home/tjakobi/work/projects/circRNA/encode_paper/circtest/fdr_0.05_default/k562_enrichment_total.csv -b /mnt/misc/stringtie_latest/ -G /mnt/big_data/genomes/GRCh38_85/GRCh38.85.gtf -o /home/tjakobi/work/projects/circRNA/encode_paper/circtools_exon/new_06_2018/k562/


Here we have the DCC data located in the folder ``DCC/``, the STAR mapping are stored in ``mapping/``, the experiment had 4 conditions, listed via ``-l HepG2-,HepG2+,K562-,K562+``, the samples in the DCC data file are sorted in the the order specified via ``-g 1,2,1,2,1,2,3,4,3,4,3,4``.

.. code-block:: bash

    Using R version 3.5.0 [/usr/bin/Rscript]
    Loading required packages
    Done loading packages
    Loading CircRNACount
    Loading CircCoordinates
    Starting ballgown processing
    Sun Jun 17 21:17:47 2018
    Sun Jun 17 21:17:47 2018: Reading linking tables
    Sun Jun 17 21:17:48 2018: Reading intron data files
    Sun Jun 17 21:17:52 2018: Merging intron data
    Sun Jun 17 21:17:54 2018: Reading exon data files
    Sun Jun 17 21:18:00 2018: Merging exon data
    Sun Jun 17 21:18:02 2018: Reading transcript data files
    Sun Jun 17 21:18:05 2018: Merging transcript data
    Wrapping up the results
    Sun Jun 17 21:18:05 2018
    Preparing necessary data structures
    Setting treatment and conditions
    Found 11031 multi exon genes
    Found 1574 single exon genes
    Starting dispersion estimation
    Fitting model...
    Writing bed files...
    Writing DCC prediction BED file
    Reading and integrating CircTest results
    Writing back splice junction enriched BED file
    Writing Excel file
    Writing additional CSV output
    Exon analysis finished


``circtools`` takes some time to process the data and prints out information on its progress.


Output produced by ``circtools exon``
-----------------------------------------

exon_analysis_bsj_enrichment.csv
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
circRNA-centric view of the exon results in CSV format. Shown are significantly enriched circRNAs merged with the results from the ballgown package.

exon_analysis_exon_enrichment.csv
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Exon-centric view of the exon results in CSV format. Shown are differentially spliced exons merged with the circRNA detection and circtest step.

exon_analysis_diff_exon_enrichment.xlsx
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
An xlsx Excel file containing 4 work sheets:

* Exon FDR 1% (ballgown): differentially spliced exons, 1% FDR
* enriched BSJ FDR 1% (CircTest): enriched circRNAs, 1% FDR
* Other BSJ FDR 1%: non-annotated circRNAs
* Exon events: all exons 

exon_analysis_dcc_bsj_enriched_track.bed
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A BED file with containing *only* circRNAs predicted by the ``circtools detect`` module that **also** pass the `circtools circtest`` statistical test. Can be displayed in all common visualization tools like IGV. 

exon_analysis_dcc_predictions_track.bed
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A BED file with containing *all* circRNAs predicted by the ``circtools detect`` module. Can be displayed in all common visualization tools like IGV. 

exon_analysis_exon_fc_track.bedgraph
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A BEDgraph file with fold changes of all differentially spliced exons. Can be displayed in all common visualization tools like IGV. 

exon_analysis_exon_pval_track.bedgraph
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A BEDgraph file with p-values of all differentially spliced exons. Can be displayed in all common visualization tools like IGV. 

