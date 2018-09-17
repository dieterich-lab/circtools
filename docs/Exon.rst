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
                     CIRCTEST_FILE [-s {mm,rn,hs}] [-H HAS_HEADER]
                     [-o OUTPUT_DIRECTORY] [-n OUTPUT_PREFIX]

    circular RNA exon usage analysis

    optional arguments:
      -h, --help            show this help message and exit

    Required:
      -d DCC_DIR, --DCC DCC_DIR
                            Path to the detect/DCC data directory
      -l CONDITION_LIST, --condition-list CONDITION_LIST
                            Comma-separated list of conditions which should be
                            comparedE.g. "RNaseR +","RNaseR -"
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
      -s {mm,rn,hs}, --species {mm,rn,hs}
                            Organism of the study (used for primer BLASTing), rn =
                            Rattus norvegicus, mm = Mus musculus, hs = Homo
                            sapiens

    Additional options:
      -H HAS_HEADER, --has-header HAS_HEADER
                            Do the CircTest result files have a header? [Default:
                            No]

    Output options:
      -o OUTPUT_DIRECTORY, --output-directory OUTPUT_DIRECTORY
                            The output directory for files created by circtools
                            [Default: .]
      -n OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                            The output name (prefix) for files created by
                            circtools [Default: exon_analysis]

Generating necessary ballgown data files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to perform per-exon analyses, the circtools exon module requires additional data generated with `StringTie <https://ccb.jhu.edu/software/stringtie/>`_. The tutorial assumes, that stringtie has been installed and is available via the ``$PATH`` environment.

.. code-block:: bash

    # download wrapper for Stringtie
    wget https://raw.githubusercontent.com/dieterich-lab/bioinfo-scripts/master/slurm_stringtie.sh
    chmod 755 slurm_stringtie.sh
    mkdir stringtie/

    # obtain the annotation of the mouse genome for splice junctions
    wget ftp://ftp.ensembl.org/pub/release-90/gtf/mus_musculus/Mus_musculus.GRCm38.90.gtf.gz
    gzip -d Mus_musculus.GRCm38.90.gtf.gz

    cd star/

    # run stringtie for all samples
    parallel ../slurm_stringtie.sh {}/Aligned.noS.bam ../Mus_musculus.GRCm38.90.gtf ../stringtie/{}_StringTieBallgown/ ::: ALL*


circtools exon module call
^^^^^^^^^^^^^^^^^^^^^^^^^^^

After generating all necessary input data, the circtools exon module can now be run via the following command:

.. code-block:: bash

    circtools exon -d 01_detect/ -r1,1,2,2,3,3,4,4 -l minus,plus -c 4,5,6,7,8,9,10,11 -g1,2,1,2,1,2,1,2 -C 04_circtest/circtest.csv -b ../stringtie/ -G ../Mus_musculus.GRCm38.90.gtf -o 05_exon/ -s mm


Here we have the DCC data located in the folder ``01_detect/``, the stringtie data are stored in ``../stringtie/``, the experiment had 2 conditions, listed alternating via ``-l minus,plus``, the samples in the circtools detect data file are sorted in the the order specified via `` -g1,2,1,2,1,2,1,2`` and columns 10-15 are used for the analysis, as specified via ``-c 4,5,6,7,8,9,10,11``. The genome annotation has to be supplied with the ``-G ../Mus_musculus.GRCm38.90.gtf`` flag. Significantly enriched circRNAs from the ``circtest`` module have to be passed via ``-C 04_circtest/circtest.csv``, the species for the internal gene ID conversion has been set via ``-s mm`` to mouse, the output will be stored in ``-o 05_exon/``.

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
