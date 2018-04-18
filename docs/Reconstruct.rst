Reconstruction module
********************************************************

FUCHS (FUll circular RNA CHaracterization from RNA-Seq) is a Python pipeline designed to fully characterize circular RNAs. It uses a list of circular RNAs and reads spanning the back-splice junction as well as a BAM file containing the mapping of all reads (alternatively of all chimeric reads).

The reads from one circle are extracted by FUCHS and saved in an individual BAM file. Based on these BAM files, FUCHS will detect alternative splicing within the same circle boundaries, summarize different circular isoforms from the same host-gene and generates coverage plots for each circle. It will also cluster circles based on their coverage profile. These results can be used to identify potential false positive circles.

Manual installation
-------------------

Required tools and packages
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

FUCHS dependes on **bedtools (>= 2.27.0)**, **samtools (>= 1.3.1)**,  **Python (> 2.7; pysam>=0.9.1.4, pybedtools>=0.7.8, numpy>=1.11.2, pathos>=0.2.1)**, and **R(>= 3.2.0; amap, Hmisc, gplots)**. All Python an R dependencies will be installed automatically when installing FUCHS. Please make sure to have the correct versions of bedtools and samtools in your ``$PATH``.

Getting FUCHS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Clone the repository and install FUCHS using setup.py:

.. code-block:: bash

  $ git clone git@github.com:dieterich-lab/FUCHS.git

  $ cd FUCHS

  $ python setup.py install --user

  # This will install a FUCHS binary in $HOME/.local/bin/
  # make sure this folder is in your $PATH

  # Check the installation:

  $ FUCHS --help


Usage
--------------
To characterize circRNAs from RNA-seq data you have to:

1. Map RNAseq data from quality checked fastq files with either STAR , BWA, TopHat-Fusion.

2. Detect circRNAs using DCC, CIRI, CIRCfinder or CIRCexplorer depending on the program you used for mapping.

3. Run FUCHS (right now only the combination STAR + DCC has been tested; other setups are under development)


Step by step tutorial
--------------
In this tutorial we will be using HEK293 data available in this repository and use STAR with DCC to detect circular RNAs


1. Mapping of RNA-Seq data
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Map RNA-seq data with `STAR <https://github.com/alexdobin/STAR>`_ (Dobin et al., 2013). Note that ``--alignSJoverhangMin`` and ``--chimJunctionOverhangMin`` should use the same value, to make the circRNA expression and linear gene expression level comparable.
Note that STARlong is not mapping chimeric reads correctly.



* Note: The joined pair mapping should be run first. If the data are paired end, two additional separate mate mappings are recommended This step is not mandatory, but will increase the sensitivity of DCC detection, because it collect small circRNAs which appear with one chimeric junction point at each read mate. If the data is single end, only one mapping step is needed. In this case, PE sequencing data was used.

.. code-block:: bash

  $ STAR --readFilesCommand zcat --runThreadN 18 
         --genomeDir [genome] 
         --outSAMtype BAM SortedByCoordinate 
         --readFilesIn [sample]_1.fastq.gz ([sample]_2.fastq.gz) 
         --outFileNamePrefix [sample]  
         --quantMode GeneCounts 
         --genomeLoad NoSharedMemory 
         --outReadsUnmapped Fastx 
         --outSJfilterOverhangMin 15 15 15 15 
         --alignSJoverhangMin 15 
         --alignSJDBoverhangMin 10 
         --outFilterMultimapNmax 20 
         --outFilterScoreMin 1   
         --outFilterMismatchNmax 999 
         --outFilterMismatchNoverLmax 0.05 
         --outFilterMatchNminOverLread 0.7 
         --alignIntronMin 20 
         --alignIntronMax 1000000 
         --alignMatesGapMax 1000000  
         --chimSegmentMin 15  
         --chimScoreMin 15   
         --chimScoreSeparation 10  
         --chimJunctionOverhangMin 15 
         --twopassMode Basic 
         --alignSoftClipAtReferenceEnds No 
         --outSAMattributes NH HI AS nM NM MD jM jI XS  
         --sjdbGTFfile [annotation].gtf



1.1. Mates separate mapping (optional for PE data)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Note: the mate assignments should be consistent throughout the mapping and circular RNA detection process. In the following case, SamplePairedRead_1.fastq.gz is the first mate which also was the first mate in the STAR call.

.. code-block:: bash

  # remap unmapped reads as single end to obtain double breakpoint fragments
  
  $ gzip sample/Unmapped.out.mate1
  $ mv sample/Unmapped.out.mate1.gz sample/Unmapped_out_mate1.fastq.gz
  $ STAR --readFilesCommand zcat --runThreadN 18 --genomeDir [genome] --outSAMtype BAM SortedByCoordinate --readFilesIn [sample]/Unmapped_out_mate1.fastq.gz --outFileNamePrefix [sample].mate1.  --quantMode GeneCounts --genomeLoad NoSharedMemory --outReadsUnmapped Fastx --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 10 --outFilterMultimapNmax 20 --outFilterScoreMin 1   --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.05 --outFilterMatchNminOverLread 0.7 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000  --chimSegmentMin 15  --chimScoreMin 15   --chimScoreSeparation 10  --chimJunctionOverhangMin 15 --twopassMode Basic --alignSoftClipAtReferenceEnds No --outSAMattributes NH HI AS nM NM MD jM jI XS  --sjdbGTFfile [annotation].gtf

  $ gzip sample/Unmapped.out.mate2
  $ mv sample/Unmapped.out.mate2.gz sample/Unmapped_out_mate2.fastq.gz
  $ STAR --readFilesCommand zcat --runThreadN 18 --genomeDir [genome] --outSAMtype BAM SortedByCoordinate --readFilesIn [sample]/Unmapped_out_mate2.fastq.gz --outFileNamePrefix [sample].mate2.  --quantMode GeneCounts --genomeLoad NoSharedMemory --outReadsUnmapped Fastx --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 10 --outFilterMultimapNmax 20 --outFilterScoreMin 1   --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.05 --outFilterMatchNminOverLread 0.7 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000  --chimSegmentMin 15  --chimScoreMin 15   --chimScoreSeparation 10  --chimJunctionOverhangMin 15 --twopassMode Basic --alignSoftClipAtReferenceEnds No --outSAMattributes NH HI AS nM NM MD jM jI XS  --sjdbGTFfile [annotation].gtf

2. Detection of circRNAs from chimeric.out.junction files with DCC
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Acquiring suitable GTF files for repeat masking
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- It is strongly recommended to specify a repetitive region file in GTF format for filtering.

- A suitable file can for example be obtained through the `UCSC table browser <http://genome.ucsc.edu/cgi-bin/hgTables>`_ . After choosing the genome, a group like **Repeats** or **Variation and Repeats** has to be selected. For the track, we recommend to choose **RepeatMasker** together with **Simple Repeats** and combine the results afterwards.

- **Note**: the output file needs to comply with the GTF format specification. Additionally it may be the case that the names of chromosomes from different databases differ, e.g. **1** for chromosome 1 from ENSEMBL compared to **chr1** for chromosome 1 from UCSC. Since the chromosome names are important for the correct functionality of DCC a sample command for converting the identifiers may be ``sed -i 's/^chr//g' your_repeat_file.gtf``


Preparation of files containing the paths to required ``chimeric.out.junction`` files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* ``samplesheet`` file, containing the paths to the jointly mapped ``chimeric.out.junction`` files

 .. code-block:: bash

  $ cat samplesheet
  /path/to/STAR/sample/joint_mapping/chimeric.out.junction

* ``mate1`` file, containing the paths to ``chimeric.out.junction`` files of the separately mapped first read of paired-end data

 .. code-block:: bash

  $ cat mate2
  /path/to/STAR/sample.mate1/joint_mapping/chimeric.out.junction



* ``mate2`` file, containing the paths to ``chimeric.out.junction`` files of the separately mapped first read of paired-end data

 .. code-block:: bash

  $ cat mate2
  /path/to/STAR/sample.mate2/joint_mapping/chimeric.out.junction


Running DCC
^^^^^^^^^^^^
After performing all preparation steps DCC can now be started:

.. code-block:: bash

  # Call DCC to detect circRNAs, using HEK293 data as example.

  $ DCC @samplesheet \ # @ is generally used to specify a file name
        -mt1 @mate1 \ # mate1 file containing the mate1 independently mapped chimeric.junction.out files
        -mt2 @mate2 \ # mate2 file containing the mate1 independently mapped chimeric.junction.out files
        -D \ # run in circular RNA detection mode
        -R [Repeats].gtf \ # regions in this GTF file are masked from circular RNA detection
        -an [Annotation].gtf \ # annotation is used to assign gene names to known transcripts
        -Pi \ # run in paired independent mode, i.e. use -mt1 and -mt2
        -F \ # filter the circular RNA candidate regions
        -M \ # filter out candidates from mitochondrial chromosomes
        -Nr 2 2 \ minimum number of replicates the candidate is showing in [1] and minimum count in the replicate [2]
        -fg \ # candidates are not allowed to span more than one gene
        -G \ # also run host gene expression
        -A [Reference].fa \ # name of the fasta genome reference file; must be indexed, i.e. a .fai file must be present

  # For details on the parameters please refer to the help page of DCC:
  $ DCC -h

**Notes:**

* By default, DCC assumes that the data is stranded. For non-stranded data the ``-N`` flag should be used.

* Although not mandatory, we strongly recommend to the ``-F`` filtering step

Output files generated by DCC
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The output of DCC consists of the following four files: CircRNACount, CircCoordinates, LinearCount and CircSkipJunctions.

- **CircRNACount:** a table containing read counts for circRNAs detected. First three columns are chr, circRNA start, circRNA end. From fourth column on are the circRNA read counts, one sample per column, shown in the order given in your samplesheet.

- **CircCoordinates:** circular RNA annotations in BED format. The columns are chr, start, end, genename, junctiontype (based on STAR; 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT), strand, circRNA region (startregion-endregion), overall regions (the genomic features circRNA coordinates interval covers).

- **LinearCount:** host gene expression count table, same setup with CircRNACount file.

- **CircSkipJunctions:** circSkip junctions. The first three columns are the same as in LinearCount/CircRNACount, the following columns represent the circSkip junctions found for each sample. circSkip junctions are given as chr:start-end:count, e.g. chr1:1787-6949:10. It is possible that for one circRNA multiple circSkip junctions are found due to the fact the the circular RNA may arise from different isoforms. In this case, multiple circSkip junctions are delimited with semicolon. A 0 implies that no circSkip junctions have been found for this circRNA.


3. Prepare input data for FUCHS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The files  ``chimeric.sam``, ``mate1.chimeric.sam``, and ``mate2.chimeric.sam`` files for FUCHS have to be merged (not necessary if circles were detected using BWA/CIRI)

.. code-block:: bash

  # convert SAM to BAM
  $ samtools view -Sb -o sample sample/Chimeric.out.sam
  $ samtools view -Sb -o sample.1 sample.1/Chimeric.out.sam
  $ samtools view -Sb -o sample.2 sample.2/Chimeric.out.sam

  # sort both BAM files
  $ samtools sort -o sample.sorted.bam sample.bam
  $ samtools sort -o sample.1.sorted.bam sample.1.bam
  $ samtools sort -o sample.2.sorted.bam sample.2.bam

  # create an index for both BAM files
  $ samtools index sample.sorted.bam
  $ samtools index sample.1.sorted.bam
  $ samtools index sample.2.sorted.bam

  # merge both mate BAM files into one new BAM file
  $ samtools merge merged_sample.bam sample.sorted.bam sample.1.sorted.bam sample.2.sorted.bam

  # re-index the newly aggregated BAM file
  $ samtools index merged_sample.bam


4. Running FUCHS
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Run FUCHS to start the pipeline which will extract reads, check mate status, detect alternative splicing events, classify different isoforms, run_primer_design coverage profiles and cluster circRNAs based on coverage profiles

.. code-block:: bash

  # using STAR/DCC Input
  $ FUCHS -r 2 -q 2 -p ensembl -e 2 -T ~/tmp 
	  -D CircRNACount 
	  -J sample/Chimeric.out.junction 
	  -F sample.1/Chimeric.out.junction 
	  -R sample.2/Chimeric.out.junction.fixed 
	  -B merged_sample.sorted.bam 
	  -A [annotation].bed 
	  -N sample 

  # if BWA/CIRI was used, use -C to specify the circIDS list (omit -D, -J, -F and -R)
  # For details on the parameters please refer to the help page of FUCHS:
  $ FUCHS --help

5. Optional FUCHS modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Run the additional module guided_denovo_circle_structure_parallel.py to obtain a more refined circle reconstruction based on intron signals. The circRNA seperated bamfiles (step 2) are the only input needed for the module. If you supply an annotation file, unsupported exons will be reported with a score of 0, if you do not supply an annotation file, unsupported will not be reported.

.. code-block:: bash

  $ guided_denovo_circle_structure_parallel -c 18 -A [annotatation].bed -I FUCHS/output/folder -N sample

  # FUCHS/output/folder corresponds to the output directory of the FUCHS pipeline
  # sample corresponds to your sample name, just as specified for the pipeline


**That's all folks**



Required input data
--------------------

**circIDs:**

==================== ==========================================================================================
 circID               read1,read2,read3
==================== ==========================================================================================
 1:3740233\|3746181  MISEQ:136:000000000-ACBC6:1:2107:10994:20458,MISEQ:136:000000000-ACBC6:1:1116:13529:8356
 1:8495063\|8614686  MISEQ:136:000000000-ACBC6:1:2118:9328:9926
==================== ==========================================================================================


The first column contains the circle id formated as folllowed **chr:start|end**. The second column is a comma separated list of read names spanning the back-splice junction.

**bamfile:** Alignment file produced by any mapper. This file must contain all chimerically mapped reads and may contain also linearly mapped reads.

**bedfile:**

====   ===========    =============     ===================================   =======  ======
Chr      Start            End               Name                               Score   Strand
====   ===========    =============     ===================================   =======  ======
 1      67092175        67093604         NR_075077_exon_0_0_chr1_67092176_r     0       \-
 1      67096251        67096321         NR_075077_exon_1_0_chr1_67096252_r     0       \-
 1      67103237        67103382         NR_075077_exon_2_0_chr1_67103238_r     0       \-
====   ===========    =============     ===================================   =======  ======

Normal BED file in BED6 format. The name should contain a gene name or gene ID and the exon_number. You can specify how the name should be processed using -p (platform), -s (character used to separate name and exon number) and -e (exon_index).


Output produced by FUCHS
--------------

**hek293.alternative_splicing.txt:**

This file summarizes the relationship of different circRNAs derived from the same host-gene.

=============  ============================================================    =========================================  =========   ===========  =============================================
Transcript      circles                                                        same_start                                 same_end    overlapping  within
=============  ============================================================    =========================================  =========   ===========  =============================================
NM_016287	1:20749723-20773610                                            .                                           .          .            .
NM_005095	1:35358925-35361789,1:35381259-35389082,1:35381259-35390098    1:35381259-35389082|1:35381259-35390098,    .          .            .
NM_001291940    1:236803428-236838599,1:236806144-236816543                    .                                           .          .            1:236803428-236838599|1:236806144-236816543,
=============  ============================================================    =========================================  =========   ===========  =============================================

| *Transcript*: Transcript name as defined by the bed-annotation file
| *circles*: Comma-separated list of circRNA ids derived from this transcript
| *same_start*: Comma-seprated list of circRNA pairs separated by |. Pairs in this column share the same start coordinates. A "." indicates that there are no circle pairs that share the same start coordinates.
| *same_end*: Same as *same_start*, only now, circle pairs share the same end coordinates.
| *overlapping*: Comma-seprated list of circRNA pairs separated by |. Pairs in this column share neither start nor end coordinates, but their relation is such that: start.x < start.y && end.x < end.y && start.y < end.x
| *within*: Same as *overlapping*, only now, circle pairs have the follwoing relation: start.x < start.y && end.x > end.y
|

**hek293.exon_counts.bed:**
This file is a bed-formatted file that describes the exon-structure and can be loaded into any genome browser. Each line corresponds to a circRNA.

=====  ============  =============    ============    =============    =======   ======== =========  ======= ===========  ==============  =====================
Chr    Circle Start   Circle  End      Transcript     Num of Reads     Strand      Start   End        Color  Num of Exon  Exon Lengths     Relative Exon Starts
=====  ============  =============    ============    =============    =======   ======== =========  ======= ===========  ==============  =====================
chr1    35358925        35361789        NM_005095       9               \+       35358925 35361789   0,255,0  3           521,61,170      0,2269,2694
chr1    20749723        20773610        NM_016287       4               \-       20749723 20773610   0,255,0  4           159,90,143,159  0,7443,21207,23728
=====  ============  =============    ============    =============    =======   ======== =========  ======= ===========  ==============  =====================

| *Chr*: Chromosome of circRNA
| *Circle Start*: The 5' site of the chimeric junction. This is relative to the reference strand, i.e. start < end! The location is 1-index based
| *Cirlce End*: The 3' site of the chimeric junction. This is relative to the reference strand, i.e. start < end! The location is 0-index based
| *Transcript*: Transcript name as defined by the bed-annotation file
| *Num of Reads* : Number of reads supporting this chimeric junction, in other words, reads that are chimerically mapped to this junction
| *Strand*: Strand of the host-gene
| *Start*: Copied *Circle Start* to stay conform with BED12 format
| *End*: Copied *Circle End* to stay conform with BED12 format
| *Color*: pre defined color the exons will show up in the genome viewer (0,255,0 -> green)
| *Num of Exon*: Number of exons in this circRNA consists of
| *Exon Lengths*: Comma-seprated list of the length of each exon
| *Relative Exon Starts*: Comma-separated list of the relative starting positions of the exons within the circle boundaries.
|
**hek293.exon_counts.txt:**
This file contains similar information as the previous file, just more detailed inforamtion on the exons. Each line corresponds to one exon.

======= =====================  ================ ============  ========== =====  ============   ============= ======= =============   ==============  ===========     ========= ========
sample   circle_id               transcript_id   other_ids       exon_id chr     start           end          strand  exon_length     unique_reads    fragments       number\+ number\-
======= =====================  ================ ============  ========== =====  ============   ============= ======= =============   ==============  ===========     ========= ========
hek293   1:35358925-35361789     NM_005095       NM_005095       2       1       35358924        35359446        \+       522          9               9               4        5
hek293   1:35358925-35361789     NM_005095       NM_005095       3       1       35361193        35361255        \+       62           3               3               1        2
hek293   1:35358925-35361789     NM_005095       NM_005095       4       1       35361618        35361789        \+       171          9               9               4        5
hek293   1:20749723-20773610     NM_016287       NM_016287       3       1       20749722        20749882        \-       160          4               4               4        0
hek293   1:20749723-20773610     NM_016287       NM_016287       4       1       20757165        20757256        \-       91           1               1               1        0
hek293   1:20749723-20773610     NM_016287       NM_016287       5       0       0               0               \0       0            0               0               0        0
hek293   1:20749723-20773610     NM_016287       NM_016287       6       0       0               0               \0       0            0               0               0        0
hek293   1:20749723-20773610     NM_016287       NM_016287       7       1       20770929        20771073        \-       144          1               1               1        0
hek293   1:20749723-20773610     NM_016287       NM_016287       8       1       20773450        20773610        \-       160          4               4               4        0
======= =====================  ================ ============  ========== =====  ============   ============= ======= =============   ==============  ===========     ========= ========

| *sample*: Sample name as specified by the user. This is useful if the user wants to merge files from different samples
| *circle_id*: circRNA-ID. The circleID is formatted to be copy and pasted to a genome browser for easy access
| *transcript_id*: Transcript name as defined by the bed-annotation file. This is the best fitting transcript. i.e. the splicing variants that contains the most exons that are actually covered
| *other_ids*: Alternative Transcript names that are either just as fitting, or contain more or less exons as supported by reads
| *exon_id*: Exon number relative to the host-gene of the circularized exon. One circle may have more than one exon. These will be listed as consecutive lines
| *chr*: Chromosome the circRNA is located on
| *start*: 5' start of the exon, relative to the reference strand, 0-based
| *end*: 3' end of the exon, relative to the reference start, 0-based
| *strand*: Strand of the host-gene
| *exon_length*: Length of the current exon
| *unique_reads*: Number of unique reads associated with the chimeric junction. When the data is paired end, then both ends are considered as separate reads.
| *fragments*: Number of broken fragments aligning to the circle
| *number\+*: Number of reads spanning the chimeric junction on the forward strand
| *number\-*: Number of reads spanning the chimeric junction on the reverse strand (if reads are only from one strand, it could indicate, that there is a sequencing bias.)
|

**hek293.mate_status.txt:**
This output file contains the results of analysing the amount of how often each fragment spans a chimeric junction. A fragment can either span the chimeric junction once (single), only one end spans the junction,
twice (double) both ends span the chimeric junction, or more than twice (undefined).

=====================  ================ =============   ============   ============    ======= ======== ==========
circle_id               transcript_ids  num_reads       min_length      max_length      single  double  undefined
=====================  ================ =============   ============   ============    ======= ======== ==========
1_20749723_20773610     NM_016287       4               790              790             4       0       0
1_35358925_35361789     NM_005095       9               754              754             9       0       0
=====================  ================ =============   ============   ============    ======= ======== ==========

| *circle_id*:
| *transcript_ids*:
| *num_reads*:
| *min_length*:
| *max_length*:
| *single*:
| *double*:
| *undefined*:
|

**hek293.skipped_exons.bed:**

=====  ==============  ============    ==============  ======= ======= =============== ============   ========= ========== ============ =============
Chr     Circle-Start    Circle-End      Transcript      Ratio  Strand   Intron-Start    Intron-End     Color    NumExon\-2 IntronLength RelativeStart
=====  ==============  ============    ==============  ======= ======= =============== ============   ========= ========== ============ =============
chr5    178885614       178931326       NM_030613       60.0    .       178913072       178931236      255,0,0  3          1,146,1      0,30950,45711
chr6    161034259       161049979       NM_001291958    40.0    .       161049332       161049852      255,0,0  3          1,520,1      0,15073,15719
=====  ==============  ============    ==============  ======= ======= =============== ============   ========= ========== ============ =============


**hek293.skipped_exons.txt:**

=====================   ==============  ======================  =============================================   ======================================================================================================================================   =============   ===========
circle_id               transcript_id   skipped_exon            intron                                          read_names                                                                                                                               splice_reads    exon_reads
=====================   ==============  ======================  =============================================   ======================================================================================================================================   =============   ===========
5_178885614_178931326   NM_030613       5:178916564-178916710   set\(\[\(\'5\', 178913072, 178931236\)\]\)      MISEQ:136:000000000-ACBC6:1:2103:10044:24618,MISEQ:136:000000000-ACBC6:1:2115:19571:6931,MISEQ:136:000000000-ACBC6:1:1119:25537:8644     3               5
6_161034259_161049979   NM_001291958    6:161049332-161049852   set\(\[\(\'6\', 161049332, 161049852\)\]\)      MISEQ:136:000000000-ACBC6:1:1113:25288:9067,MISEQ:136:000000000-ACBC6:1:2116:11815:3530                                                  2               5
=====================   ==============  ======================  =============================================   ======================================================================================================================================   =============   ===========


--------------------

**hek293_exon_chain_inferred_12.bed:**

--------------------

**hek293_exon_chain_inferred_6.bed**

--------------------

**hek293:**

1_35358925_35361789_9reads.sorted.bam
1_35358925_35361789_9reads.sorted.bam.bai
1_20749723_20773610_4reads.sorted.bam
1_20749723_20773610_4reads.sorted.bam.bai

--------------------

**hek293.coverage_pictures:**

1_35358925_35361789_NM_005095.png
1_20749723_20773610_NM_016287.png
cluster_means_all_circles.png

--------------------

**hek293.coverage_profiles:**

1_35358925_35361789.NM_005095.txt
1_20749723_20773610.NM_016287.txt
coverage.clusters.all_circles.pdf
coverage_profiles.all_circles.pdf
