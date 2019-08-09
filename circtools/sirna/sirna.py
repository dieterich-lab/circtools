import circ_module.circ_template

import os
import sys
import subprocess
import pybedtools
pybedtools.set_bedtools_path("/biosw/bedtools/2.27.1/bin")
import tempfile
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
import pandas as pd
from pandas import DataFrame
from IPython.display import HTML

class Sirna(circ_module.circ_template.CircTemplate): 

    def __init__(self, argparse_arguments, program_name, version):

        #create short and long mode
        
        self.cli_params = argparse_arguments
        self.program_name = program_name
        self.version = version

        self.organism = self.cli_params.organism
        self.homo_sapiens_blast_db = "GPIPE/9606/current/rna"
        self.mus_musculus_blast_db = "GPIPE/10090/current/rna"
        self.rattus_norvegicus_blast_db = "GPIPE/10116/current/rna"
        self.other_blast_db = "nt"
        
        self.dcc_file = self.cli_params.dcc_file

        self.fasta_file = self.cli_params.fasta_file
        self.gtf_file = self.cli_params.gtf_file
        
        self.gene_list = self.cli_params.gene_list
        self.id_list = self.cli_params.id_list
        self.input_circRNA = self.cli_params.sequence_file
        
        self.experiment_title = self.cli_params.experiment_title
        self.temp_dir = self.cli_params.global_temp_dir
        self.output_dir = self.cli_params.output_dir
        
        self.target = self.cli_params.target

        # define cache dicts
        self.exon_cache = {}
        self.flanking_exon_cache = {}
        self.siRNA_to_circ_cache = {}
        self.siRNA_data_cache = {}
        self.siRNA_blast_cache = {}
        self.siRNA_findParameter_cache = {}

        self.blast_input_file = ""
        self.blast_xml_tmp = self.temp_dir + "circtools_blast_results.xml"
        self.delNum_blast = 0
        self.hitlist_size = 10
        self.sirnacsv = self.temp_dir + "sirna.csv"

        #define user-defined values
        self.find_parameter = self.cli_params.findParameter
        self.overlap_parameter = self.cli_params.overlapParameter
        self.G_length = self.cli_params.GLength
        self.T_length = self.cli_params.TLength
        self.A_length = self.cli_params.ALength
        self.delNum_length = 0
        
        self.no_blast = self.cli_params.blast
        self.mismatch_tolerance = self.cli_params.mismatchTolerance
        self.mismatch_threshhold = self.cli_params.mismatchThreshhold
        self.seed_mismatch = self.cli_params.seedMismatch
        self.overhang_parameter = self.cli_params.overhangParameter


        # define reference dicts
        # make these static and belong to the class instead of the object?
        self.enthalpyDict = {'AA': -6.6, 'UU': -6.6, 'AU': -5.7, 'UA': -8.1, 'CA': -10.5, 'UG': -10.5,
                       'CU': -7.6, 'AG': -7.6, 'GA': -13.3, 'UC': -13.3, 'GU': -10.2, 'AC': -10.2,
                       'CG': -8.0, 'GC': -14.2, 'GG': -12.2, 'CC': -12.2}
        self.entropyDict = {'AA': -18.4, 'UU': -18.4, 'AU': -15.5, 'UA': -22.6, 'CA': -27.8, 'UG': -27.8,
                       'CU': -19.2, 'AG': -19.2, 'GA': -35.5, 'UC': -35.5, 'GU': -26.2, 'AC': -26.2,
                       'CG': -19.4, 'GC': -34.9, 'GG': -29.7, 'CC': -29.7}

    @staticmethod
    def read_annotation_file(annotation_file, entity="gene", string=False):
        """Reads a GTF file
        Will halt the program if file not accessible
        Returns a BedTool object only containing gene sections
        """

        try:
            file_handle = open(annotation_file)
        except PermissionError:
            message = ("Input file " + str(annotation_file) + " cannot be read, exiting.")
            sys.exit(message)
        else:

            with file_handle:
                line_iterator = iter(file_handle)
                bed_content = ""
                print("Start parsing GTF file")
                for line in line_iterator:
                    # we skip any comment lines
                    if line.startswith("#"):
                        continue

                    # split up the annotation line
                    columns = line.split('\t')

                    if not (columns[2] == entity):
                        continue

                    # we do not want any 0-length intervals -> bedtools segfault
                    if int(columns[4]) - int(columns[3]) == 0:
                        continue

                    # extract chromosome, start, stop, score(0), name and strand
                    # we hopefully have a gene name now and use this one for the entry

                    entry = [
                        columns[0],
                        columns[3],
                        columns[4],
                        "name",
                        str(0),
                        columns[6],
                    ]

                    # concatenate lines to one string
                    bed_content += '\t'.join(entry) + "\n"

            if not bed_content:
                exit(-1)

            if string:
                return bed_content
            else:
                # create a "virtual" BED file
                virtual_bed_file = pybedtools.BedTool(bed_content, from_string=True)
                print("Start merging GTF file")


                # we trust that bedtools >= 2.27 is installed. Otherwise this merge will probably fail
                return virtual_bed_file.sort().merge(s=True,  # strand specific
                                                     c="4,5,6",  # copy columns 5 & 6
                                                     o="distinct,distinct,distinct")  # group

    def extractSequence(self):
        if self.id_list and os.access(self.id_list[0], os.R_OK):
            print("Detected supplied circRNA ID file.")
            with open(self.id_list[0]) as f:
                lines = f.read().splitlines()
            self.id_list = lines

        #let's first check if the temporary directory exists
        if not (os.access(self.temp_dir, os.W_OK)):
            print("Temporary directory %s not writable." % self.temp_dir)
            # exit with -1 error if we can't use it
            exit(-1)

        # let's first check if the output directory exists
        if not (os.access(self.output_dir, os.W_OK)):
            print("Output directory %s not writable." % self.output_dir)
            # exit with -1 error if we can't use it
            exit(-1)

        circ_rna_number = 0
        
        if self.input_circRNA:
            from Bio import SeqIO
            for record in SeqIO.parse(self.input_circRNA, "fasta"):

                # from the FASTA file we cannot tell the coordinates of the circRNA
                name = str(record.id)+"_0_0_"+str(len(record.seq))+"_0"
                
                #data_store.write("\t".join([name, str(record.seq), "", "\n"]))
                seq = str(record.seq)
                bsj = seq[-31:] + seq[0:31]
                self.exon_cache[name] = {1: str(record.seq), 2: "", 3: bsj}
                
        else:
            exons = self.read_annotation_file(self.gtf_file, entity="exon")

            dcc_file = self.dcc_file
            with open(dcc_file) as fp:

                for line in fp:

                    # make sure we remove the header
                    if line.startswith('Chr\t'):
                        continue

                    line = line.rstrip()
                    current_line = line.split('\t')

                    if current_line[3] == "not_annotated":
                        continue

                    if self.gene_list and not self.id_list and current_line[3] not in self.gene_list:
                        continue

                    sep = "_"
                    name = sep.join([current_line[3],
                                     current_line[0],
                                     current_line[1],
                                     current_line[2],
                                     current_line[5]])

                    if self.id_list and not self.gene_list and name not in self.id_list:
                        continue

                    circrna_length = int(current_line[2]) - int(current_line[1])

                    self.flanking_exon_cache[name] = {}

                    sep = "\t"
                    bed_string = sep.join([current_line[0],
                                           current_line[1],
                                           current_line[2],
                                           current_line[3],
                                           str(0),
                                           current_line[5]])

                    virtual_bed_file = pybedtools.BedTool(bed_string, from_string=True)
                    result = exons.intersect(virtual_bed_file, s=True)

                    fasta_bed_line_start = ""
                    fasta_bed_line_stop = ""

                    start = 0
                    stop = 0

                    for result_line in str(result).splitlines():
                        bed_feature = result_line.split('\t')

                        # this is a single-exon circRNA
                        if bed_feature[1] == current_line[1] and bed_feature[2] == current_line[2]:
                            fasta_bed_line_start += result_line + "\n"
                            start = 1
                            stop = 1

                        if bed_feature[1] == current_line[1] and start == 0:
                            fasta_bed_line_start += result_line + "\n"
                            start = 1

                        if bed_feature[2] == current_line[2] and stop == 0:
                            fasta_bed_line_stop += result_line + "\n"
                            stop = 1

                        # these exons are kept for correctly drawing the circRNAs later
                        # not used for primer design
                        if bed_feature[1] > current_line[1] and bed_feature[2] < current_line[2]:
                            self.flanking_exon_cache[name][bed_feature[1] + "_" + bed_feature[2]] = 1

                    virtual_bed_file_start = pybedtools.BedTool(fasta_bed_line_start, from_string=True)
                    virtual_bed_file_stop = pybedtools.BedTool(fasta_bed_line_stop, from_string=True)

                    virtual_bed_file_start = virtual_bed_file_start.sequence(fi=self.fasta_file)
                    virtual_bed_file_stop = virtual_bed_file_stop.sequence(fi=self.fasta_file)

                    if stop == 0 or start == 0:
                        print("Could not identify the exact exon-border of the circRNA.")
                        print("Will continue with non-annotated, manually extracted sequence.")

                        # we have to manually reset the start position

                        fasta_bed_line = "\t".join([current_line[0],
                                                    current_line[1],
                                                    current_line[2],
                                                    current_line[5]])

                        virtual_bed_file_start = pybedtools.BedTool(fasta_bed_line, from_string=True)
                        virtual_bed_file_start = virtual_bed_file_start.sequence(fi=self.fasta_file)
                        virtual_bed_file_stop = ""

                    exon1 = ""
                    exon2 = ""

                    # double check the naming of the exons (for example the "\n" stuff)
                    if virtual_bed_file_start:
                        exon1 = open(virtual_bed_file_start.seqfn).read().split("\n", 1)[1].rstrip()

                    if virtual_bed_file_stop:
                        exon2 = open(virtual_bed_file_stop.seqfn).read().split("\n", 1)[1].rstrip()

                    circ_rna_number += 1
                    print("extracting flanking exons for circRNA #", circ_rna_number, name)

                    if exon2 and not exon1:
                        exon1 = exon2
                        exon2 = ""

                    # sequence from which potential siRNAs can be extracted 
                    bsj = exon1[-31:] + exon2[0:31]

                    self.exon_cache[name] = {1: exon1, 2: exon2, 3: bsj}

        if not self.exon_cache:
            print("Could not find any circRNAs matching your criteria, exiting.")
            exit(-1)

    def reverseComplement(self, circ):
        sequence = self.exon_cache[circ][3]
        reverseComplement = ""
        for x in sequence:
            if x == "A":
                reverseComplement = "U" + reverseComplement
            if x == "T":
                reverseComplement = "A" + reverseComplement
            if x == "G":
                reverseComplement = "C" + reverseComplement
            if x == "C":
                reverseComplement = "G" + reverseComplement
        self.exon_cache[circ][3] = reverseComplement
        
    def DNAtoRNA(self, circ):
        rna = ""
        sequence = self.exon_cache[circ][3]
        for x in sequence:
            if x == "A":
                rna = rna + "A"
            if x == "T":
                rna = rna + "U" 
            if x == "G":
                rna = rna + "G" 
            if x == "C":
                rna = rna + "C"
        self.exon_cache[circ][3] = rna
    
    def reverseComplementRNA(self, sequence):
        reverseComplement = ""
        for x in sequence:
            if x == "A":
                reverseComplement = "U" + reverseComplement
            if x == "U":
                reverseComplement = "A" + reverseComplement
            if x == "G":
                reverseComplement = "C" + reverseComplement
            if x == "C":
                reverseComplement = "G" + reverseComplement
        return reverseComplement
    
    def complementRNA(self, sequence):
        complement = ""
        for x in sequence:
            if x == "A":
                complement = complement + "U"
            if x == "U":
                complement = complement + "A"
            if x == "G":
                complement = complement + "C" 
            if x == "C":
                complement = complement + "G"
        return complement

    def findsiRNAsUiTei(self, circ):

        siRNAList = []
        startPosition = 0
        endPosition = 19
        BSJPOSITION = 29
        sampleString = self.exon_cache[circ][3]
        for x in sampleString:
            if x == "A" or x == "U":
                y = sampleString[(startPosition + 18):endPosition]
                if y == "G" or y == "C":
                    myRNA = sampleString[startPosition:endPosition]

                    if startPosition < (BSJPOSITION-self.overlap_parameter+1) and endPosition > (BSJPOSITION+self.overlap_parameter):
                        if len(myRNA) == 19:
                            #myRNA = self.complementRNA(myRNA)
                            siRNAList.append(myRNA)
                            #Assuming all siRNAs for all circs are unique, no overwriting of the same key should occur (Length of 19 should be sufficient to ensure uniqueness)
                            self.siRNA_findParameter_cache[myRNA] = {0}
            startPosition = startPosition + 1
            endPosition = startPosition + 19
        self.siRNA_to_circ_cache[circ] = {1: siRNAList, 2: 'Ui-Tei'}

    def findsiRNAsReynolds(self, circ):

        siRNAList = []
        startPosition = 0
        endPosition = 19
        BSJPOSITION = 29
        sampleString = self.exon_cache[circ][3]
        for x in sampleString:
            if x == "A" or x == "U":
                myRNA = sampleString[startPosition:endPosition]
                if startPosition < (BSJPOSITION-self.overlap_parameter+1) and endPosition > (BSJPOSITION+self.overlap_parameter):
                    if len(myRNA) == 19:
                        #myRNA = self.complementRNA(myRNA)
                        siRNAList.append(myRNA)
                        self.siRNA_findParameter_cache[myRNA] = {1}
            startPosition = startPosition + 1
            endPosition = startPosition + 19
        self.siRNA_to_circ_cache[circ] = {1: siRNAList, 2: 'Reynolds'}

    def findsiRNAs_multiLength(self, circ):
        #for circ in self.exon_cache:

        siRNAList = []
        startPosition = 0
        endPosition = 19
        BSJPOSITION = 29
        sampleString = self.exon_cache[circ][3]
        for x in sampleString:
            if x == "A" or x == "U":
                y = sampleString[(startPosition + 18):endPosition]
                if y == "G" or y == "C":
                    myRNA = sampleString[startPosition:endPosition]
                    # make more efficient by substringing the exon_storage entry
                    # (have to change BSJPosition location too)
                    # the start/end position limitations were arbitrary, need to figure that out too
                    if startPosition < (BSJPOSITION-self.overlap_parameter+1) and endPosition > (BSJPOSITION+self.overlap_parameter):
                        #myRNA = self.complementRNA(myRNA)
                        siRNAList.append(myRNA)
                        self.siRNA_findParameter_cache[myRNA] = {2}
            startPosition = startPosition + 1
            endPosition = startPosition + 19

        for x in sampleString:
            if x == "A" or x == "U":
                substring = sampleString[endPosition:(endPosition+9)]
                for y in substring:
                    if y == "G" or y == "C":
                        if startPosition < (BSJPOSITION - 3) and endPosition > (BSJPOSITION + 4):
                            myRNA = sampleString[startPosition:(endPosition+1)]
                            
                            siRNAList.append(myRNA)
                        endPosition = endPosition + 1
            startPosition = startPosition + 1
            endPosition = startPosition + 19
        self.siRNA_to_circ_cache[circ] = {1: siRNAList, 2: 'Ui-Tei'}

    def findsiRNAs(self, circ):
        if self.find_parameter == 0:
            self.findsiRNAsUiTei(circ)
        if self.find_parameter == 1:
            self.findsiRNAsReynolds(circ)
        if self.find_parameter == 2:
            self.findsiRNAs_multiLength(circ)

    def deleteRepeats(self, circ):
        siRNAList = self.siRNA_to_circ_cache[circ][1]
        G_Parameter = ""
        for i in range(0, self.G_length + 1, 1):
            G_Parameter += "G"  

        for i in range(len(siRNAList) - 1, -1, -1):
            x = siRNAList[i].find(G_Parameter)
            if x != -1:
                del siRNAList[i]
                self.delNum_length += 1

        T_Parameter = ""
        for i in range(0, self.T_length + 1, 1):
            T_Parameter += "T"
        for i in range(len(siRNAList) - 1, -1, -1):
            x = siRNAList[i].find(T_Parameter)
            if x != -1:
                del siRNAList[i]  
                self.delNum_length += 1

        A_Parameter = ""
        for i in range(0, self.A_length + 1, 1):
            A_Parameter += "G"
        for i in range(len(siRNAList) - 1, -1, -1):
            x = siRNAList[i].find(A_Parameter)
            if x != -1:
                del siRNAList[i]
                self.delNum_length += 1

        self.siRNA_to_circ_cache[circ][1] = siRNAList

    def calculateNNEnthalpy(self, sequence):
        sequencePairs = []
        for n in range(0, len(sequence) - 1, 1):
            pair = sequence[n:n+2]
            sequencePairs.append(pair)
        enthalpySequence = 0
        for p in sequencePairs:
            enthalpySequence += self.enthalpyDict[p]
        return enthalpySequence

    def calculateNNEntropy(self, sequence):
        sequencePairs = []
        for n in range(0, len(sequence) - 1, 1):
            pair = sequence[n:n+2]
            sequencePairs.append(pair)
        entropySequence = 0
        for p in sequencePairs:
            entropySequence += self.entropyDict[p]
        return entropySequence 

    def calculateTm(self, sequence):
        enthalpy = self.calculateNNEnthalpy(sequence)
        entropy = self.calculateNNEntropy(sequence)
        Tm = (1000*enthalpy)/(entropy-31.8555) - 289.75
        return Tm

    def calculateSeedStability(self, sequence):
        seedTm = self.calculateTm(sequence[1:8])
        return seedTm

    #@staticmethod
    def calculateGCContent(self, rna_string):
        myString = rna_string
        gcNum = 0;
        for x in myString:
            if x == "G" or x == "C":
                gcNum = gcNum + 1
        stringLength = len(myString)
        gcContent = gcNum / stringLength
        return gcContent
    
    #@staticmethod
    def calculateGCStretch(self, rna_string):
        myString = rna_string
        stretch = 0
        
        i = 0
        while i < len(myString):
            newStretch = 0
            if myString[i:i+1] == "G" or myString[i:i+1] == "C":
                newStretch += 1
                while myString[i+1:i+2] == "G" or myString[i+1:i+2] == "C":
                    newStretch += 1
                    i += 1
            if newStretch > stretch:
                stretch = newStretch
            i += 1
        return stretch

    #@staticmethod
    def calculateAUContent(self, rna_string):
        myString = rna_string
        gcNum = 0;
        for x in myString:
            if x == "A" or x == "U":
                gcNum = gcNum + 1
        stringLength = len(myString)
        gcContent = gcNum / stringLength
        return gcContent

    #@staticmethod
    def calculateScoreUiTei(self, rna_string):
        rnaString = rna_string
        score = 0
        #double check these values
        thirdString = rnaString[0:6]
        thirdAUContent = self.calculateAUContent(thirdString)
        if thirdAUContent > .5:
            score = score + 5
        if thirdAUContent > .6:
            score = score + 10
        if thirdAUContent > .7:
            score = score + 5

        #double check these values
        twoThirdString = rnaString[7:14]
        twoThirdGCContent = self.calculateGCContent(twoThirdString)
        if twoThirdGCContent > .5:
            score = score + 5
        if twoThirdGCContent > .6:
            score = score + 10
        if twoThirdGCContent > .7:
            score = score + 5
        GCContent = self.calculateGCContent(rnaString)
        if GCContent > .32 and GCContent < .50:
            score = score + 25
        positionTen = rnaString[9:10]
        if positionTen == "A" or positionTen == "U":
            score = score + 10
        GCStretch = self.calculateGCStretch(rnaString)
        if GCStretch >= 10:
            score = score - 25
        return score

    #@staticmethod
    def calculateScoreReynolds(self, rna_string):
        score = 0
        rnaString = rna_string
        GCContent = self.calculateGCContent(rnaString)
        if GCContent > .32 and GCContent < .50:
            score = score + 25
        initialAUContent = self.calculateAUContent(rnaString[0:5])
        if initialAUContent >= .60:
            score = score + 25
        if rnaString[0:1] == "A":
            score = score + 5
        positionTen = rnaString[9:10]
        if positionTen == "U":
            score = score + 5
        positionSeventeen = rnaString[16:17]
        if positionSeventeen == "A":
            score = score + 5
        positionSeven = rnaString[6:7]
        if positionSeven == "G":
            score = score - 5
        return score

    #@staticmethod
    def calculateScore(self, rna_string, x):
        score = 0
        if x == 0:
            score = self.calculateScoreUiTei(rna_string)
        if x == 1:
            score = self.calculateScoreReynolds(rna_string)
        if x == 2:
            score = self.calculateScoreUiTei(rna_string)
        return score

    def scoreSiRNAs(self, circ):
        siRNAList = self.siRNA_to_circ_cache[circ][1]
        scoreList = []
        for a in siRNAList:
            scoreFindParameter = self.siRNA_findParameter_cache[a]
            print("circ: " + circ + "findParameter: " + str(scoreFindParameter))
            tempScore = self.calculateScore(a, scoreFindParameter)
            scoreList.append(tempScore)
        scores = {'siRNA': siRNAList, 'Silencing Score': scoreList}
        scoresDF = DataFrame(scores, columns=['siRNA', 'Silencing Score'])
        scoresDF = scoresDF.sort_values(['Silencing Score'], ascending=0)
        orderedsiRNAList = scoresDF['siRNA'].tolist()
        self.siRNA_to_circ_cache[circ][1] = scoresDF
            
    def calculateSeedDuplex(self, circ):
        siRNAList = self.siRNA_to_circ_cache[circ][1]['siRNA'].tolist()
        Tmlist = []
        for i in range(len(siRNAList) - 1, -1, -1):
            Tm = self.calculateSeedStability(siRNAList[i])
            mysiRNA = siRNAList[i]
            #Code to filter out siRNAs based on seed-Dupled Tm:
            #if (Tm > self.seedDuplex_threshhold):
            #    self.siRNA_to_circ_cache[circ] = self.siRNA_to_circ_cache[circ][self.siRNA_to_circ_cache[circ].siRNA != mysiRNA]
            #else:
            Tmlist.insert(0, Tm)
        self.siRNA_to_circ_cache[circ][1]['Seed-Duplex stability (Tm C)'] = Tmlist
        
    def calculateThermodynamicStability(self, circ):
        siRNAList = self.siRNA_to_circ_cache[circ][1]['siRNA'].tolist()
        Thermolist = []
        for i in range(len(siRNAList) - 1, -1, -1):
            fivePrime = siRNAList[i][0:4]
            fiveTherm = self.calculateTm(fivePrime)
            #fiveTherm = self.calculateNNEnthalpy(fivePrime)
            threePrime = siRNAList[i][-4:]
            threeTherm = self.calculateTm(threePrime)
            #threeTherm = self.calculateNNEnthalpy(threePrime)
            therm = fiveTherm - threeTherm
            Thermolist.insert(0, therm)
        self.siRNA_to_circ_cache[circ][1]['Thermodynamic Stability'] = Thermolist

    def call_blast(self, input_file, organism):
        blast_db = "nt"
        if organism == "mm":
            blast_db = self.mus_musculus_blast_db
        elif organism == "hs":
            blast_db = self.homo_sapiens_blast_db
        elif organism == "rn":
            blast_db = self.rattus_norvegicus_blast_db

        print("Running blast...")
        # double check the other arguments in the qblast call
        return_handle2 = NCBIWWW.qblast("blastn",
                                       blast_db,
                                       input_file,
                                       hitlist_size=self.hitlist_size,
                                       expect=1000,
                                       word_size=7,
                                       gapcosts="5 2"
                                       )
        return return_handle2
    
    def createBlastInputFile(self, circ):
        siRNAList = self.siRNA_to_circ_cache[circ][1]['siRNA'].tolist()
        for i in range(len(siRNAList) - 1, -1, -1):
            mysiRNA = siRNAList[i]
            self.blast_input_file += "\n>" + mysiRNA + "$" + circ + "\n" + mysiRNA

    def runBlast(self):

        #for circ in self.siRNA_to_circ_cache:
        #siRNAList = self.siRNA_to_circ_cache[circ]['siRNA'].tolist()
        #siRNANumber = 0
        #for i in range(len(siRNAList) - 1, -1, -1):
        #    mysiRNA = siRNAList[i]
        #    siRNANumber += 1
        #    print("Sending siRNA #" + str(siRNANumber) + " for " + circ + " to blast...")
        #    return_handle = self.call_blast(mysiRNA, self.organism)
        
        
        result_handle = self.call_blast(self.blast_input_file, self.organism)
        with open(self.blast_xml_tmp, "w") as out_handle:
                out_handle.write(result_handle.read())

        result_handle.close()
        result_handle = open(self.blast_xml_tmp)
        blast_records = NCBIXML.parse(result_handle)

        for blast_record in blast_records:

            queryTitle = blast_record.query.split('$')
            mysiRNA = queryTitle[0]
            circ = queryTitle[1]
            blastCount = 0
            # E-value threshhold not used currently, correct threshhold value based on length of query and size of database?
            E_VALUE_THRESH = 10
            #low complexity filter?
            #for alignment in blast_record.alignments:
            #    for hsp in alignment.hsps:
            #        if hsp.expect < E_VALUE_THRESH:
            #            blastCount = blastCount + 1

            mismatchCount = 0
            blastStorage = []
            alignStorage = []
            blastHits = 0
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    blastHits += 1
                    mismatches = len(mysiRNA) - hsp.identities
                    ##count indels too
                    seedMismatches = 0
                    #hsp.match that's empty - potential bug
                    for match in hsp.match[0:9]:
                        if match != "|":
                            seedMismatches += 1
                    
                    if self.seed_mismatch == True:
                        if seedMismatches <= self.mismatch_tolerance:
                            mismatchCount += 1
                    
                    if self.seed_mismatch == False:
                        if mismatches <= self.mismatch_tolerance:
                            mismatchCount += 1
                            
                    blastString = ""
                    alignString = ""
                    blastString += '****Alignment****\n'
                    blastString += 'sequence:' + alignment.title + '\n'
                    blastString += 'length:' + str(alignment.length) + '\n'
                    blastString += 'e value:' + str(hsp.expect) + '\n'
                    alignString += hsp.query[0:75] + '\n'
                    alignString += hsp.match[0:75] + '\n'
                    alignString += hsp.sbjct[0:75] + '\n'
                    blastStorage.append(blastString)
                    alignStorage.append(alignString)
                    ##show position of mismatches
            blasts = {'Alignment':alignStorage, 'Blast Info':blastStorage}
            blastDf = DataFrame(blasts, columns=['Alignment', 'Blast Info'])
            #blastDf.to_html('blastTest' + mysiRNA + '.html', index = False)
            
            #add weighting of mismatches - mismatches closer to 5' end are more significant
            if mismatchCount >= self.mismatch_threshhold:
                #siRNA_to_circ_cache[circ] = siRNA_to_circ_cache[circ].set_index("siRNA")
                #siRNA_t o_circ_cache[circ].drop([mysiRNA], axis = 0)
                self.siRNA_to_circ_cache[circ][1] = self.siRNA_to_circ_cache[circ][1][self.siRNA_to_circ_cache[circ][1].siRNA != mysiRNA]
                self.delNum_blast += 1
            else:
                blastHits2 = blastHits - self.hitlist_size
                self.siRNA_blast_cache[mysiRNA] = blastHits2
            

            # this statement might need to be changed since if the blast_records is too large,
            # converting to a list can break the program

            #if blastCount > 0:
            #    del siRNAList[i]
        #siRNA_to_circ_cache[circ] = siRNAList

    #create output method that shows which rules were satisfied by each siRNA etc. 
    def showOutput(self, circ):
        #for circ in self.siRNA_to_circ_cache:
        scoresDF = self.siRNA_to_circ_cache[circ][1]
        if scoresDF.empty:
            print("Could not find any siRNAs targeting " + circ)
            print("Try running in a different mode or relaxing the input parameters")
        else:
            print("Showing siRNAs and their scores for " + circ)
    
    #writes output to an HTML
    #add method for converting to a csv or pdf file
    def createOutput(self, circ, dfList):
        
        #empty = False
        #for circ in self.siRNA_to_circ_cache:
        scoresDf = self.siRNA_to_circ_cache[circ][1]
        totalRows = len(scoresDf.index)
        if self.no_blast == False:
            blastList = []
            for x in scoresDf['siRNA']:
                blastNum = self.siRNA_blast_cache[x]
                blastList.append(blastNum)
        if self.no_blast == True:
            blastList = []
            for i in range (0, totalRows, 1):
                blastList.append("N/A")
        newStrandList = []
        if self.overhang_parameter == 0:
            for x in scoresDf['siRNA']:
                newString = x+"UU"
                tempString = x+"UU"
                sense = self.reverseComplementRNA(x)
                sense = sense+"UU"
                newString  = "Guide: " + newString + "<br>" + "Passenger: " + sense
                newStrandList.append(newString)
                #scoresDf = scoresDf.replace(x, newString)
                #self.siRNA_to_circ_cache[circ][1] = self.siRNA_to_circ_cache[circ][1].replace(x, tempString)
        elif self.overhang_parameter == 1:
            for x in scoresDf['siRNA']:
                newString = x+"TT"
                tempString =x+"TT"
                sense = self.reverseComplementRNA(x)
                sense = sense+"TT"
                newString  = "Guide: " + newString + "<br>" + "Passenger: " + sense
                newStrandList.append(newString)
                #scoresDf = scoresDf.replace(x, newString)
                #self.siRNA_to_circ_cache[circ][1] = self.siRNA_to_circ_cache[circ][1].replace(x, tempString)
        #if scoresDf.empty:
            #print("Could not find any siRNAs targeting " + circ)
            #empty = True
            #continuei
        else:
            for x in scoresDf['siRNA']:
                sense = self.reverseComplementRNA(x)
                newString = "Guide: " + x + "<br>" + "Passenger: " + sense
                newStrandList.append(newString)
        geneParams = circ.split('_')
        
        strandList = []
        stopList = []
        startList = []
        chrList = []
        annotList = []
        ruleList = []
        
        for i in range(0, totalRows, 1):
            strandList.append(geneParams[4])
            stopList.append(geneParams[3])
            startList.append(geneParams[2])
            chrList.append(geneParams[1])
            annotList.append(geneParams[0])
            ruleList.append(self.siRNA_to_circ_cache[circ][2])
        scoresDf.insert(0, 'Strand', strandList)
        scoresDf.insert(0, 'Stop', stopList)
        scoresDf.insert(0, 'Start', startList)
        scoresDf.insert(0, 'Chr', chrList)
        scoresDf.insert(0, 'Annot', annotList)
        scoresDf.insert(7, 'Rule', ruleList)
        scoresDf.insert(8, 'Blast', blastList)
        scoresDf.insert(5, 'newsiRNA', newStrandList)
        #scoresDf.rename(columns={'siRNA':'Anti-sense Guide siRNA'}, inplace=True)
        dfList.append(scoresDf)
       # if empty:
        #    print("Try running in a different mode or relaxing the input parameters")

        
    def writeOutput(self, totalDf):
        # let's first check if the temporary directory exists
       # if not (os.access(self.temp_dir, os.W_OK)):
        #    print("Temporary directory %s not writable." % self.temp_dir)
            # exit with -1 error if we can't use it
        #    exit(-1)
            
        testCsv = totalDf.to_csv()
        with open(self.sirnacsv, "w") as siRNA_file:
            siRNA_file.write(testCsv)
        blast_r = "True"
        if self.no_blast:
            blast_r = "False"
        subprocess.check_output(['Rscript', "../scripts/circtools_sirna_formatter.R", self.sirnacsv, blast_r, self.output_dir, self.experiment_title])
    
    def drawsiRNA(self, circ):
        if circ in self.exon_cache:
            circ_info = circ.split('_')
            circrna_length = int(circ_info[3]) - int(circ_info[2])
            exon1_length = len(self.exon_cache[circ][1])
            
            exon2_length = len(self.exon_cache[circ][2])
            exon2_colour = "#ffac68"
            
            if exon2_length == 0:
                exon1_length = int(len(self.exon_cache[circular_rna_id][1])/2)+1
                exon2_length = int(len(self.exon_cache[circular_rna_id][1])/2)
                exon2_colour = "#ff6877"
            siRNAList = self.siRNA_to_circ_cache[circ][1]['siRNA']
            for siRNA in siRNAList:
                pos = self.exon_cache[circ][3].find(siRNA)
                siRNA_start = exon2_length-30+ pos 
                
                siRNA_length = len(siRNA)
                siRNA_end = siRNA_start + siRNA_length
               
                
                
                gdd = GenomeDiagram.Diagram('circRNA siRNA diagram')
                gdt_features = gdd.new_track(1, greytrack=True, name="", )
                gds_features = gdt_features.new_set()
                
                feature = SeqFeature(FeatureLocation(exon2_length, exon1_length+exon2_length), strand=-1)
                gds_features.add_feature(feature, name="Exon 1", label=False, color="#ff6877", label_size=22)
                
                feature = SeqFeature(FeatureLocation(0, exon2_length), strand=-1)
                gds_features.add_feature(feature, name="Exon 2", label=False, color=exon2_colour, label_size=22)

                #siRNA
                feature = SeqFeature(FeatureLocation(siRNA_start, siRNA_end), strand=+1)
                gds_features.add_feature(feature, name="siRNA_2", label=False, color="#32CD32", label_size=22)
                
                #more labels
                feature = SeqFeature(FeatureLocation(round(exon2_length/2), round(exon2_length/2)+1), strand =+1)
                gds_features.add_feature(feature, name="exon 2", color="#F5F5F5",label=True, label_size=20)
                feature = SeqFeature(FeatureLocation(round(exon1_length/2)+exon2_length, round(exon1_length/2)
                                                     +exon2_length +1), strand =+1)
                gds_features.add_feature(feature, name="exon 1", color="#F5F5F5",label=True, label_size=20)
                
                feature = SeqFeature(FeatureLocation(siRNA_start, siRNA_start+1), strand =+1)
                gds_features.add_feature(feature, name="siRNA", color="#32CD32",label=True, label_size=22)
                
                feature = SeqFeature(FeatureLocation(exon2_length, exon2_length+1))
                gds_features.add_feature(feature, name="BSJ", label=True, color="white", label_size=22)
                
                gdd.draw(format='linear', pagesize=(600, 600), track_size=0.3, tracklines=0, x=0.00, 
                         y=0.00, start=0, end=exon1_length+exon2_length)
                gdd.write(self.output_dir+circ+siRNA+".svg","svg")
                #if self.overhang_parameter == 0:
                #    gdd.write(self.output_dir + circ + siRNA + "UU"+".svg", "svg")
                #elif self.overhang_parameter == 1:
                #    gdd.write(self.output_dir + circ + siRNA + "TT"+".svg", "svg")
                #else:
                #    gdd.write(self.ouput_dir + circ + siRNA +".svg", "svg")
            
            
    def run_module(self):
        #get the circRNA sequence
        self.extractSequence()
        
        #convert the circRNA sequence into a sequence from which potential siRNAs can be extracted
        if self.target == "anti-sense":
            for circ in self.exon_cache:
                self.reverseComplement(circ)
        if self.target == "sense":
            for circ in self.exon_cache:
                self.DNAtoRNA(circ)
        print("Finding siRNAs...")
        
        #extract potential siRNA sequences
        for circ in self.exon_cache:
            self.findsiRNAs(circ)
            if not self.siRNA_to_circ_cache[circ][1] and self.find_parameter == 0:
                print("Could not find any siRNAs for " + circ + " using the Ui-Tei rule")
                print("Finding siRNAs using Reynolds rule")
                self.findsiRNAsReynolds(circ)
                
        #delete siRNAs with repeats
        for circ in self.siRNA_to_circ_cache:
            self.deleteRepeats(circ)
        print(str(self.delNum_length) + " siRNAs were eliminated due to repeats")
        
        #score the silencing potential of each siRNA
        print("Scoring siRNAs...")
        for circ in self.siRNA_to_circ_cache:
            self.scoreSiRNAs(circ)
        for circ in self.siRNA_to_circ_cache:
            self.calculateSeedDuplex(circ)
        for circ in self.siRNA_to_circ_cache:
            self.calculateThermodynamicStability(circ)
        for circ in self.siRNA_to_circ_cache:
            self.createBlastInputFile(circ)
        if self.no_blast == False:
            self.runBlast()
            print(str(self.delNum_blast) + " siRNAs were eliminated in the blast step")
        empty = False
        dfList = []
        for circ in self.siRNA_to_circ_cache:
            scoresDf = self.siRNA_to_circ_cache[circ][1]
            if scoresDf.empty:
                empty = True
                print("Could not find any siRNAs targeting " + circ)
            else:   
                self.createOutput(circ, dfList)
                self.drawsiRNA(circ)
        
        if empty:
            print("Try running in a different mode or relaxing the input parameters")
        
        if dfList:
            print("Writing siRNAs and their scores to an html file")
            totalDf = pd.concat(dfList)
            self.writeOutput(totalDf)
            
            #fix threshhold value of 30
            print("A score above 50 represents an effective siRNA")

#method for outputing good features (silencing score, etc.) of an input siRNA
 
