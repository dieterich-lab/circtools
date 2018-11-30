#! /usr/bin/env python3

# Copyright (C) 2017 Tobias Jakobi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either self.version 3 of the License, or
# (at your option) any later self.version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import circ_module.circ_template

import os
import sys

import pybedtools
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram


class Primex(circ_module.circ_template.CircTemplate):
    def __init__(self, argparse_arguments, program_name, version):

        # get the user supplied options
        self.cli_params = argparse_arguments
        self.program_name = program_name
        self.version = version
        self.command = 'Rscript'
        self.temp_dir = self.cli_params.global_temp_dir
        self.gtf_file = self.cli_params.gtf_file
        self.fasta_file = self.cli_params.fasta_file
        self.dcc_file = self.cli_params.dcc_file
        self.output_dir = self.cli_params.output_dir
        self.organism = self.cli_params.organism
        self.gene_list = self.cli_params.gene_list
        self.id_list = self.cli_params.id_list
        self.product_range = self.cli_params.product_size
        self.junction = self.cli_params.junction
        self.no_blast = self.cli_params.blast
        self.experiment_title = self.cli_params.experiment_title
        self.input_circRNA = self.cli_params.sequence_file

        if self.id_list and self.gene_list:
            print("Please specify either host genes via -G or circRNA IDs via -i.")
            sys.exit(-1)

        self.homo_sapiens_blast_db = "GPIPE/9606/current/rna"
        self.mus_musculus_blast_db = "GPIPE/10090/current/rna"
        self.rattus_norvegicus_blast_db = "GPIPE/10116/current/rna"

        self.other_blast_db = "nt"

    def module_name(self):
        """"Return a string representing the name of the module."""
        return self.program_name

    # Register an handler for the timeout
    def handler(self, signum, frame):
        raise Exception("Maximum execution time for remote BLAST reached. Please try again later.")

    def call_blast(self, input_file, organism):

        blast_db = "nt"
        if organism == "mm":
            blast_db = self.mus_musculus_blast_db
        elif organism == "hs":
            blast_db = self.homo_sapiens_blast_db
        elif organism == "rn":
            blast_db = self.rattus_norvegicus_blast_db

        return_handle = NCBIWWW.qblast("blastn",
                                       blast_db,
                                       input_file,
                                       hitlist_size=10,
                                       expect=1000,
                                       word_size=7,
                                       gapcosts="5 2"
                                       )
        return return_handle

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

    def run_module(self):

        if self.id_list and os.access(self.id_list[0], os.R_OK):
            print("Detected supplied circRNA ID file.")
            with open(self.id_list[0]) as f:
                lines = f.read().splitlines()
            self.id_list = lines

        # let's first check if the temporary directory exists
        if not (os.access(self.temp_dir, os.W_OK)):
            print("Temporary directory %s not writable." % self.temp_dir)
            # exit with -1 error if we can't use it
            exit(-1)

        # let's first check if the temporary directory exists
        if not (os.access(self.output_dir, os.W_OK)):
            print("Output directory %s not writable." % self.output_dir)
            # exit with -1 error if we can't use it
            exit(-1)

        circ_rna_number = 0

        # define temporary files
        exon_storage_tmp = self.temp_dir + "circtools_flanking_exons.tmp"
        blast_storage_tmp = self.temp_dir + "circtools_blast_results.tmp"
        blast_xml_tmp = self.temp_dir + "circtools_blast_results.xml"

        output_html_file = self.output_dir + self.experiment_title.replace(" ", "_") + ".html"

        # erase old contents
        open(exon_storage_tmp, 'w').close()

        # define cache dicts
        exon_cache = {}
        flanking_exon_cache = {}
        primer_to_circ_cache = {}

        if self.input_circRNA:
            from Bio import SeqIO
            with open(exon_storage_tmp, 'a') as data_store:
                for record in SeqIO.parse(self.input_circRNA, "fasta"):

                    # from the FASTA file we cannot tell the coordinates of the circRNA
                    name = str(record.id)+"_0_0_"+str(len(record.seq))+"_0"

                    data_store.write("\t".join([name, str(record.seq), "", "\n"]))
                    exon_cache[name] = {1: str(record.seq), 2: ""}

        else:
            exons = self.read_annotation_file(self.gtf_file, entity="exon")

            with open(self.dcc_file) as fp:

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

                    flanking_exon_cache[name] = {}

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
                            flanking_exon_cache[name][bed_feature[1] + "_" + bed_feature[2]] = 1

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

                    if virtual_bed_file_start:
                        exon1 = open(virtual_bed_file_start.seqfn).read().split("\n", 1)[1].rstrip()

                    if virtual_bed_file_stop:
                        exon2 = open(virtual_bed_file_stop.seqfn).read().split("\n", 1)[1].rstrip()

                    circ_rna_number += 1
                    print("extracting flanking exons for circRNA #", circ_rna_number, name, end="\n", flush=True)

                    if exon2 and not exon1:
                        exon1 = exon2
                        exon2 = ""

                    exon_cache[name] = {1: exon1, 2: exon2}

                    with open(exon_storage_tmp, 'a') as data_store:
                        data_store.write("\t".join([name, exon1, exon2, "\n"]))

        if not exon_cache:
            print("Could not find any circRNAs matching your criteria, exiting.")
            exit(-1)

        # need to define path top R wrapper
        primer_script = 'circtools_primex_wrapper.R'

        # ------------------------------------ run script and check output -----------------------

        script_result = os.popen(primer_script + " " +
                                 exon_storage_tmp + " " +
                                 str(self.product_range[0]) + "," + str(self.product_range[1]) + " " +
                                 self.junction).read()

        # this is the first time we look through the input file
        # we collect the primer sequences and unify everything in one blast query

        blast_object_cache = {}
        blast_result_cache = {}

        blast_input_file = ""

        if circ_rna_number < 50:

            for line in script_result.splitlines():
                entry = line.split('\t')
                circular_rna_id = entry[0].split('_')

                if entry[1] == "NA":
                    continue

                # only blast 1
                elif entry[2] in blast_object_cache and not entry[1] in blast_object_cache:
                    blast_input_file += "\n>" + entry[1] + "\n" + entry[1]
                    blast_object_cache[entry[1]] = 1
                    primer_to_circ_cache[entry[1]] = circular_rna_id[0]

                # only blast 2
                elif entry[1] in blast_object_cache and not entry[2] in blast_object_cache:
                    blast_input_file += "\n>" + entry[2] + "\n" + entry[2]
                    blast_object_cache[entry[2]] = 1
                    primer_to_circ_cache[entry[2]] = circular_rna_id[0]

                # seen both already, skip
                elif entry[1] in blast_object_cache and entry[2] in blast_object_cache:
                    continue

                # nothing seen yet, blast both
                else:
                    blast_input_file += "\n>" + entry[1] + "\n" + entry[1] + "\n>" + entry[2] + "\n" + entry[2]
                    blast_object_cache[entry[1]] = 1
                    blast_object_cache[entry[2]] = 1
                    primer_to_circ_cache[entry[1]] = circular_rna_id[0]
                    primer_to_circ_cache[entry[2]] = circular_rna_id[0]
        else:
            print("Too many circRNAs selected, skipping BLAST step.")

        if self.no_blast:
            print("User disabled BLAST search, skipping.")

        run_blast = 0

        # check if we have to blast
        if not self.no_blast and blast_input_file:

            try:
                print("Sending " + str(len(blast_object_cache)) + " primers to BLAST")
                print("This may take a few minutes, please be patient.")
                result_handle = self.call_blast(blast_input_file, self.organism)
                run_blast = 1
            except Exception as exc:
                print(exc)
                exit(-1)

            with open(blast_xml_tmp, "w") as out_handle:
                out_handle.write(result_handle.read())

            result_handle.close()
            result_handle = open(blast_xml_tmp)

            blast_records = NCBIXML.parse(result_handle)

            for blast_record in blast_records:

                if blast_record.query not in blast_result_cache:
                    blast_result_cache[blast_record.query] = []

                for description in blast_record.descriptions:

                    # filter out the host gene we're in now
                    # also filter out all "PREDICTED" stuff
                    if description.title.find(primer_to_circ_cache[blast_record.query]) == -1 and\
                            description.title.find("PREDICTED") == -1:
                        blast_result_cache[blast_record.query].append(description.title)

        # if we encounter NAs nothing has been blasted, we manually set the values now

        blast_result_cache["NA"] = ["Not blasted, no primer pair found"]

        primex_data_with_blast_results = ""

        for line in script_result.splitlines():
            entry = line.split('\t')

            # split up the identifier for final plotting
            line = line.replace("_", "\t")

            if run_blast == 1:
                left_result = "No hits"
                right_result = "No hits"
            else:
                left_result = "Not blasted, no primer pair found"
                right_result = left_result

            if entry[1] in blast_result_cache:
                left_result = ";".join(blast_result_cache[entry[1]])

            if entry[2] in blast_result_cache:
                right_result = ";".join(blast_result_cache[entry[2]])

            # update line
            primex_data_with_blast_results += line + "\t" + left_result + "\t" + right_result + "\n"

        with open(blast_storage_tmp, 'w') as data_store:
            data_store.write(primex_data_with_blast_results)

        # need to define path top R wrapper
        primer_script = 'circtools_primex_formatter.R'

        # ------------------------------------ run script and check output -----------------------

        primex_data_formatted = os.popen(primer_script + " " +
                                         blast_storage_tmp + " "
                                         + "\"" + self.experiment_title + "\""
                                         ).read()

        with open(output_html_file, 'w') as data_store:
            data_store.write(primex_data_formatted)

        print("Writing results to "+output_html_file)

        # here we create the circular graphics for primer visualisation
        for line in primex_data_with_blast_results.splitlines():
            entry = line.split('\t')

            # no primers, no graphics
            if entry[6] == "NA":
                continue

            circular_rna_id = "_".join([entry[0],
                                        entry[1],
                                        entry[2],
                                        entry[3],
                                        entry[4]])

            if circular_rna_id in exon_cache:

                circular_rna_id_isoform = circular_rna_id + "_" + entry[5]

                circrna_length = int(entry[3]) - int(entry[2])

                exon1_length = len(exon_cache[circular_rna_id][1])
                exon2_length = len(exon_cache[circular_rna_id][2])

                exon2_colour = "#ffac68"

                if exon2_length == 0:
                    exon1_length = int(len(exon_cache[circular_rna_id][1])/2)+1
                    exon2_length = int(len(exon_cache[circular_rna_id][1])/2)
                    exon2_colour = "#ff6877"

                forward_primer_start = int(entry[8].split(',')[0]) + circrna_length - exon2_length
                forward_primer_length = int(entry[8].split(',')[1])

                reverse_primer_start = int(entry[9].split(',')[0]) - exon2_length
                reverse_primer_length = int(entry[9].split(',')[1])

                product_size = entry[14]

                gdd = GenomeDiagram.Diagram('circRNA primer diagram')
                gdt_features = gdd.new_track(1, greytrack=True, name="", )
                gds_features = gdt_features.new_set()

                feature = SeqFeature(FeatureLocation(0, exon1_length), strand=+1)
                gds_features.add_feature(feature, name="Exon 1", label=False, color="#ff6877", label_size=22)

                feature = SeqFeature(FeatureLocation(circrna_length - exon2_length, circrna_length), strand=+1)
                gds_features.add_feature(feature, name="Exon 2", label=False, color=exon2_colour, label_size=22)

                feature = SeqFeature(FeatureLocation(forward_primer_start, circrna_length), strand=-1)
                gds_features.add_feature(feature, name="Product", label=False, color="#6881ff")

                feature = SeqFeature(FeatureLocation(0, reverse_primer_start), strand=-1)
                gds_features.add_feature(feature, name="Product: " + product_size + "bp", label=False, color="#6881ff",
                                         label_size=22, label_position="middle")

                if self.junction == "f":

                    feature = SeqFeature(FeatureLocation(reverse_primer_start - reverse_primer_length, reverse_primer_start),
                                         strand=-1)
                    gds_features.add_feature(feature, name="Reverse", label=False, sigil="BIGARROW", color="#75ff68",
                                             arrowshaft_height=0.3, arrowhead_length=0.1, label_size=22)

                    # the primer spans the BSJ, therefore we have to draw it in two pieces:
                    # piece 1: primer start to circRNA end
                    # piece 2: remaining primer portion beginning from 0

                    # piece 1:
                    feature = SeqFeature(FeatureLocation(forward_primer_start, circrna_length))
                    gds_features.add_feature(feature, name="Forward", label=False, sigil="BIGARROW", color="#75ff68",
                                             arrowshaft_height=0.3, arrowhead_length=0.1, label_size=22)

                    # piece 2:
                    feature = SeqFeature(
                        FeatureLocation(0, forward_primer_length - (circrna_length - forward_primer_start)))
                    gds_features.add_feature(feature, name="Forward", label=False, sigil="BIGARROW", color="#75ff68",
                                             arrowshaft_height=0.3, arrowhead_length=0.1, label_size=22)
                elif self.junction == "r":
                    # the primer spans the BSJ, therefore we have to draw it in two pieces:
                    # piece 1: primer start of circRNA to circRNA end
                    # piece 2: remaining primer portion beginning from 0

                    # piece 1:
                    feature = SeqFeature(FeatureLocation(circrna_length - reverse_primer_start, circrna_length), strand=-1)
                    gds_features.add_feature(feature, name="Reverse", label=False, sigil="BIGARROW", color="#75ff68",
                                             arrowshaft_height=0.3, arrowhead_length=0.1, label_size=22)

                    # piece 2:
                    feature = SeqFeature(
                        FeatureLocation(0, reverse_primer_start), strand=-1)
                    gds_features.add_feature(feature, name="Reverse", label=False, sigil="BIGARROW", color="#75ff68",
                                             arrowshaft_height=0.3, arrowhead_length=0.1, label_size=22)

                    feature = SeqFeature(
                        FeatureLocation(forward_primer_start, forward_primer_start + forward_primer_length))
                    gds_features.add_feature(feature, name="Forward", label=False, sigil="BIGARROW", color="#75ff68",
                                             arrowshaft_height=0.3, arrowhead_length=0.1, label_size=22)
                else:
                    feature = SeqFeature(
                        FeatureLocation(reverse_primer_start - reverse_primer_length, reverse_primer_start),
                        strand=-1)
                    gds_features.add_feature(feature, name="Reverse", label=False, sigil="BIGARROW", color="#75ff68",
                                             arrowshaft_height=0.3, arrowhead_length=0.1, label_size=22)

                    feature = SeqFeature(
                        FeatureLocation(forward_primer_start, forward_primer_start + forward_primer_length))
                    gds_features.add_feature(feature, name="Forward", label=False, sigil="BIGARROW", color="#75ff68",
                                             arrowshaft_height=0.3, arrowhead_length=0.1, label_size=22)

                feature = SeqFeature(FeatureLocation(0, 1))
                gds_features.add_feature(feature, name="BSJ", label=True, color="white", label_size=22)

                if circular_rna_id in flanking_exon_cache:
                    for exon in flanking_exon_cache[circular_rna_id]:
                        exon_start, exon_stop = exon.split('_')

                        exon_start = int(exon_start) - int(entry[2])
                        exon_stop = int(exon_stop) - int(entry[2])

                        feature = SeqFeature(FeatureLocation(exon_start, exon_stop), strand=+1)
                        gds_features.add_feature(feature, name="Exon", label=False, color="grey", label_size=22)

                gdd.draw(format='circular', pagesize=(600, 600), circle_core=0.6, track_size=0.3, tracklines=0, x=0.00,
                         y=0.00, start=0, end=circrna_length-1)

                gdd.write(self.output_dir + "/" + circular_rna_id_isoform + ".svg", "svg")
