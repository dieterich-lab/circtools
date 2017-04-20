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

import logging
import time
import os
import sys
import re
import multiprocessing
import functools

import pybedtools
import circ_module.circ_template


class EnrichmentModule(circ_module.circ_template.CircTemplate):
    def __init__(self, argparse_arguments, program_name, version):

        # get the user supplied options
        self.cli_params = argparse_arguments
        self.program_name = program_name
        self.version = version

    def run_module(self):

        # set time format
        time_format = time.strftime("%Y_%m_%d__%H_%M")
        # time_format = ""

        # set up the multiprocessing pool for multi-threading
        mp_pool = multiprocessing.Pool(processes=self.cli_params.num_processes)

        # let's first check if the output directory exists
        if not (os.access(self.cli_params.output_directory, os.W_OK)):
            print("Output directory %s not writable." % self.cli_params.output_directory)
            # exit with -1 error if we can't use it
            exit(-1)

        # let's first check if the temporary directory exists
        if not (os.access(self.cli_params.tmp_directory, os.W_OK)):
            print("Temporary directory %s not writable." % self.cli_params.tmp_directory)
            # exit with -1 error if we can't use it
            exit(-1)

        # set temporary directory for pybedtools
        pybedtools.set_tempdir(self.cli_params.tmp_directory)

        # starting up logging system
        logging.basicConfig(filename=os.path.join(self.cli_params.output_directory,
                                                  self.program_name + "__" + time_format + ".log"),
                            filemode="w",
                            level=logging.DEBUG,
                            format="%(asctime)s %(message)s"
                            )

        logging.info("%s %s started" % (self.program_name, self.version))
        logging.info("%s command line: %s" % (self.program_name, " ".join(sys.argv)))

        # check if input files are in place, we don't wan't to fail later
        self.check_input_files([self.cli_params.bed_input,
                                self.cli_params.circ_rna_input,
                                self.cli_params.annotation,
                                self.cli_params.genome_file])

        # read in CLIP peaks
        supplied_bed = self.read_bed_file(self.cli_params.bed_input)

        # read in annotation
        annotation_bed = self.read_annotation_file(self.cli_params.annotation)

        gene_annotation_file = self.cli_params.output_directory + '/' + \
                               os.path.basename(self.cli_params.annotation) + \
                               '_genes.bed'

        annotation_bed.saveas(gene_annotation_file)

        # read in circular RNA predictions from DCC
        circ_rna_bed = self.read_circ_rna_file(self.cli_params.circ_rna_input, annotation_bed)

        # do circle saves
        circle_annotation_file = self.cli_params.output_directory + '/' + os.path.basename(
            self.cli_params.circ_rna_input) + '_circles.bed'

        circ_rna_bed.saveas(circle_annotation_file)

        # create list of shuffled peaks
        shuffled_peaks_linear = mp_pool.map(functools.partial(self.shuffle_peaks_through_genome,
                                                              bed_file=supplied_bed,
                                                              annotation=gene_annotation_file,
                                                              genome_file=self.cli_params.genome_file),
                                            range(self.cli_params.num_iterations))

        shuffled_peaks_linear.append(supplied_bed)

        # create list of shuffled peaks
        shuffled_peaks_circular = mp_pool.map(functools.partial(self.shuffle_peaks_through_genome,
                                                                bed_file=supplied_bed,
                                                                annotation=circle_annotation_file,
                                                                genome_file=self.cli_params.genome_file),
                                              range(self.cli_params.num_iterations))

        shuffled_peaks_circular.append(supplied_bed)

        # do the intersections
        results = mp_pool.map(functools.partial(self.random_sample_step,
                                                circ_rna_bed=circ_rna_bed,
                                                annotation_bed=annotation_bed,
                                                shuffled_peaks_linear=shuffled_peaks_linear,
                                                shuffled_peaks_circular=shuffled_peaks_circular
                                                ), range(self.cli_params.num_iterations + 1))

        result_table = self.generate_count_table(results)

        result_file = self.cli_params.output_directory + "/output_" + time_format + ".csv"

        with open(result_file, "w") as text_file:
            print(result_table, end="", file=text_file)

        # ------------------------------------- Function definitions start here ---------------------------------------

    def read_circ_rna_file(self, circ_rna_input, annotation_bed):
        """Reads a CircCoordinates file produced by DCC
        Will halt the program if file not accessible
        Returns a BedTool object
        """
        self.log_entry("Parsing circular RNA input file...")
        try:
            file_handle = open(circ_rna_input)
        except PermissionError:
            message = ("Input file " + str(circ_rna_input) + " cannot be read, exiting.")
            logging.info(message)
            sys.exit(message)
        else:
            with file_handle:
                line_iterator = iter(file_handle)
                # skip first line with the header
                # we assume it's there (DCC default)
                next(line_iterator)
                bed_content = ""
                bed_entries = 0
                bed_peak_sizes = 0
                for line in line_iterator:
                    columns = line.split('\t')

                    # extract chromosome, start, stop, gene name, and strand
                    entry = [self.strip_chr_name(columns[0]), columns[1], columns[2], columns[3], "0", columns[5]]
                    # print(test)

                    # concatenate lines to one string
                    bed_content += '\t'.join(entry) + "\n"
                    # bed_content += str(test)

                    # sys.stdout.write("Processing circRNA # %s \r" % bed_entries)
                    # sys.stdout.flush()

                    bed_entries += 1
                    bed_peak_sizes += (int(columns[2]) - int(columns[1]))

            # create a "virtual" BED file
            virtual_bed_file = pybedtools.BedTool(bed_content, from_string=True)
            print(virtual_bed_file)

            test = annotation_bed.intersect(virtual_bed_file, s=True, wb=True)

        self.log_entry("Done parsing circular RNA input file:")
        self.log_entry("=> %s circular RNAs, %s nt average (theoretical unspliced) length" %
                       (bed_entries, round(bed_peak_sizes / bed_entries)))

        return test

    def read_bed_file(self, bed_input):
        """Reads a BED file
        Will halt the program if file not accessible
        Returns a BedTool object
        """
        self.log_entry("Parsing BED input file...")
        try:
            file_handle = open(bed_input)
        except PermissionError:
            message = ("Input file " + str(bed_input) + " cannot be read, exiting.")
            logging.info(message)
            sys.exit(message)
        else:
            with file_handle:
                line_iterator = iter(file_handle)
                bed_content = ""
                bed_entries = 0
                bed_peak_sizes = 0
                for line in line_iterator:
                    columns = line.split('\t')
                    # extract chr, start, stop, score(0), name, strand
                    entry = [self.strip_chr_name(columns[0]), columns[1], columns[2], columns[3], columns[4],
                             columns[5]]
                    # concatenate lines to one string
                    bed_content += '\t'.join(entry) + "\n"
                    bed_entries += 1
                    bed_peak_sizes += (int(columns[2]) - int(columns[1]))

            # create a "virtual" BED file
            virtual_bed_file = pybedtools.BedTool(bed_content, from_string=True)

        self.log_entry("Done parsing BED input file:")
        self.log_entry("=> %s peaks, %s nt average width" %
                       (bed_entries, round(bed_peak_sizes / bed_entries)))

        return virtual_bed_file

    def extract_gene_name_from_gtf(self, annotation_string):
        """Gets the 8th column of a GTF file
        Returns the gene name
        """
        # create a group match in the "..." of gene_name
        pattern = re.compile(".*gene_name\s\"([a-zA-Z0-9_\-./]+)\"")
        match = pattern.match(annotation_string)
        if match:
            return match.group(1)
        else:
            self.log_entry("Problem parsing gene name; offending entry is: %s" % annotation_string)
            return "GENE_NAME_PARSE_ERROR"

    @staticmethod
    def strip_chr_name(chromosome_string):
        """Removes "chr" in an non case-sensitive manner
        Returns only the chromosome number
        """

        match = re.compile(re.escape('chr'), re.IGNORECASE)
        chromosome_name = match.sub('', chromosome_string)

        return chromosome_name

    def read_annotation_file(self, annotation_file):
        """Reads a GTF file
        Will halt the program if file not accessible
        Returns a BedTool object only containing gene sections
        """
        self.log_entry("Parsing annotation...")

        try:
            file_handle = open(annotation_file)
        except PermissionError:
            message = ("Input file " + str(annotation_file) + " cannot be read, exiting.")
            logging.info(message)
            sys.exit(message)
        else:
            with file_handle:
                line_iterator = iter(file_handle)
                bed_content = ""
                gene_entry = 1
                for line in line_iterator:
                    # we skip any comment lines
                    if line.startswith("#"):
                        continue

                    # split up the annotation line
                    columns = line.split('\t')

                    # we only want the coordinates of the gene entries
                    if not (columns[2] == "exon"):
                        continue

                    # columns 8 is the long annotation string from GTF
                    gene_name = self.extract_gene_name_from_gtf(columns[8])

                    # extract chromosome, start, stop, score(0), name and strand
                    # we hopefully have a gene name now and use this one for the entry
                    entry = [self.strip_chr_name(columns[0]), columns[3], columns[4], gene_name, str(0), columns[6], ]

                    # concatenate lines to one string
                    bed_content += '\t'.join(entry) + "\n"

                    sys.stdout.write("Processing exon # %s \r" % gene_entry)
                    sys.stdout.flush()

                    gene_entry += 1

                # "escape the \r from counting output"
                sys.stdout.write("\n")

                # count will be increased one more time even if done - so we subtract 1
                self.log_entry("Processed %s genes" % (gene_entry - 1))

            # create a "virtual" BED file
            virtual_bed_file = pybedtools.BedTool(bed_content, from_string=True)

        self.log_entry("Done parsing annotation")

        return virtual_bed_file

    def shuffle_peaks_through_genome(self, iteration, bed_file, annotation, genome_file):
        """Gets a (virtual) BED files and shuffle its contents throughout the supplied genome
        Will only use supplied annotation for features (in our case only transcript regions)
        """

        self.log_entry("Starting shuffling thread %d" % iteration)
        shuffled_bed = bed_file.shuffle(g=genome_file, chrom=True, incl=annotation)
        self.log_entry("Finished shuffling thread %d" % iteration)

        return shuffled_bed

    @staticmethod
    def do_intersection(query_bed, base_bed):
        """Gets two bed files (supplied peaks and circle coordinates) and does an intersection
        """

        # we employ the c=true parameter to directly get the counts as part of the results
        intersect_return = base_bed.intersect(query_bed, c=True)

        return intersect_return

    @staticmethod
    def process_intersection(intersection_output, tag):
        """Gets two bed files (supplied peaks and circle coordinates) and does an intersection
        """
        # initialise empty dict
        count_table = {}

        feature_iterator = iter(intersection_output)
        for bed_feature in feature_iterator:
            key = tag + "\t" + bed_feature.name +\
                  "\t" + bed_feature.chrom + "_" +\
                  str(bed_feature.start) + "_" +\
                  str(bed_feature.stop) +\
                  bed_feature.strand
            count_table[key] = bed_feature[6]  # [6] == number of "hits"

        return count_table

    @staticmethod
    def generate_count_table(count_table):
        """Gets two bed files (supplied peaks and circle coordinates) and does an intersection
        """

        all_keys = []

        for iteration in range(len(count_table)):
            all_keys += list(count_table[iteration][0])
            all_keys += list(count_table[iteration][1])

        all_keys_unique = set(all_keys)

        result_table = ""

        for key in all_keys_unique:
            result_table += key + "\t"

            for i in range(len(count_table)):

                if "lin" in key:
                    if key in count_table[i][0]:
                        result_table += count_table[i][0][key] + "\t"
                    else:
                        result_table += "0\t"

                if "circ" in key:
                    if key in count_table[i][1]:
                        result_table += count_table[i][1][key] + "\t"
                    else:
                        result_table += "0\t"

            result_table += "\n"

        return result_table

    def random_sample_step(self, iteration, circ_rna_bed, annotation_bed, shuffled_peaks_linear,
                           shuffled_peaks_circular):
        """Logs to log file and prints on screen
        """
        self.log_entry("Starting intersection thread %d" % iteration)

        # get circular and linear intersect
        circular_intersect = self.do_intersection(shuffled_peaks_circular[iteration], circ_rna_bed)
        linear_intersect = self.do_intersection(shuffled_peaks_linear[iteration], annotation_bed)

        # process results of the intersects
        linear_count_table = self.process_intersection(linear_intersect, "lin")
        circular_count_table = self.process_intersection(circular_intersect, "circ")

        self.log_entry("Finished intersection thread %d" % iteration)

        return linear_count_table, circular_count_table

    def module_name(self):
        """"Return a string representing the name of the module."""
        return self.program_name
