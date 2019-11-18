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

import functools
import logging
import multiprocessing
import os
import re
import sys
import time

import circ_module.circ_template
import pybedtools

FILE_TYPE_GTF = ".gtf"
FILE_TYPE_BED = ".bed"


class EnrichmentModule(circ_module.circ_template.CircTemplate):
    def __init__(self, argparse_arguments, program_name, version):

        # get the user supplied options
        self.cli_params = argparse_arguments
        self.program_name = program_name
        self.version = version
        self.observed_counts = []
        self.results = []
        self.tmp_dict = {}
        self.phase_storage = {}
        self.virtual_inclusion_file_path = "all"
        self.virtual_inclusion_object = None
        self.circRNA_buddies = {}
        self.whitelist_fg = ""
        self.whitelist_bg = ""

    def run_module(self):

        # set time format
        time_format = time.strftime("%Y_%m_%d__%H_%M")

        # let's first check if the temporary directory exists
        if not os.path.exists(self.cli_params.tmp_directory):
            os.makedirs(self.cli_params.tmp_directory)

        # let's first check if the temporary directory exists
        if not os.path.exists(self.cli_params.output_directory):
            os.makedirs(self.cli_params.output_directory)

        # let's first check if the output directory exists
        if not (os.access(self.cli_params.output_directory, os.W_OK)):
            self.log_entry("Output directory %s not writable." % self.cli_params.output_directory)
            # exit with -1 error if we can't use it
            exit(-1)

        # let's first check if the temporary directory exists
        if not (os.access(self.cli_params.tmp_directory, os.W_OK)):
            self.log_entry("Temporary directory %s not writable." % self.cli_params.tmp_directory)
            # exit with -1 error if we can't use it
            exit(-1)

        # holds the temporary bed file content
        temp_bed = ""

        # starting up logging system
        print(os.path.join(self.cli_params.output_directory,
                                                  self.cli_params.output_filename + "_" + time_format + ".log"))

        logging.basicConfig(filename=os.path.join(self.cli_params.output_directory,
                                                  self.cli_params.output_filename + "_" + time_format + ".log"),
                            filemode="w",
                            level=logging.DEBUG,
                            format="%(asctime)s %(message)s"
                            )

        logging.info("%s %s started" % (self.program_name, self.version))
        logging.info("%s command line: %s" % (self.program_name, " ".join(sys.argv)))

        import subprocess

        try:
            bedtools_version = subprocess.check_output("bedtools --version",
                                                       shell=True,
                                                       stderr=subprocess.DEVNULL).decode("utf-8")
            bedtools_version = bedtools_version.rstrip()

            self.log_entry("%s detected" % bedtools_version)

            bedtools_version = int(bedtools_version.split(".")[1])

        except subprocess.CalledProcessError:
            self.log_entry("Bedtools binary not found in your path, exiting")
            exit(-1)

        # (default is ["all"])
        if self.cli_params.include_features:

            # for each of the user supplied include features
            for feature_type in self.cli_params.include_features:
                temp_bed += self.read_annotation_file(self.cli_params.annotation, entity=feature_type, string=True)

            # we create a bed file on disk for all features
            self.virtual_inclusion_object = pybedtools.BedTool(temp_bed, from_string=True)

            if bedtools_version > 26:
                self.virtual_inclusion_object = self.virtual_inclusion_object.sort().merge(s=True,  # strand specific
                                                                                       c="4,5,6",  # copy columns 5 & 6
                                                                                       o="distinct,distinct,distinct")  # group
            else:
                self.virtual_inclusion_object = self.virtual_inclusion_object.sort().merge(s=True,  # strand specific
                                                                                       c="5,6",  # copy columns 5 & 6
                                                                                       o="distinct,distinct")  # group
            # set the path where to store it
            self.virtual_inclusion_file_path = self.cli_params.output_directory + \
                                               '/' + self.cli_params.output_filename + \
                                               "_" + \
                                               os.path.basename(self.cli_params.annotation) + \
                                               '_features.bed'
            # save temporary bedtools file as "real" file
            self.virtual_inclusion_object.saveas(self.virtual_inclusion_file_path)
        else:
            # fixed keyword to disable inclusion of features
            self.virtual_inclusion_file_path = "all"

        # set temporary directory for pybedtools
        pybedtools.set_tempdir(self.cli_params.tmp_directory)

        # check if input files are in place, we don't wan't to fail later
        self.check_input_files([self.cli_params.bed_input,
                                self.cli_params.circ_rna_input,
                                self.cli_params.annotation,
                                self.cli_params.genome_file])

        # read in CLIP peaks
        supplied_bed = self.read_bed_file(self.cli_params.bed_input)

        # read in annotation
        annotation_bed = self.read_annotation_file(self.cli_params.annotation)

        gene_annotation_file = self.cli_params.output_directory +\
                               '/' + self.cli_params.output_filename +\
                               "_" + \
                               os.path.basename(self.cli_params.annotation) +\
                               '_genes.bed'

        # if we are in "feature" mode we have to get all exons for the linear genes
        if self.virtual_inclusion_file_path != "all":
            annotation_bed = self.virtual_inclusion_object.intersect(annotation_bed, s=True, wb=True)

        annotation_bed.saveas(gene_annotation_file)
        # read in circular RNA predictions from DCC
        circ_rna_bed = self.read_circ_rna_file(self.cli_params.circ_rna_input,
                                               annotation_bed)

        # do circle saves
        circle_annotation_file = self.cli_params.output_directory +\
                                 '/' + self.cli_params.output_filename +\
                                 "_" +\
                                 os.path.basename(self.cli_params.circ_rna_input) +\
                                 '_circles.bed'

        # if we are in "feature" mode we have to get all exons for the circRNAs
        if self.virtual_inclusion_file_path != "all":
            circ_rna_bed = self.virtual_inclusion_object.intersect(circ_rna_bed, s=True, wb=True)

        circ_rna_bed.saveas(circle_annotation_file)

        if self.cli_params.whitelist:
            # self.whitelist = self.generate_location_hash(self.cli_params.whitelist)
            self.whitelist_fg = pybedtools.BedTool(self.cli_params.whitelist)
            self.whitelist_fg = self.virtual_inclusion_object.intersect(self.whitelist_fg, f=1.0, wb=True)

            tmp = ""
            for line in str(self.whitelist_fg).splitlines():
                bed_feature = line.split('\t')
                bed_feature[9] = bed_feature[3]
                tmp += "\t".join(bed_feature)+"\n"

            self.whitelist_fg = pybedtools.BedTool(tmp, from_string=True)

            self.whitelist_fg.saveas("test2.bed")
            # exit(0)
            # generate inverse list of exons: those that are not enriched: i.e. background exons
            # self.whitelist_bg = circ_rna_bed #.intersect(self.whitelist_fg, v=True)


        # set up the multiprocessing pool for multi-threading
        mp_pool = multiprocessing.Pool(processes=self.cli_params.num_processes)

        # create list of shuffled peaks
        self.log_entry("Starting random shuffling of input peaks")
        shuffled_peaks = (mp_pool.map(functools.partial(self.shuffle_peaks_through_genome,
                                                        bed_file=supplied_bed,
                                                        genome_file=self.cli_params.genome_file,
                                                        include=self.virtual_inclusion_file_path),
                                      range(self.cli_params.num_iterations)))

        shuffled_peaks.insert(0, supplied_bed)

        # shuffled_peaks.append(supplied_bed)
        # do the intersections
        self.log_entry("Starting data acquisition from samplings")
        self.results = mp_pool.map(functools.partial(self.random_sample_step,
                                                     circ_rna_bed=circ_rna_bed,
                                                     annotation_bed=annotation_bed,
                                                     shuffled_peaks=shuffled_peaks,
                                                     ), range(self.cli_params.num_iterations + 1))

        self.tmp_dict = self.process_intersection(self.results[0][0])

        self.observed_counts = (
            self.tmp_dict,
            self.process_intersection(self.results[0][1], linear_start=True)
        )
        import pprint
        pp = pprint.PrettyPrinter(indent=4)
        pp.pprint(self.observed_counts[0])
        print("------------")
        pp.pprint(self.observed_counts[1])

        # how many iterations do we want to do before cleaning up?
        iterations_per_phase = int(self.cli_params.num_iterations / self.cli_params.num_processes)

        # we probably have some leftover iterations we have to take care of
        leftover_iterations = (self.cli_params.num_iterations % self.cli_params.num_processes)

        # how many phases do we have?
        num_phases = self.cli_params.num_processes

        # add one to compensate for remaining permutations
        if leftover_iterations > 0:
            num_phases += 1

        # we'll store the final results here
        self.phase_storage = {}

        # we need some manual garbage collection to keep memory profile lower
        # import gc

        self.log_entry("Cleaning up... just a second")

        # before starting permutation test, clean up
        # gc.collect()

        # for each phase we perform iterations_per_phase tests
        # the remainder is done in an extra run
        for phase in range(0, num_phases):

            self.log_entry("Starting permutation test phase %d" % (phase+1))

            start_id = (phase*iterations_per_phase)+1
            stop_id = ((phase+1)*iterations_per_phase)+1

            # compensate for the first round so that we start correctly at ID 0
            if phase == 0:
                start_id -= 1

            # results of one phase of the computation
            intermediate_result = mp_pool.map(functools.partial(self.run_permutation_test), range(start_id, stop_id))

            # how many results did we get back?
            num_results = len(intermediate_result)

            # we now have to convert the temporary dicts into out main dict to save memory
            for current_num in range(0, num_results):

                # for each gene per permutation
                for gene in intermediate_result[current_num]:

                    # initialize with 2 empty dicts for linear and circular RNA
                    if gene not in self.phase_storage:
                        self.phase_storage[gene] = {0: {}, 1: {}}

                    # for both RNA types
                    for rna in range(0, 2):
                        if rna in intermediate_result[current_num][gene]:
                            # get the location
                            for location_key in self.observed_counts[rna][gene]:

                                if location_key not in self.phase_storage[gene][rna]:
                                    self.phase_storage[gene][rna][location_key] = 0

                                # only if the temporary data says true (= we saw more peaks than observed)
                                if location_key in intermediate_result[current_num][gene][rna]:
                                    self.phase_storage[gene][rna][location_key] += 1

            self.log_entry("Cleaning up... just a second")

        # generate the result table
        result_table = self.print_results()

        # and print it to a file
        result_file = self.cli_params.output_directory + "/" +\
                      self.cli_params.output_filename + "_" +\
                      str(self.cli_params.num_iterations) + "_" +\
                      time_format + ".csv"

        with open(result_file, 'w') as text_file:
            text_file.write(result_table)

        if not self.cli_params.keep_temp:
            self.clean_up_temp_files()

        # ------------------------------------- Function definitions start here ---------------------------------------

    def read_circ_rna_file(self, circ_rna_input, annotation_bed):
        """Reads a CircCoordinates file produced by DCC
        Will halt the program if file not accessible
        Returns a BedTool object
        """
        self.log_entry("Parsing circular RNA input file...")

        # set temporary directory for pybedtools
        pybedtools.set_tempdir(self.cli_params.tmp_directory)

        try:
            file_handle = open(circ_rna_input)
        except PermissionError:
            message = ("Input file " + str(circ_rna_input) + " cannot be read, exiting.")
            logging.info(message)
            sys.exit(message)
        else:
            with file_handle:
                line_iterator = iter(file_handle)

                bed_content = ""
                bed_entries = 0
                bed_peak_sizes = 0
                last_buddy = "NULL"
                last_key = "NULL"

                for line in line_iterator:
                    columns = line.split('\t')

                    # skip first line with the header
                    # we assume it's there (DCC default)
                    if columns[0] == "Chr" and columns[1] == "Start":
                        continue

                    # extract chromosome, start, stop, gene name, and strand
                    entry = [self.strip_chr_name(columns[0]), columns[1], columns[2], columns[3], "0", columns[5]]

                    # concatenate lines to one string
                    bed_content += '\t'.join(entry) + "\n"

                    bed_entries += 1
                    bed_peak_sizes += (int(columns[2]) - int(columns[1]))

                    # column7 -> reserved for flanking intron detection
                    # update: make sure column7 has the format we expect

                    if len(columns) > 7:


                        components = columns[7].split("_")
                        # we have a normal key without feature information
                        if len(components) == 4 and len(columns) <= 7:

                            # remove \n
                            columns[7] = columns[7].rstrip()

                            # generate key for this buddy
                            this_key = self.strip_chr_name(columns[0]) + "_" + columns[1] + "_" + columns[2] + "_" + columns[5] + "_" + str(int(columns[2])-int(columns[1])) + "_1"

                            if last_buddy != "NULL" and last_buddy == columns[7]:

                                # save the buddy key pair
                                self.circRNA_buddies[this_key] = last_key
                                self.circRNA_buddies[last_key] = this_key

                            # reset buddy
                            last_buddy = columns[7]

                            # generate key for this buddy
                            last_key = self.strip_chr_name(columns[0]) + "_" + columns[1] + "_" + columns[2] + "_" + columns[5] + "_" + str(int(columns[2])-int(columns[1])) + "_1"

            self.log_entry("Done parsing circular RNA input file:")
            self.log_entry("=> %s circular RNAs, %s nt average (theoretical unspliced) length" %
                           (bed_entries, round(bed_peak_sizes / bed_entries)))

            # create a "virtual" BED file
            virtual_bed_file = pybedtools.BedTool(bed_content, from_string=True)

            # if user supplied a feature we are in "feature mode"
            # this means some different treatment downstream

            # But for now "all" means we use gene as feature -> classic mode
            if self.virtual_inclusion_file_path == "all":
                return annotation_bed.intersect(virtual_bed_file, s=True)

            # This is the feature mode
            # we don't intersect with the annotation since we do this in a different step
            else:
                return virtual_bed_file

    def read_bed_file(self, bed_input):
        """Reads a BED file
        Will halt the program if file not accessible
        Returns a BedTool object
        """

        # set temporary directory for pybedtools
        pybedtools.set_tempdir(self.cli_params.tmp_directory)

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

    def read_annotation_file(self, annotation_file, entity="gene", string=False):
        """Reads a GTF file
        Will halt the program if file not accessible
        Returns a BedTool object only containing gene sections
        """

        # set temporary directory for pybedtools
        pybedtools.set_tempdir(self.cli_params.tmp_directory)

        self.log_entry("Parsing annotation...")

        try:
            file_handle = open(annotation_file)
        except PermissionError:
            message = ("Input file " + str(annotation_file) + " cannot be read, exiting.")
            logging.info(message)
            sys.exit(message)
        else:
            # check file type: BED or GTF allowed
            file_type = os.path.splitext(annotation_file)[1].lower()
            if file_type != FILE_TYPE_GTF and file_type != FILE_TYPE_BED:
                message = ("Input file " + str(annotation_file) + " is neither BED nor GTF format!")
                logging.info(message)
                sys.exit(message)
            with file_handle:
                entity_list = []
                line_iterator = iter(file_handle)
                bed_content = ""
                entity_count = 1
                for line in line_iterator:
                    # we skip any comment lines
                    if line.startswith("#"):
                        continue

                    # split up the annotation line
                    columns = line.split('\t')

                    # we have no "type" in BED
                    if file_type == FILE_TYPE_GTF:
                        if columns[2] not in entity_list:
                            entity_list.append(columns[2])

                    # we have no "type" in BED
                    if file_type == FILE_TYPE_GTF:
                        # we only want the coordinates of the gene entries
                        if not (columns[2] == entity):
                            continue

                    # we do not want any 0-length intervals -> bedtools segfault
                    if file_type == FILE_TYPE_GTF and int(columns[4]) - int(columns[3]) == 0:
                            continue
                    if file_type == FILE_TYPE_BED and int(columns[1]) - int(columns[2]) == 0:
                            continue

                    # columns 8 is the long annotation string from GTF
                    if file_type == FILE_TYPE_GTF:
                        gene_name = self.extract_gene_name_from_gtf(columns[8])
                    else:
                        gene_name = columns[3]  # somewhat easier in BED...

                    # extract chromosome, start, stop, score(0), name and strand
                    # we hopefully have a gene name now and use this one for the entry

                    if file_type == FILE_TYPE_GTF:
                        entry = [
                            self.strip_chr_name(columns[0]),
                            columns[3],
                            columns[4],
                            gene_name,
                            str(0),
                            columns[6],
                        ]
                    else:
                        entry = [
                            self.strip_chr_name(columns[0]),
                            columns[1],
                            columns[2],
                            gene_name,
                            str(0),
                            columns[5]
                        ]

                    # concatenate lines to one string
                    bed_content += '\t'.join(entry) + "\n"

                    entity_count += 1

                # count will be increased one more time even if done - so we subtract 1
                if file_type is FILE_TYPE_GTF:
                    self.log_entry("Found %s entries of type %s" % (entity_count - 1, entity))
                else:
                    self.log_entry("Found %s entries " % (entity_count -1))

            self.log_entry("Done parsing annotation")

            if not bed_content:
                self.log_entry("No features of type \"%s\" found in annotation file." % entity)
                self.log_entry("Available entities found in the annotation file:\n%s" % entity_list)
                self.log_entry("Exiting now.")
                exit(-1)

            if string:
                return bed_content
            else:
                # create a "virtual" BED file
                virtual_bed_file = pybedtools.BedTool(bed_content, from_string=True)
                return virtual_bed_file

    def shuffle_peaks_through_genome(self, iteration, bed_file, genome_file, include):
        """Gets a (virtual) BED files and shuffle its contents throughout the supplied genome
        Will only use supplied annotation for features (in our case only transcript regions)
        """
        # set temporary directory for pybedtools
        pybedtools.set_tempdir(self.cli_params.tmp_directory)

        self.log_entry("Processing shuffling thread %d" % (iteration+1))

        if include == "all":
            shuffled_bed = bed_file.shuffle(g=genome_file)
        else:
            shuffled_bed = bed_file.shuffle(g=genome_file, incl=include)

        return shuffled_bed

    def do_intersection(self, query_bed, base_bed):
        """Gets two bed files (supplied peaks and circle coordinates) and does an intersection
        """
        # set temporary directory for pybedtools
        pybedtools.set_tempdir(self.cli_params.tmp_directory)

        # we employ the c=true parameter to directly get the counts as part of the results

        if self.virtual_inclusion_file_path != "all":
            # we need to do a second intersect to get the amount of peaks that are located within circular RNAs _AND_
            # also part of a features: e.g. an exon of a circRNA

            intersect_return = base_bed.intersect(query_bed, c=True, s=True)
            intersect_return = self.pre_process_intersection(intersect_return)
        else:
            intersect_return = base_bed.intersect(query_bed, c=True, s=True)

        return intersect_return

    def pre_process_intersection(self, intersection_input):
        """Processes intersection output generated by do_intersection() in case of feature-based shuffling
        Returns a count table for the given intersection
        """
        # initialise empty dict
        return_string = ""

        # holds the different circRNA isoforms for an annotated host gene
        # once we see another host gene annotation, all isoforms are written in the bed string
        isoform_net_length = {}
        isoform_num_features = {}
        isoform_count = {}
        isoform_name = {}

        for line in str(intersection_input).splitlines():
            bed_feature = line.split('\t')

            # we add up the length of each feature that is part of the "uber" feature, e.g. sum up exon length
            # as sub features of a gene

            key = bed_feature[6] + "_" + bed_feature[7] + "_" + bed_feature[8] + "_" + bed_feature[11]

            # we also count how many features are part of this entity
            if key not in isoform_net_length:
                isoform_net_length[key] = (int(bed_feature[2]) - int(bed_feature[1]))
                isoform_num_features[key] = 1
                isoform_count[key] = int(bed_feature[12])
                isoform_name[key] = bed_feature[9]

            elif bed_feature[1] != bed_feature[7] or bed_feature[2] != bed_feature[8]:
                isoform_net_length[key] += (int(bed_feature[2]) - int(bed_feature[1]))
                isoform_num_features[key] += 1
                isoform_count[key] += int(bed_feature[12])

        for feature_key in isoform_net_length:
            data = self.decode_location_key(feature_key)

            new_line = data["chr"] + "\t" \
                       + str(data["start"]) + "\t" \
                       + str(data["stop"]) + "\t" \
                       + str(isoform_name[feature_key]) + "\t" \
                       + str(isoform_net_length[feature_key]) + "_" + \
                       str(isoform_num_features[feature_key]) + "\t" + \
                       data["strand"] + "\t" \
                       + str(isoform_count[feature_key])

            return_string += new_line + "\n"

        # reset stuff
        isoform_net_length.clear()
        isoform_num_features.clear()
        isoform_count.clear()
        isoform_name.clear()

        return pybedtools.BedTool(return_string, from_string=True)

    def process_intersection(self, intersection_input, normalize=False, linear_start=False):
        """Processes intersection output generated by do_intersection()
        Returns a count table for the given intersection
        """

        # initialise empty dict
        count_table = {}

        # wrapper function to transparently allow for normalization
        def normalize_count(start, stop, count):
            if normalize and count > 0:
                return count / (stop - start)
            elif normalize and count == 0 * 1000:
                return 0
            else:
                return count

        for line in str(intersection_input).splitlines():

            bed_feature = line.split('\t')
            # key has the form: chromosome_start_stop[strand]
            key = bed_feature[0] + "_" + \
                str(bed_feature[1]) + "_" + \
                str(bed_feature[2]) + "_" + \
                bed_feature[5]

            gene_name = bed_feature[3]

            # in feature mode, we extend the key by count and feature length
            if self.virtual_inclusion_file_path != "all":
                key += "_" + str(bed_feature[4])

            # print(bed_feature)

            if linear_start and gene_name in self.tmp_dict:
                # need to get list of circRNAs linked to this host gene

                for location_key_circular in self.tmp_dict[gene_name]:
                    tmp_data = self.decode_location_key(location_key_circular)

                    circ_count = self.tmp_dict[gene_name][location_key_circular]

                    key = (key + "_" +
                           tmp_data["chr"] + "_" +
                           str(tmp_data["start"]) + "_" +
                           str(tmp_data["stop"]) + "_" +
                           tmp_data["strand"] + "_" +
                           str(tmp_data["feature_length"]) + "_" +
                           str(tmp_data["feature_count"])
                           )

                    # we have to create the nested dictionaries if not already existing
                    if gene_name not in count_table:
                        count_table[gene_name] = {}
                        # sometime gene names are not unique
                        # we simply check if the last key is at least of the same chromosome
                    else:
                        first_chr = self.decode_location_key(next(iter(count_table[gene_name].keys())))["chr"]

                        # check if chromosomes match
                        if first_chr != bed_feature[0]:
                            continue
                        #     gene_name = gene_name + "_" + bed_feature[0] + "_" + bed_feature[1]
                        #     count_table[gene_name] = {}

                    if key not in count_table[gene_name]:
                        count_table[gene_name][key] = {}

                    # set the appropriate dict entry
                    count_table[gene_name][key] = normalize_count(bed_feature[1],
                                                                       bed_feature[2],
                                                                       int(bed_feature[6])
                                                                       ) - circ_count
                    # key has the form: chromosome_start_stop[strand]
                    key = bed_feature[0] + "_" + \
                          str(bed_feature[1]) + "_" + \
                          str(bed_feature[2]) + "_" + \
                          bed_feature[5]

                    # in feature mode, we extend the key by count and feature length
                    if self.virtual_inclusion_file_path != "all":
                        key += "_" + str(bed_feature[4])

            else:

                # we have to create the nested dictionaries if not already existing
                if gene_name not in count_table:
                    count_table[gene_name] = {}
                    # sometime gene names are not unique
                    # we simply check if the last key is at least of the same chromosome
                else:
                    first_chr = self.decode_location_key(next(iter(count_table[gene_name].keys())))["chr"]

                    # check if chromosomes match
                    if first_chr != bed_feature[0]:
                        continue
                    # gene_name = gene_name + "_" + bed_feature[0] + "_" + bed_feature[1]
                    # count_table[gene_name] = {}

                if key not in count_table[gene_name]:
                    count_table[gene_name][key] = {}

                # set the appropriate dict entry
                count_table[gene_name][key] = normalize_count(bed_feature[1],
                                                                   bed_feature[2],
                                                                   int(bed_feature[6])
                                                                   )
        # return one unified count table
        return count_table

    @staticmethod
    def decode_location_key(key):
        """Decodes a given location key
        Returns a dictionary of decoded values
        """
        key_components = key.split("_")
        # we have a normal key without feature information
        if len(key_components) == 4:

            data = {
                    "chr": key_components[0],
                    "start": int(key_components[1]),
                    "stop": int(key_components[2]),
                    "strand": key_components[3],
                    "length": (int(key_components[2])-int(key_components[1])),
                    }
        # we have an extended key with feature information and will return also those
        elif len(key_components) == 6:
            data = {
                    "chr": key_components[0],
                    "start": int(key_components[1]),
                    "stop": int(key_components[2]),
                    "strand": key_components[3],
                    "length": (int(key_components[2])-int(key_components[1])),
                    "feature_length": int(key_components[4]),
                    "feature_count": int(key_components[5])
                    }

        else:
            data = {
                    "chr": key_components[0],
                    "start": int(key_components[1]),
                    "stop": int(key_components[2]),
                    "strand": key_components[3],
                    "length": (int(key_components[2])-int(key_components[1])),
                    "feature_length": int(key_components[4]),
                    "feature_count": int(key_components[5]),
                    "chr_circ": key_components[6],
                    "start_circ": int(key_components[7]),
                    "stop_circ": int(key_components[8]),
                    "strand_circ": key_components[9],
                    "feature_length_circ": int(key_components[10]),
                    "feature_count_circ": int(key_components[11]),
                    "circ_data": key_components[6] + "_" +
                                 key_components[7] + "_" +
                                 key_components[8] + "_" +
                                 key_components[9] + "_" +
                                 key_components[10] + "_" +
                                 key_components[11]
                    }
        return data

    def linear_length_wo_circ(self, key_circ, key_linear):
        """Computes the remaining linear RNA length given the circular RNA location key
        Returns a list of circular RNA length [0] and linear RNA length [1]
        Only used in non-feature mode
        """
        circ = self.decode_location_key(key_circ)
        linear = self.decode_location_key(key_linear)
        length_circ = (circ["stop"] - circ["start"])
        length_lin = (linear["stop"] - linear["start"]) - length_circ
        return length_circ, length_lin

    def get_extended_key_data(self, key):
        """Computes the remaining linear RNA length given the circular RNA location key
        Returns a list of circular RNA length [0] and linear RNA length [1]
        """
        circ = self.decode_location_key(key)

        return {"feature_length": int(circ["feature_length"]), "feature_count": int(circ["feature_count"])}

    @staticmethod
    def normalize_count(length, count):
        """Normalizes a given absolute count of peaks to a given length
        Returns the length-normalized count
        """
        if count > 0 and length > 0:
            return (count / length) * 1000
        else:
            return 0

    def print_results(self):
        """Generates the final CSV output
        Returns a string containing the whole CSV file content to be written to disk
        """
        # import method for binomial test (tip of @Alexey)
        from statsmodels.stats.proportion import proportion_confint

        if not self.cli_params.whitelist:

            # construct header of the CSV output file
            result_string = "circRNA_host_gene\t" \
                            "chr\t" \
                            "start\t" \
                            "stop\t" \
                            "strand\t" \
                            "p-val_circular\t" \
                            "raw_count_circ_rna\t" \
                            "observed_input_peaks_circ_rna\t" \
                            "length_circ_rna\t" \
                            "length_normalized_count_circ_rna\t" \
                            "number_of_features_intersecting_circ\t" \
                            "circ_rna_confidence_interval_0.05\t" \
                            "p-val_linear\t" \
                            "raw_count_host_gene\t" \
                            "observed_input_peaks_host_gene\t" \
                            "length_host_gene_without_circ_rna\t" \
                            "length_normalized_count_host_gene\t" \
                            "number_of_features_intersecting_linear\t" \
                            "host_gene_confidence_interval_0.05\t" \
                            "distance_normalized_counts\n"
        else:
            # construct header of the CSV output file
            result_string = "circRNA_host_gene\t" \
                            "chr\t" \
                            "start\t" \
                            "stop\t" \
                            "strand\t" \
                            "p-val_enriched_circular\t" \
                            "raw_count_enriched_circ_rna\t" \
                            "observed_input_peaks_enriched_circ_rna\t" \
                            "length_enriched_circ_rna\t" \
                            "length_normalized_count_enriched_circ_rna\t" \
                            "number_of_features_intersecting_enriched_circ\t" \
                            "circ_rna_confidence_interval_0.05\t" \
                            "p-val_non_enriched_circRNA\t" \
                            "raw_count_not_enriched_circ_rna\t" \
                            "observed_input_peaks_not_enriched\t" \
                            "NOT_USED\t" \
                            "length_normalized_count_not_enriched\t" \
                            "number_of_features_intersecting__not_enriched\t" \
                            "not_enriched_confidence_interval_0.05\t" \
                            "distance_normalized_counts\n"

        # print(self.observed_counts[0])
        # for all genes we have seen

        with open("f11.txt", 'w') as text_file:
            text_file.write(str(self.observed_counts[0]))

        with open("f22.txt", 'w') as text_file:
            text_file.write(str(self.observed_counts[1]))


        import pprint
        pp1 = pprint.PrettyPrinter(stream=open("f1.txt",'w'),indent=4)
        pp2 = pprint.PrettyPrinter(stream=open("f2.txt",'w'),indent=4)

        pp1.pprint(self.observed_counts[0])
        pp2.pprint(self.observed_counts[1])

        for gene in self.observed_counts[1]:
            # make sure we found a circular RNA
            if gene in self.observed_counts[0]:

                # get the location key of the linear host RNA
                for location_key_linear in self.observed_counts[1][gene]:

                        # for each location key of the circRNA
                        location_key_circular = self.decode_location_key(location_key_linear)["circ_data"]

                        if self.decode_location_key(location_key_circular)["chr"] == \
                                self.decode_location_key(location_key_linear)["chr"] and \
                                self.observed_counts[0][gene][location_key_circular] > 0:

                            # get the count of simulated peaks > than observed peaks
                            count_circular = 0

                            if gene in self.phase_storage and location_key_circular in self.phase_storage[gene][0]:
                                count_circular = self.phase_storage[gene][0][location_key_circular]

                            # get length of the host gene without the circRNA annotation
                            # this is a list: entry 0: circular RNA, entry 1: linear RNA

                            length = []

                            if self.virtual_inclusion_file_path != "all":
                                length.append(self.get_extended_key_data(location_key_circular)["feature_length"])
                                length.append(self.get_extended_key_data(location_key_linear)["feature_length"] -
                                              self.get_extended_key_data(location_key_circular)["feature_length"])
                                location_data_circ = self.decode_location_key(location_key_circular)
                                location_data_linear = self.decode_location_key(location_key_linear)

                            elif self.virtual_inclusion_file_path != "all" and self.cli_params.whitelist:
                                length.append(self.get_extended_key_data(location_key_circular)["feature_length"])
                                length.append(self.get_extended_key_data(location_key_linear)["feature_length"] -
                                              self.get_extended_key_data(location_key_circular)["feature_length"])
                                location_data_circ = self.decode_location_key(location_key_circular)
                                location_data_linear = self.decode_location_key(location_key_linear)
                            else:
                                length = self.linear_length_wo_circ(location_key_circular, location_key_linear)
                                location_data_circ = self.decode_location_key(location_key_circular)
                                location_data_circ["feature_count"] = 0
                                location_data_linear = self.decode_location_key(location_key_linear)
                                location_data_linear["feature_count"] = 0

                            # hint: count means: we saw simulated data that had more peaks
                            # in the region than we observed experimentally

                            count_linear = 0
                            # get the count of simulated peaks > than observed peaks
                            if gene in self.phase_storage and location_key_linear in self.phase_storage[gene][1]:
                                count_linear = self.phase_storage[gene][1][location_key_linear]

                            # get the length-normalized count for the linear RNA
                            count_linear_normalized = self.normalize_count(length[1], count_linear - count_circular)

                            # get the length-normalized count for the circular RNA
                            count_circular_normalized = self.normalize_count(length[0], count_circular)

                            # how many experimental peaks did we see?
                            observed_count_circular = self.observed_counts[0][gene][location_key_circular]
                            observed_count_linear = self.observed_counts[1][gene][location_key_linear]

                            # compute simple p-val
                            p_val_linear = 0
                            p_val_circular = 0

                            if count_linear > 0:
                                p_val_linear = count_linear / self.cli_params.num_iterations

                            if count_circular > 0:
                                p_val_circular = count_circular / self.cli_params.num_iterations

                            # compute a 0.05 confidence interval
                            confidence_interval_circular = proportion_confint(count_circular,
                                                                              self.cli_params.num_iterations,
                                                                              method="beta")

                            confidence_interval_linear = proportion_confint(count_linear,
                                                                            self.cli_params.num_iterations,
                                                                            method="beta")

                            # check that the host gene length is not 0
                            # and that we are above the user-defined threshold
                            # also:
                            # we only want to see entries where the count is lower for the circ RNA

                            # if (length[1] > 0) and p_val_linear <= self.cli_params.pval \
                            #         and observed_count_circular > self.cli_params.threshold\
                            #         and count_circular_normalized < count_linear_normalized:

                            # this distance is a kind of measure how far apart linear and circular RNA are
                            distance = count_linear_normalized - count_circular_normalized

                            # construct the result line
                            result_string += (
                                "%s\t%s\t%d\t%d\t%s\t%f\t%d\t%d\t%d\t%f\t%d\t%s\t%f\t%d\t%d\t%d\t%f\t%d\t%s\t%f\n" %
                                (
                                    gene,
                                    str(location_data_circ["chr"]),
                                    location_data_circ["start"],
                                    location_data_circ["stop"],
                                    location_data_circ["strand"],
                                    p_val_circular,
                                    count_circular,
                                    observed_count_circular,
                                    length[0],
                                    count_circular_normalized,
                                    location_data_circ["feature_count"],
                                    confidence_interval_circular,
                                    p_val_linear,
                                    count_linear,
                                    observed_count_linear,
                                    length[1],
                                    count_linear_normalized,
                                    location_data_linear["feature_count"],
                                    confidence_interval_linear,
                                    distance
                                )
                            )
        # return the data string
        return result_string

    def run_permutation_test(self, current_iteration):
        """Computes the probability that a given number of peaks is significantly higher than the simulated random
        distribution obtained by the random shuffling.
        Returns a dictionary containing the result for a given iteration
        """
        gene_dict = {}

        if current_iteration >= self.cli_params.num_iterations:
            return gene_dict

        self.log_entry("Permutation test iteration %d" % (current_iteration+1))

        # we pre-process the first entry of the list because here we stored the observed counts
        # order of tuples: pos 0: circular rna, pos 1: linear rna

        processed_counts = []

        for rna_type in range(0, 2):
            processed_counts.append(self.process_intersection(self.results[current_iteration][rna_type]))

        for rna_type in range(0, 2):

            for gene, nested_dict in processed_counts[rna_type].items():

                # we need to get the observed count for this gene before we start
                observed_value_dict = self.observed_counts[rna_type][gene]
                # for each location key (for linear that's only one anyway. for circular it may me multiple)
                for location_key, shuffled_value in nested_dict.items():

                    # location_data = self.decode_location_key(location_key)
                    #
                    #
                    # length = self.decode_location_key(location_key)["length"]
                    #
                    #
                    # if self.normalize_count(length, shuffled_value) >= \
                    #         self.normalize_count(length, observed_value_dict[location_key]):

                    # let's test if we observed a higher count in this iteration than web observed experimentally
                    # first make sure the location exists.. should always be true for linear rna but not for
                    # circular RNAs

                    # we need the observed counts from the current circle

                    if rna_type == 1 and gene in self.observed_counts[0]:
                        for location_key_circular in self.observed_counts[0][gene]:
                            tmp_data = self.decode_location_key(location_key_circular)

                            location_key_new = location_key + "_" + \
                                               tmp_data["chr"] + "_" + \
                                               str(tmp_data["start"]) + "_" \
                                               + str(tmp_data["stop"]) + "_" + \
                                               tmp_data["strand"] + "_" + \
                                               str(tmp_data["feature_length"]) + "_" + \
                                               str(tmp_data["feature_count"])
                            circ_count = 0

                            if location_key_circular in processed_counts[0][gene]:
                                circ_count = processed_counts[0][gene][location_key_circular]
                                # for intron mode we have to also subtract / add here

                                # we got the first circRNA of the buddy pair, we now need to get the second one
                                count_correction = 0
                                observed_count_correction = 0

                                if location_key_circular in self.circRNA_buddies:

                                    if self.circRNA_buddies[location_key_circular] in processed_counts[0][gene]:

                                        count_correction = processed_counts[0][gene][self.circRNA_buddies[location_key_circular]]

                                    if self.circRNA_buddies[location_key_circular] in self.observed_counts[0][gene]:

                                        observed_count_correction = self.observed_counts[0][gene][self.circRNA_buddies[location_key_circular]]

                            if shuffled_value-circ_count-count_correction > observed_value_dict[location_key_new] - observed_count_correction:

                                # Yes, it's higher, so we update the count of "more than observed" for this gene
                                if gene not in gene_dict:
                                    # initialize new dict entry
                                    gene_dict[gene] = {}

                                # look if we already have circ/linear rna entries
                                if rna_type not in gene_dict[gene]:
                                    # first time we see a higher shuffled value
                                    gene_dict[gene][rna_type] = {}

                                if location_key_new not in gene_dict[gene][rna_type]:
                                    gene_dict[gene][rna_type][location_key_new] = True

                    else:
                        # if circ rna get sister rna and get counts for obs and raw

                        count_correction = 0
                        observed_count_correction = 0

                        if location_key in self.circRNA_buddies:

                            if self.circRNA_buddies[location_key] in processed_counts[0][gene]:
                                count_correction = processed_counts[0][gene][self.circRNA_buddies[location_key]]
                                observed_count_correction = self.observed_counts[0][gene][self.circRNA_buddies[location_key]]

                        if shuffled_value + count_correction > observed_value_dict[location_key] + observed_count_correction:

                            # Yes, it's higher, so we update the count of "more than observed" for this gene
                            if gene not in gene_dict:
                                # initialize new dict entry
                                gene_dict[gene] = {}

                            # look if we already have circ/linear rna entries
                            if rna_type not in gene_dict[gene]:
                                # first time we see a higher shuffled value
                                gene_dict[gene][rna_type] = {}

                            if location_key not in gene_dict[gene][rna_type]:
                                gene_dict[gene][rna_type][location_key] = True

        return gene_dict

    def random_sample_step(self, iteration, circ_rna_bed, annotation_bed, shuffled_peaks):
        """Supervises the intersection of the shuffled peaks of a given iteration
        with circular and linear RNA annotation
        Returns a list of bedtools intersection outputs: circular [0] and linear RNA [1] intersection counts
        """

        if not self.cli_params.whitelist:
            # get circular and linear intersect
            circular_intersect = self.do_intersection(shuffled_peaks[iteration], circ_rna_bed)
            linear_intersect = self.do_intersection(shuffled_peaks[iteration], annotation_bed)
        else:
            # get circular and circular intersect for enriched / non enriched exons

            # the enriched exons
            circular_intersect = self.do_intersection(shuffled_peaks[iteration], self.whitelist_fg)

            # the "normal" exons are running as linear from here on
            linear_intersect = self.do_intersection(shuffled_peaks[iteration], circ_rna_bed)

        # process results of the intersects
        intersects = [circular_intersect, linear_intersect]

        self.log_entry("Processed intersections for iteration %s" % (iteration+1))

        return intersects

    def clean_up_temp_files(self):
        """Delete temporary files created by pybedtools
        """
        self.log_entry("Cleaning up temporary files")

        import glob
        for tmp_file in glob.glob(self.cli_params.tmp_directory+"/"+"pybedtools*"):
            os.remove(tmp_file)
            self.log_entry("Deleting " + tmp_file)

        try:
            os.rmdir(self.cli_params.tmp_directory)
        except OSError:
            pass

        self.log_entry("Done")

    def module_name(self):
        """"Return a string representing the name of the module."""
        return self.program_name
