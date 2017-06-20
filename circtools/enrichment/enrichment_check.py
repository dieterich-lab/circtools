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
        self.observed_counts = []
        self.results = []
        self.phase_storage = {}

    def run_module(self):

        # set time format
        time_format = time.strftime("%Y_%m_%d__%H_%M")

        # set up the multiprocessing pool for multi-threading
        mp_pool = multiprocessing.Pool(processes=self.cli_params.num_processes)

        # let's first check if the temporary directory exists
        if not os.path.exists(self.cli_params.tmp_directory):
            os.makedirs(self.cli_params.tmp_directory)

        # let's first check if the temporary directory exists
        if not os.path.exists(self.cli_params.output_directory):
            os.makedirs(self.cli_params.output_directory)

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
                                                  self.cli_params.output_filename + "_" + time_format + ".log"),
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

        gene_annotation_file = self.cli_params.output_directory +\
                               '/' + self.cli_params.output_filename +\
                               "_" + \
                               os.path.basename(self.cli_params.annotation) +\
                               '_genes.bed'

        annotation_bed.saveas(gene_annotation_file)

        # read in circular RNA predictions from DCC
        circ_rna_bed = self.read_circ_rna_file(self.cli_params.circ_rna_input,
                                               annotation_bed,
                                               self.cli_params.has_header)

        # do circle saves
        circle_annotation_file = self.cli_params.output_directory +\
                                 '/' + self.cli_params.output_filename +\
                                 "_" +\
                                 os.path.basename(self.cli_params.circ_rna_input) +\
                                 '_circles.bed'

        circ_rna_bed.saveas(circle_annotation_file)

        # create list of shuffled peaks
        shuffled_peaks = (mp_pool.map(functools.partial(self.shuffle_peaks_through_genome,
                                                       bed_file=supplied_bed,
                                                       genome_file=self.cli_params.genome_file),
                                     range(self.cli_params.num_iterations)))

        shuffled_peaks.insert(0, supplied_bed)

        # shuffled_peaks.append(supplied_bed)
        # do the intersections
        self.results = mp_pool.map(functools.partial(self.random_sample_step,
                                                circ_rna_bed=circ_rna_bed,
                                                annotation_bed=annotation_bed,
                                                shuffled_peaks=shuffled_peaks,
                                                ), range(self.cli_params.num_iterations+1))

        self.observed_counts = (
            self.process_intersection(self.results[0][0]),
            self.process_intersection(self.results[0][1])
        )

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

            # results of one phase of the computation
            intermediate_result = mp_pool.map(functools.partial(self.run_permutation_test),
                                              range(phase*iterations_per_phase+1, (phase+1)*iterations_per_phase)
                                              )

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

            # clean up here is important to keep memory down
            # gc.collect()

        # generate the result table
        result_table = self.print_results()

        # and print it to a file
        result_file = self.cli_params.output_directory + "/" +\
                      self.cli_params.output_filename + "_" +\
                      str(self.cli_params.num_iterations) + "_" +\
                      time_format + ".csv"

        with open(result_file, 'w') as text_file:
            text_file.write(result_table)

        self.clean_up_temp_files()

        # ------------------------------------- Function definitions start here ---------------------------------------

    def read_circ_rna_file(self, circ_rna_input, annotation_bed, has_header):
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
                # skip first line with the header
                # we assume it's there (DCC default)
                if has_header:
                    next(line_iterator)
                bed_content = ""
                bed_entries = 0
                bed_peak_sizes = 0
                for line in line_iterator:
                    columns = line.split('\t')

                    # extract chromosome, start, stop, gene name, and strand
                    entry = [self.strip_chr_name(columns[0]), columns[1], columns[2], columns[3], "0", columns[5]]

                    # concatenate lines to one string
                    bed_content += '\t'.join(entry) + "\n"

                    bed_entries += 1
                    bed_peak_sizes += (int(columns[2]) - int(columns[1]))

            # create a "virtual" BED file
            virtual_bed_file = pybedtools.BedTool(bed_content, from_string=True)
            # Todo: figure out what this code was supposed to do
            test = annotation_bed.intersect(virtual_bed_file, s=True)

        self.log_entry("Done parsing circular RNA input file:")
        self.log_entry("=> %s circular RNAs, %s nt average (theoretical unspliced) length" %
                       (bed_entries, round(bed_peak_sizes / bed_entries)))

        return test

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

    def read_annotation_file(self, annotation_file, entity="gene"):
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
            with file_handle:
                line_iterator = iter(file_handle)
                bed_content = ""
                entity_count = 1
                for line in line_iterator:
                    # we skip any comment lines
                    if line.startswith("#"):
                        continue

                    # split up the annotation line
                    columns = line.split('\t')

                    # we only want the coordinates of the gene entries
                    if not (columns[2] == entity):
                        continue

                    # columns 8 is the long annotation string from GTF
                    gene_name = self.extract_gene_name_from_gtf(columns[8])

                    # extract chromosome, start, stop, score(0), name and strand
                    # we hopefully have a gene name now and use this one for the entry
                    entry = [self.strip_chr_name(columns[0]), columns[3], columns[4], gene_name, str(0), columns[6], ]

                    # concatenate lines to one string
                    bed_content += '\t'.join(entry) + "\n"

                    # sys.stdout.write("Processing %s # %s \r" % (entity, entity_count))
                    # sys.stdout.flush()

                    entity_count += 1

                # "escape the \r from counting output"
                sys.stdout.write("\n")

                # count will be increased one more time even if done - so we subtract 1
                self.log_entry("Processed %s entries" % (entity_count - 1))

            # create a "virtual" BED file
            virtual_bed_file = pybedtools.BedTool(bed_content, from_string=True)

        self.log_entry("Done parsing annotation")

        return virtual_bed_file

    def shuffle_peaks_through_genome(self, iteration, bed_file, genome_file):
        """Gets a (virtual) BED files and shuffle its contents throughout the supplied genome
        Will only use supplied annotation for features (in our case only transcript regions)
        """
        # set temporary directory for pybedtools
        pybedtools.set_tempdir(self.cli_params.tmp_directory)

        self.log_entry("Processing shuffling thread %d" % (iteration+1))
        shuffled_bed = bed_file.shuffle(g=genome_file)

        return shuffled_bed

    def do_intersection(self, query_bed, base_bed):
        """Gets two bed files (supplied peaks and circle coordinates) and does an intersection
        """
        # set temporary directory for pybedtools
        pybedtools.set_tempdir(self.cli_params.tmp_directory)

        # we employ the c=true parameter to directly get the counts as part of the results
        intersect_return = base_bed.intersect(query_bed, c=True)
        return intersect_return

    @staticmethod
    def process_intersection(intersection_input, normalize=False):
        """Gets two bed files (supplied peaks and circle coordinates) and does an intersection
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

            # we have to create the nested dictionaries if not already existing
            if bed_feature[3] not in count_table:
                count_table[bed_feature[3]] = {}

            if key not in count_table[bed_feature[3]]:
                count_table[bed_feature[3]][key] = {}

            # set the appropriate dict entry
            count_table[bed_feature[3]][key] = normalize_count(bed_feature[1],
                                                               bed_feature[2],
                                                               int(bed_feature[6])
                                                               )
        # return one unified count table
        return count_table

    @staticmethod
    def decode_location_key(key):
        tmp = key.split("_")
        data = {"chr": tmp[0], "start": int(tmp[1]), "stop": int(tmp[2]), "strand": tmp[2]}
        return data

    def linear_length_wo_circ(self, key_circ, key_linear):
        circ = self.decode_location_key(key_circ)
        linear = self.decode_location_key(key_linear)
        length_circ = (circ["stop"] - circ["start"])
        length_lin = (linear["stop"] - linear["start"]) - length_circ
        return length_circ, length_lin

    @staticmethod
    def normalize_count(length, count):
        if count > 0 and length > 0:
            return (count / length) * 100000
        else:
            return 0

    def print_results(self):

        # import method for binomial test (tip of @Alexey)
        from statsmodels.stats.proportion import proportion_confint

        # construct header of the CSV output file
        result_string = "gene\t" \
                        "location\t" \
                        "p-val_circular\t" \
                        "raw_count_circ_rna\t" \
                        "observed_input_peaks_circ_rna\t" \
                        "length_circ_rna\t" \
                        "length_normalized_count_circ_rna\t" \
                        "circ_rna_confidence_interval_0.05\t" \
                        "p-val_linear\t" \
                        "raw_count_host_gene\t" \
                        "observed_input_peaks_host_gene\t" \
                        "length_host_gene_without_circ_rna\t" \
                        "length_normalized_count_host_gene\t" \
                        "host_gene_confidence_interval_0.05\t" \
                        "distance_normalized_counts\n"

        # for all genes we have seen
        for gene in self.observed_counts[1]:
            # make sure we found a circular RNA
            if gene in self.observed_counts[0]:

                # get the location key of the linear host RNA
                for location_key_linear in self.observed_counts[1][gene]:

                        # for each location key of the circRNA
                        for location_key_circular in self.observed_counts[0][gene]:

                            if self.observed_counts[0][gene][location_key_circular] > 0:

                                # get the count of simulated peaks > than observed peaks
                                count_circular = 0

                                if gene in self.phase_storage and location_key_circular in self.phase_storage[gene][0]:
                                    count_circular = self.phase_storage[gene][0][location_key_circular]

                                # get length of the host gene without the circRNA annotation
                                # this is a list: entry 0: circular RNA, entry 1: linear RNA
                                length = self.linear_length_wo_circ(location_key_circular, location_key_linear)

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
                                    "%s\t%s\t%f\t%d\t%d\t%d\t%f\t%s\t%f\t%d\t%d\t%d\t%f\t%s\t%f\n" %
                                                (gene,
                                                 location_key_circular,
                                                 p_val_circular,
                                                 count_circular,
                                                 observed_count_circular,
                                                 length[0],
                                                 count_circular_normalized,
                                                 confidence_interval_circular,
                                                 p_val_linear,
                                                 count_linear,
                                                 observed_count_linear,
                                                 length[1],
                                                 count_linear_normalized,
                                                 confidence_interval_linear,
                                                 distance
                                                 )
                                )
        # return the data string
        return result_string

    def run_permutation_test(self, current_iteration):

        gene_dict = {}

        if current_iteration >= self.cli_params.num_iterations:
            return gene_dict

        self.log_entry("Permutation test iteration %d" % (current_iteration+1))

        # we pre-process the first entry of the list because here we stored the observed counts
        # order of tuples: pos 0: circular rna, pos 1: linear rna

        # iterate through circular and linear intersection
        for rna_type in range(0, 2):

            processed_counts = self.process_intersection(self.results[current_iteration][rna_type])

            for gene, nested_dict in processed_counts.items():

                # we need to get the observed count for this gene before we start
                observed_value_dict = self.observed_counts[rna_type][gene]
                # for each location key (for linear that's only one anyway. for circular it may me multiple)
                for location_key, shuffled_value in nested_dict.items():

                    # let's test if we observed a higher count in this iteration than web observed experimentally
                    # first make sure the location exists.. should always be true for linear rna but not for
                    # circular RNAs
                    if shuffled_value >= observed_value_dict[location_key]:

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
        """Logs to log file and prints on screen
        """

        # get circular and linear intersect
        circular_intersect = self.do_intersection(shuffled_peaks[iteration], circ_rna_bed)
        linear_intersect = self.do_intersection(shuffled_peaks[iteration], annotation_bed)

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
            print(tmp_file)

        os.rmdir(self.cli_params.tmp_directory)

        self.log_entry("Done")

    def module_name(self):
        """"Return a string representing the name of the module."""
        return self.program_name
