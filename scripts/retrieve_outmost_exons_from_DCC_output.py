#!/usr/bin/env python3

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

import argparse
import pybedtools
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import os
import signal


# Register an handler for the timeout
def handler(signum, frame):
    raise Exception("Maximum execution time for remote BLAST reached. Please try again later.")


def call_blast(input_file):
    return_handle = NCBIWWW.qblast("blastn",
                                   "GPIPE/9606/current/rna",
                                   input_file,
                                   hitlist_size=10,
                                   expect=1000,
                                   word_size=7,
                                   gapcosts="5 2"
                                   )
    return return_handle


# Register the signal function handler
signal.signal(signal.SIGALRM, handler)

# Define a timeout for your function
signal.alarm(120)


def generate(input_file, exons_bed, fasta_file, tmp_file):

    exons = pybedtools.example_bedtool(exons_bed)

    line_number = 0

    open(tmp_file, 'w').close()

    with open(input_file) as fp:

            for line in fp:

                # make sure we remove the header
                if line.startswith('Chr'):
                    continue

                line_number += 1

                line = line.rstrip()
                current_line = line.split('\t')

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

                    if bed_feature[1] == current_line[1] and start == 0:
                        fasta_bed_line_start += result_line + "\n"
                        start = 1

                    if bed_feature[2] == current_line[2] and stop == 0:
                        fasta_bed_line_stop += result_line + "\n"
                        stop = 1

                virtual_bed_file_start = pybedtools.BedTool(fasta_bed_line_start, from_string=True)
                virtual_bed_file_stop = pybedtools.BedTool(fasta_bed_line_stop, from_string=True)

                fasta = pybedtools.example_filename(fasta_file)

                virtual_bed_file_start = virtual_bed_file_start.sequence(fi=fasta)
                virtual_bed_file_stop = virtual_bed_file_stop.sequence(fi=fasta)

                if stop == 0 or start == 0:

                    continue

                else:
                    exon1 = open(virtual_bed_file_start.seqfn).read().split("\n", 1)[1].rstrip()
                    exon2 = open(virtual_bed_file_stop.seqfn).read().split("\n", 1)[1].rstrip()
                    sep = "_"
                    name = sep.join([current_line[3],
                                     current_line[0],
                                     current_line[1],
                                     current_line[2],
                                     current_line[5]])
                    print("extracting flanking exons for circRNA #", line_number, name, end="\n", flush=True)

                    with open(tmp_file, 'a') as data_store:
                        data_store.write("\t".join([name, exon1, exon2, "\n"]))
                # TODO: remove this constraint
                if line_number == 30:
                    break

    # need to define path top R wrapper
    primer_script = 'circtools_primex_wrapper.R'

    # ------------------------------------ run script and check output -----------------------

    script_result = os.popen(primer_script + " " + tmp_file).read()

    # this is the first time we look through the input file
    # we collect the primer sequences and unify everything in one blast query

    blast_object_cache = {}
    blast_result_cache = {}

    blast_input_file = ""

    for line in script_result.splitlines():
        entry = line.split('\t')
        circular_rna_id = entry[0].split('_')

        print(entry[0])

        # both primers are already in cache
        if entry[1] in blast_object_cache and entry[2] in blast_object_cache:
            print("do not blast")
        # only blast 1
        elif entry[2] in blast_object_cache:
            blast_input_file += "\n>" + entry[1] + "\n" + entry[1]
            blast_object_cache[entry[1]] = 1
        # only blast 2
        elif entry[1] in blast_object_cache:
            blast_input_file += "\n>" + entry[2] + "\n" + entry[2]
            blast_object_cache[entry[2]] = 1
        # nothing seen yet, blast both
        else:
            blast_input_file += "\n>" + entry[1] + "\n"+entry[1]+"\n>" + entry[2] + "\n"+entry[2]
            blast_object_cache[entry[1]] = 1
            blast_object_cache[entry[2]] = 1

    # check if we have to blast

    if blast_input_file:

        try:
            print("Sending " + str(len(blast_object_cache)) + " primers to BLAST")
            result_handle = call_blast(blast_input_file)
        except Exception as exc:
            print(exc)
            exit(-1)

        with open("my_blast.xml", "w") as out_handle:
            out_handle.write(result_handle.read())

        result_handle.close()
        result_handle = open("my_blast.xml")

        blast_records = NCBIXML.parse(result_handle)

        for blast_record in blast_records:

            blast_object_cache[blast_record.query] = blast_record

    # we have our two primer blast results now in this object, either fresh or cached
    blast_obj = [blast_object_cache[entry[1]], blast_object_cache[entry[2]]]

    for blast_record in blast_obj:

        for description in blast_record.descriptions:

            # filter out the host gene we're in now
            # also filter out all "PREDICTED" stuff
            if description.title.find(circular_rna_id[0]) == -1:
                if blast_record.query not in blast_result_cache:
                    blast_result_cache[blast_record.query] = []
                blast_result_cache[blast_record.query].append(description.title)

    # go again through the primex output
    # this time we'll add the blast results as last column

    primex_data_with_blast_results = ""

    for line in script_result.splitlines():
        entry = line.split('\t')

        # split up the identifier for final plotting
        line = line.replace("_", "\t")

        left_result = "No hits"
        right_result = "No hits"

        if entry[1] in blast_result_cache:
            left_result = ";".join(blast_result_cache[entry[1]])

        if entry[2] in blast_result_cache:
            right_result = ";".join(blast_result_cache[entry[2]])

        # update line
        primex_data_with_blast_results += line + "\t" + left_result + "\t" + right_result + "\n"

    # print(primex_data_with_blast_results)

    with open(tmp_file, 'w') as data_store:
        data_store.write(primex_data_with_blast_results)

    # need to define path top R wrapper
    primer_script = 'circtools_primex_formatter.R'

    # ------------------------------------ run script and check output -----------------------

    primex_data_formatted = os.popen(primer_script + " " + tmp_file).read()

    with open("/tmp/bla.html", 'w') as data_store:
        data_store.write(primex_data_formatted)


# main script starts here


parser = argparse.ArgumentParser(description='Create ')

group = parser.add_argument_group("Input")

group.add_argument("-d",
                   "--dcc-file",
                   dest="dcc_file",
                   help="CircCoordinates file from DCC",
                   required=True
                   )

group.add_argument("-g",
                   "--gtf-file",
                   dest="gtf_file",
                   help="GTF file with all exons",
                   required=True
                   )

group.add_argument("-f",
                   "--fasta",
                   dest="fasta_file",
                   help="FASTA file with genome sequence",
                   required=True
                   )

group.add_argument("-t",
                   "--threshold",
                   dest="base_threshold",
                   help="max length of intron sequence",
                   type=int,
                   default=2000
                   )

args = parser.parse_args()

generate(args.dcc_file, args.gtf_file, args.fasta_file, "/tmp/circtools_primex.tmp")
