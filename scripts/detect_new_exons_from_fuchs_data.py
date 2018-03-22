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


def parse_fuchs_file(input_file, exon_dict):

    with open(input_file) as fp:

        for line in fp:

            if line.startswith("#"):
                continue

            line = line.rstrip()

            columns = line.split('\t')

            start = 0
            stop = 0

            for wobble in range(-10, 10):

                if columns[0] + "_" + str(int(columns[1])+wobble) in exon_dict:
                    start = 1

                if columns[0] + "_" + str(int(columns[2])+wobble) in exon_dict:
                    stop = 1

            if start == 0 and stop == 0:
                print(line)


def parse_gtf_file(input_file):
    from collections import OrderedDict
    import sys

    exon_dict = OrderedDict()

    try:
        file_handle = open(input_file)
    except PermissionError:
        message = ("Input file " + str(input_file) + " cannot be read, exiting.")
        sys.exit(message)
    else:

        with file_handle:
            line_iterator = iter(file_handle)
            for line in line_iterator:
                # we skip any comment lines
                if line.startswith("#"):
                    continue

                # split up the annotation line
                columns = line.split('\t')

                # we only want the coordinates of the gene entries
                if not (columns[2] == "exon"):
                    continue

                # we do not want any 0-length intervals -> bedtools segfault
                if int(columns[4]) - int(columns[3]) == 0:
                    continue

                start_key = str(columns[0])+"_"+str(columns[3])
                stop_key = str(columns[0])+"_"+str(columns[4])

                exon_dict[start_key] = 1
                exon_dict[stop_key] = 1

    return exon_dict

# main script starts here


parser = argparse.ArgumentParser(description='Create ')

group = parser.add_argument_group("Input")

group.add_argument("-g",
                   "--base-exons",
                   dest="base_exon_file",
                   help="Bed file holding the known base exons",
                   required=True
                   )

group.add_argument("-f",
                   "--fuchs-exons",
                   dest="fuchs_exon_file",
                   help="Bed file holding the reconstructed FUCHS exons",
                   required=True
                   )


args = parser.parse_args()

gtf_input = parse_gtf_file(args.base_exon_file)
parse_fuchs_file(args.fuchs_exon_file, gtf_input)
