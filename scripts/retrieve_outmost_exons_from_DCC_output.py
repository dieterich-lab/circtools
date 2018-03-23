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


def generate(input_file, exons_bed):

    exons = pybedtools.example_bedtool(exons_bed)

    with open(input_file) as fp:

            for line in fp:

                # make sure we remove the header
                if line.startswith('Chr'):
                    continue


                line = line.rstrip()
                current_line = line.split('\t')
                print("processing " + current_line[3])

                sep = "\t"
                bed_string = sep.join([current_line[0],
                                      current_line[1],
                                      current_line[2],
                                      current_line[3],
                                      str(0),
                                      current_line[5]])

                virtual_bed_file = pybedtools.BedTool(bed_string, from_string=True)
                result = exons.intersect(virtual_bed_file, s=True)
                #print(result)

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

group.add_argument("-t",
                   "--threshold",
                   dest="base_threshold",
                   help="max length of intron sequence",
                   type=int,
                   default=2000
                   )

args = parser.parse_args()

dcc_input = generate(args.dcc_file, args.gtf_file)
