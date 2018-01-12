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


def extract_start(string):

    split = string.split('_')
    return split[0] + "_" + split[1] + "_" + split[3]


def extract_stop(string):

    split = string.split('_')
    return split[0] + "_" + split[2] + "_" + split[3]


def print_results(dcc_dict, gtf_start_dict, gtf_stop_dict):

    for entry in dcc_dict:

        stop = extract_stop(entry)
        start = extract_start(entry)

        start_found = ""
        stop_found = ""

        #print(entry)
        if stop in gtf_stop_dict:
            stop_found = (stop + " STOP " + gtf_stop_dict[stop])

        if start in gtf_start_dict:
            start_found = (start + " START " + gtf_start_dict[start])

        if start_found and stop_found:
            print(entry)

    return


def parse_dcc_file(input_file):
    from collections import OrderedDict
    loc_list = OrderedDict()

    with open(input_file) as fp:

        for line in fp:
            current_line = line.split('\t')
            loc = current_line[0] + "_" + \
                  current_line[1] + "_" + \
                  current_line[2] + "_" + \
                  current_line[5]  # stop point of previous exon

            if loc not in loc_list:
                loc_list[loc] = 1

    return loc_list


def parse_gtf_file(input_file):
    from collections import OrderedDict
    start_list = OrderedDict()
    stop_list = OrderedDict()

    with open(input_file) as fp:

        for line in fp:
            current_line = line.split('\t')
            stop = current_line[0] + "_" + str(int(current_line[3]) - 1) + "_" + current_line[6]  # stop point of previous exon
            start = current_line[0] + "_" + str(int(current_line[4]) + 1) + "_" + current_line[6]  # start point of next exon

            loc = current_line[0] + "_" + \
                  current_line[3] + "_" + \
                  current_line[4] + "_" + \
                  current_line[6]  # stop point of previous exon

            # create key
            if stop not in stop_list:
                stop_list[stop] = loc

            if start not in start_list:
                start_list[start] = loc

    return start_list, stop_list


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
                   help="GTF file with all introns",
                   required=True
                   )

args = parser.parse_args()

dcc_input = parse_dcc_file(args.dcc_file)

gtf_input = parse_gtf_file(args.gtf_file)

print_results(dcc_input, gtf_input[0], gtf_input[1])


