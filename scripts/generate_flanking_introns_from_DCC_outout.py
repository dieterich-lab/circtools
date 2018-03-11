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


def return_bed_line(five_prime, three_prime, gene_name, comment, threshold):

    three = three_prime.split('_')
    five = five_prime.split('_')

    three[1] = int(three[1])
    three[2] = int(three[2])

    five[1] = int(five[1])
    five[2] = int(five[2])

    if three[2] > three[1] + threshold:
        three[2] = three[1] + threshold

    if five[1] < five[2] - threshold:
        five[1] = five[2] - threshold

    start = five[1]
    stop = three[2]

    return five[0] + "\t" +\
        str(start) + "\t" +\
        str(stop) + "\t" +\
        gene_name + "\t" +\
        "0" + "\t" + \
        five[3] + "\t" +\
        "none" + "\t" +\
        comment


def extract_start(string):

    split = string.split('_')
    return split[0] + "_" + split[1] + "_" + split[3]


def return_wobble(string, wobble):

    split = string.split('_')
    return split[0] + "_" + str(int(split[1])+wobble) + "_" + split[2]


def extract_stop(string):

    split = string.split('_')
    return split[0] + "_" + split[2] + "_" + split[3]


def print_results(dcc_dict, gtf_start_dict, gtf_stop_dict, threshold):

    for entry in dcc_dict:

        stop = extract_stop(entry)
        start = extract_start(entry)

        start_found = ""
        stop_found = ""

        for wobble in range(-200, 200):

            current = return_wobble(stop, wobble)

            if current in gtf_stop_dict:
                # stop_found = (stop + " STOP " + gtf_stop_dict[current])
                stop_found = gtf_stop_dict[current]

            current = return_wobble(start, wobble)

            if current in gtf_start_dict:
                # start_found = (start + " START " + gtf_start_dict[current])
                start_found = gtf_start_dict[current]

        if start_found and stop_found:
            #  print(return_bed_line(entry, dcc_dict[entry]))

            # line for downstream intron
            print(return_bed_line(start_found, stop_found, dcc_dict[entry], entry, threshold))

            # line for upstream exon
            #print(return_bed_line(stop_found, dcc_dict[entry], entry, 1, threshold))

    return


def parse_dcc_file(input_file):
    from collections import OrderedDict
    loc_list = OrderedDict()

    with open(input_file) as fp:

        for line in fp:

            # make sure we remove the header
            if line.startswith('Chr'):
                continue

            current_line = line.split('\t')
            loc = current_line[0] + "_" + \
                current_line[1] + "_" + \
                current_line[2] + "_" + \
                current_line[5]

            if loc not in loc_list:
                loc_list[loc] = current_line[3]

    return loc_list


def parse_gtf_file(input_file):
    from collections import OrderedDict
    start_list = OrderedDict()
    stop_list = OrderedDict()

    with open(input_file) as fp:

        for line in fp:
            current_line = line.split('\t')

            # stop point of previous exon
            stop = current_line[0] + "_" + str(int(current_line[3]) - 1) + "_" + current_line[6]

            # start point of next exon
            start = current_line[0] + "_" + str(int(current_line[4]) + 1) + "_" + current_line[6]

            loc = current_line[0] + "_" + \
                  current_line[3] + "_" + \
                  current_line[4] + "_" + \
                  current_line[6]  # stop point of previous exon

            # create key
            stop_list[stop] = loc

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

group.add_argument("-t",
                   "--threshold",
                   dest="base_threshold",
                   help="max length of intron sequence",
                   type=int,
                   default=2000
                   )

args = parser.parse_args()

dcc_input = parse_dcc_file(args.dcc_file)

gtf_input = parse_gtf_file(args.gtf_file)

print_results(dcc_input, gtf_input[0], gtf_input[1], args.base_threshold)


