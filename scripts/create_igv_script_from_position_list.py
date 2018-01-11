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


def build_tracks(cli_args):
    data = parse_file(cli_args.input_file)
    generate_igv_script(data, cli_args)


def generate_alternative_exon_tracks(exon_directory_list):

    # note: BED graph tracks have to be printed before all BED tracks otherwise IGV does not show them

    # BED graph tracks
    for directory in exon_directory_list:
        print("load " + directory + "exon_analysis_exon_fc_track.bedgraph")
        print("load " + directory + "exon_analysis_exon_pval_track.bedgraph")

    # BED tracks
    for directory in exon_directory_list:
        print("load " + directory + "exon_analysis_dcc_bsj_enriched_track.bed")

    # DCC track (use first list element (0) because it's the same track for all runs)
    print("load " + exon_directory_list[0] + "exon_analysis_dcc_predictions_track.bed")


def generate_raw_data_tracks(bam_file_list):

    # BAM tracks
    for bam_file in bam_file_list:
        print("load " + bam_file)


def generate_reconstruct_tracks(reconstruct_file_list):

    # BED tracks
    for bed_file in reconstruct_file_list:
        print("load " + bed_file)


def generate_header(gene_name, genome_build, output_directory):

    print("new")
    print("genome " + genome_build)
    print("snapshotDirectory " + output_directory)
    print("maxPanelHeight 50000")
    print("load /mnt/big_data/genomes/GRCh38_85/GRCh38.85.sorted.gtf")


def generate_footer(bam_file_list):

    print("#####################")


def location_zoom(location, zoom_level):
    location_bits = location.split(':')
    starts_end_bit = location_bits[1].split('-')
    chromosome = location_bits[0]
    start_bit = str(int(starts_end_bit[0]) - zoom_level)
    stop_bit = str(int(starts_end_bit[1]) + zoom_level)
    return chromosome + ":" + start_bit + "-" + stop_bit


def generate_igv_script(parsed_data, cli_args):

    max_genes_to_show = 0
    for gene_name in parsed_data:
        for location in parsed_data[gene_name]['location'].keys():

            generate_header(gene_name, cli_args.genome_build, cli_args.output_directory)

            generate_raw_data_tracks(cli_args.bam_files)
            generate_alternative_exon_tracks(cli_args.alt_exon_dirs)
            generate_reconstruct_tracks(cli_args.fuchs_files)

            generate_footer(cli_args.bam_files)

            print("region " + location)
            print("goto " + location_zoom(location, 100000))
            print("squish")

            print("snapshot " + str(parsed_data[gene_name]['rank']) + "_" + gene_name + "_" + location + "_gene.png")
            print("snapshot " + str(parsed_data[gene_name]['rank']) + "_" + gene_name + "_" + location + "_gene.svg")

            print("goto " + location_zoom(location, 5000))
            print("snapshot " + str(parsed_data[gene_name]['rank']) + "_" + gene_name + "_" + location + "_zoom.png")
            print("snapshot " + str(parsed_data[gene_name]['rank']) + "_" + gene_name + "_" + location + "_gene.svg")

        # stop if we reached the maximum number of genes to print
        max_genes_to_show += 1
        if max_genes_to_show > cli_args.max_genes:
            break


def parse_file(input_file):

    from collections import OrderedDict
    entries_list = OrderedDict()

    with open(input_file) as fp:
        rank = 1
        for line in fp:
            current_line = line.split('\t')
            host_gene = current_line[0]
            location = "chr" + current_line[1] + ":" + current_line[2] + "-" + current_line[3]

            # create key
            if host_gene not in entries_list:
                entries_list[host_gene] = {}
                entries_list[host_gene]['rank'] = rank
                entries_list[host_gene]['location'] = OrderedDict()
                entries_list[host_gene]['location'][location] = 1
                rank += 1
            else:
                if location not in entries_list[host_gene]['location']:
                    entries_list[host_gene]['location'][location] = 1

    return entries_list


# main script starts here

parser = argparse.ArgumentParser(description='Create an auto-executing IGV script')

group = parser.add_argument_group("Input")

group.add_argument("-i",
                   "--input-file",
                   dest="input_file",
                   help="A CSV/BED file from DCC / CircTest",
                   required=True
                   )

group.add_argument("-b",
                   "--bam-files",
                   dest="bam_files",
                   nargs='*',
                   help="List of one or more BAM files with read mapping data",
                   required=True
                   )

group.add_argument("-a",
                   "--alternative-exon-dirs",
                   dest="alt_exon_dirs",
                   nargs='*',
                   help="List of one or more directories containing the result files "
                        "of the exon module",
                   required=True
                   )

group.add_argument("-f",
                   "--fuchs-files",
                   dest="fuchs_files",
                   nargs='*',
                   help="List of one or more BED files containing results of "
                        "FUCHS / reconstruct module",
                   required=True
                   )

group.add_argument("-g",
                   "--genome",
                   dest="genome_build",
                   help="Which genome build to use as reference [Default: hg38]",
                   choices=("hg38", "hg19", "mm9", "mm10"),
                   default="hg38"
                   )

group.add_argument("-m",
                   "--max-number",
                   dest="max_genes",
                   help="Maximum number of genes to show [Default: 5]",
                   type=int,
                   default=5
                   )

group = parser.add_argument_group("Output")

group.add_argument("-o",
                   "--output-directory",
                   dest="output_directory",
                   help="Directory for snapshots created by IGV [Default: ./]",
                   default="./"
                   )

args = parser.parse_args()

build_tracks(args)
