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
import sys
import os.path

parser = argparse.ArgumentParser(description='Exon composition script')

parser.add_argument("-f1",
                    "--file-list-1",
                    dest="file_list1",
                    help="The files to run on",
                    required=True)

args = parser.parse_args()

def print_tracks():




    print("load /home/tjakobi/work/projects/circRNA/encode_paper/circtools_exon/k562/exon_analysis_exon_fc_track.bedgraph")
    print("load /home/tjakobi/work/projects/circRNA/encode_paper/circtools_exon/k562/exon_analysis_exon_pval_track.bedgraph")

    print("load /home/tjakobi/work/projects/circRNA/encode_paper/circtools_exon/hepg2/exon_analysis_exon_fc_track.bedgraph")
    print("load /home/tjakobi/work/projects/circRNA/encode_paper/circtools_exon/hepg2/exon_analysis_exon_pval_track.bedgraph")

    print("load /home/tjakobi/work/projects/circRNA/encode_paper/circtools_exon/hepg2/exon_analysis_dcc_bsj_enriched_track.bed")
    print("load /home/tjakobi/work/projects/circRNA/encode_paper/circtools_exon/k562/exon_analysis_dcc_bsj_enriched_track.bed")

    print("load /home/tjakobi/work/projects/circRNA/encode_paper/circtools_exon/k562/exon_analysis_dcc_predictions_track.bed")

    print("load /home/tjakobi/work/projects/circRNA/encode_paper/fuchs/K562.bam")
    print("load /home/tjakobi/work/projects/circRNA/encode_paper/fuchs/HepG2.bam")
    print("load /home/tjakobi/work/projects/circRNA/encode_paper/fuchs/output/igv/B_exon_chain_12.bed")
    print("load /home/tjakobi/work/projects/circRNA/encode_paper/fuchs/output/igv/D_exon_chain_12.bed")
    print("load /home/tjakobi/work/projects/circRNA/encode_paper/fuchs/output/igv/F_exon_chain_12.bed")
    print("load /home/tjakobi/work/projects/circRNA/encode_paper/fuchs/output/igv/H_exon_chain_12.bed")
    print("load /home/tjakobi/work/projects/circRNA/encode_paper/fuchs/output/igv/J_exon_chain_12.bed")
    print("load /home/tjakobi/work/projects/circRNA/encode_paper/fuchs/output/igv/L_exon_chain_12.bed")


def print_header(gene):
    print("new")
    print("genome hg38")
    print("snapshotDirectory /tmp/")
    print("maxPanelHeight 5000")
    print("goto "+gene)
    print_tracks()

def print_footer():
    print("collapse")
    print("squish K562.bam")
    print("squish HepG2.bam")
    print("#####################")


def zoom_location(location):
    tmp = location.split(':')
    tmp2 = tmp[1].split('-')
    chromosome = tmp[0]
    start = str(int(tmp2[0]) - 5000)
    stop = str(int(tmp2[1]) + 5000)
    return chromosome+":"+start+"-"+stop

def overview_location(location):
    tmp = location.split(':')
    tmp2 = tmp[1].split('-')
    chromosome = tmp[0]
    start = str(int(tmp2[0]) - 100000)
    stop = str(int(tmp2[1]) + 100000)
    return chromosome+":"+start+"-"+stop

def print_script(data):
    gene_num = 0
    for gene in data:
        num = 0
        for location in data[gene]['loc'].keys():

            print_header(gene)

            for cell in data[gene]['loc'][location].keys():
                for rbp in data[gene]['loc'][location][cell].keys():

                    print("load /home/tjakobi/work/data/circtools/encode_hg38_clip_peaks/" + cell + "/combined/" + rbp + "_" + cell + "_combined.bed")
                    num += 1
                    if num > 4:
                        break
            print_footer()

            print("region " + location)
            print("goto "+overview_location(location))
            print("snapshot "+str(data[gene]['rank'])+"_"+gene+"_"+location+"_gene.png")
            print("goto "+zoom_location(location))
            print("snapshot "+str(data[gene]['rank'])+"_"+gene+"_"+location+"_zoom.png")


        gene_num += 1
        if gene_num > 2:
            break



def parse_file(filename_high):
    #entries_list = {}
    from collections import OrderedDict
    entries_list = OrderedDict()

    with open(filename_high) as fp:
        rank = 1
        for line in fp:
            tmp = line.split('\t')
            name_split = tmp[0].split('_')
            gene = tmp[1]
            location = "chr"+tmp[2] + ":" + tmp[3]  + "-" + tmp[4]
            rbp = name_split[0]
            cell = name_split[1]
            # create key
            if gene not in entries_list:
                entries_list[gene] = {}
                entries_list[gene]['rank'] = rank
                entries_list[gene]['rbp_high'] = {}

                entries_list[gene]['loc'] = OrderedDict()
                entries_list[gene]['loc'][location] = OrderedDict()

                entries_list[gene]['loc'][location]['K562'] = OrderedDict()
                entries_list[gene]['loc'][location]['HepG2'] = OrderedDict()

                entries_list[gene]['loc'][location][cell][rbp] = 1

                rank += 1
            else:
                if location not in entries_list[gene]['loc']:
                    entries_list[gene]['loc'][location] = OrderedDict()
                    entries_list[gene]['loc'][location]['K562'] = OrderedDict()
                    entries_list[gene]['loc'][location]['HepG2'] = OrderedDict()
                else:
                    entries_list[gene]['loc'][location][cell][rbp] = 1


    return entries_list



entries1 = parse_file(args.file_list1)


import pprint
pp = pprint.PrettyPrinter(indent=4)
#pp.pprint(entries1)

print_script(entries1)
