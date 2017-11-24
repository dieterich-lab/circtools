#! /usr/bin/env python3

# Copyright (C) 2017 Franziska Metge, Tobias Jakobi
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


# script to permutate motif files

# define functions
def read_in_motif_file(infile):
    motifs = {}
    input_motif = open(infile)
    header = ''
    for motif in range(0, 9):
        header += input_motif.readline()
    motif = ''
    counter = 0
    for Line in input_motif:
        if 'MOTIF' in Line:
            motif = Line
            if motif in motifs:
                motif = motif + str(counter)
                counter += 1
        elif 'letter-probability' in Line:
            read_motif = True
            motifs[motif] = [Line]
        elif Line == '\n':
            read_motif = False
        elif Line.count('\t') == 4 and read_motif:
            motifs[motif] += [Line]
    input_motif.close()
    return motifs, header


def write_permutated_motifs(outfile, motif_list, header):
    output_file = open(outfile, 'w')
    output_file.write(header)
    for motif in motif_list:
        tmp = motif_list[motif][1:]
        random.shuffle(tmp)
        output_file.write('%s\n\n%s%s\n' % (motif.split('\n')[0], motif_list[motif][0], ''.join(tmp)))
    output_file.close()
    return


# run script
if __name__ == '__main__':

    # required packages
    import random
    import argparse

    # input
    parser = argparse.ArgumentParser(description='Takes a meme formatted motif file and permutates each motif n times.')

    parser.add_argument('meme_file', metavar='motifs.meme',
                        help='Meme formatted files with motifs to permutate for p-value calculations of motif enrichment.')
    parser.add_argument('output_prefix', metavar='prefix_for_path',
                        help='prefix for for all permuation files to. Permuations will be written to PREFIX.permuation.i.meme')
    # options
    parser.add_argument('-n', dest='repeats', default=100, help='number of permutations to perform')

    # parse arguments
    args = parser.parse_args()

    meme_file = args.meme_file
    output_prefix = args.output_prefix
    permuations = args.repeats

    # run
    MOTIFS, HEADER = read_in_motif_file(meme_file)
    for i in range(0, permuations):
        write_permutated_motifs('%s.permutation.%s.meme' % (output_prefix, i), MOTIFS, HEADER)
