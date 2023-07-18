#!/bin/usr/python3

# This script is modified Python3 version of Chistoph's PERL script get_circJunction_from_gff_buffered.pl
# the script takes as an input bed and genome fasta file and writes a circular RNA contigs sequences
# Author - Shubhada
# Date - 06-07-2023

########
# INPUT -> named arguments -bed the bed coordinates of BSJs, -db is the genome fasta file, -fo is the output fasta file where the sequences are stored
# OUTPUT -> the output fasta contigs

########
###
# READ ARGUMENTS
import argparse, sys, os
import time
starttime = time.time()

from itertools import groupby
parser = argparse.ArgumentParser()
parser.add_argument("-db", help="The genome fasta file")
parser.add_argument("-bed", help="BED file with BSJs coordinates")
parser.add_argument("-fo", help="BED file with BSJs coordinates")
args = parser.parse_args()
#print(args)
#print(args.db, args.bed)
# check if files provided exist
if not os.path.isfile(args.bed):
    raise ValueError("BED file not found: ", args.bed)
if not os.path.isfile(args.db):
    raise ValueError("DB FASTA file not found: ", args.db)

###
# READ THE FASTA GENOME AND STORE IT IN A HASH
fasta = {}
print("----------- Loading fasta file ----------")
f = open(args.db)
seqs = (x[1] for x in groupby(f, lambda line: line[0] == ">"))
for header in seqs:
    header = header.__next__()[1:].strip().split(" ")[0]  # skip ">"
    seq = "".join(s.strip() for s in seqs.__next__())
    print(header, len(seq))
    fasta[header] = seq
#print(fasta)
print("----------- Loading fasta done ----------")

###
# READ THE BED FILE
bed = open(args.bed).readlines()
# open the file to write the output to
fout = open(args.fo, "w")
for line in bed:
    print(line)
    line = line.strip().split("\t")
    fout.write(">"+"_".join(line)+"\n")
    strand = line[5]
    line[1] = int(line[1])
    line[2] = int(line[2])
    if (strand == "+"):
        seq1 = fasta[line[0]][line[1] : line[1]+100] 
        seq2 = fasta[line[0]][line[1]+100 : line[2]+1] 
        seq = seq2 + seq1 + "\n"
    elif (strand == "-"):
        seq1 = fasta[line[0]][line[2]-100 : line[2]]
        seq2 = fasta[line[0]][line[1] : line[2]-100+1]
        seq = seq1 + seq2 + "\n"
    fout.write(seq)

fout.close()

endtime = time.time()

print("Finished! Time reuqired: ", endtime-starttime, " seconds.")
