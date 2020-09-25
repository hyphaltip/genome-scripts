#!/usr/bin/env python3
from Bio import SeqIO
import argparse, sys, os

parser = argparse.ArgumentParser(description='Count Full Length ORFs from a CDS transcript file.')
parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), help="Input FASTA file",
    default=sys.stdin)

parser.add_argument('--format', "-f", default="fasta", help="Input Sequence format")

args = parser.parse_args()

full_length_ORF = 0
total_ORFs = 0

for record in SeqIO.parse(args.infile, args.format):
    firstcodon = record[0:3].translate().seq
    lastcodon  = record[-3:].translate().seq
    total_ORFs += 1
    if firstcodon == "M" and lastcodon == "*":
        full_length_ORF += 1
    #print(record.id,firstcodon,lastcodon)

print("There are %d ORFs and %d (%.2f%%) are full length"%(total_ORFs,full_length_ORF,100*full_length_ORF/total_ORFs))
