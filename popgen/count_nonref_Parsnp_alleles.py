#!/usr/bin/python3

""" 
Count ref or non-ref SNPs from parsnp file (alleles coded as 0,1,2)
"""
import sys
parsnpfile = "parsnp.vcf"
if len(sys.argv) > 1:
	parsnpfile = sys.argv[1]
AllelesColStart = 9
ColCount = 9
AlleleCountsRef = {}
AlleleCountsNonRef = {}
Header = {}

with open(parsnpfile,"r") as fh:
   for line in fh:
        if line.startswith('##'):
            continue
        elif line.startswith('#CHROM'):
            line = line.rstrip()
            Header = line.split("\t")
            ColCount = len(Header)
            for col in range(AllelesColStart,ColCount):
                Header[col] = Header[col].replace(".fasta","")
                AlleleCountsRef[Header[col]] = 0
                AlleleCountsNonRef[Header[col]] = 0

        else:
            row = line.split("\t")
            for n in range(AllelesColStart,ColCount):
                row[n] = row[n].rstrip()
                #print("n is ",n," ",Header[n], " ",row[n])
                if row[n] == "0":
                    AlleleCountsRef[Header[n]] += 1
                else:
                    AlleleCountsNonRef[Header[n]] += 1

print("\t".join(['Strain','RefAllele','NonRefAllele']))
for n in range(AllelesColStart,ColCount):
    strain=Header[n]
    rowcounts = [ strain,repr(AlleleCountsRef[strain]),repr(AlleleCountsNonRef[strain]) ]
    print( "\t".join(rowcounts) )
