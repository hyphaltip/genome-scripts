#!/usr/bin/env python3
import os,sys,re

# split a hmm db into pieces of individual hmmer files
inputfile = sys.argv[1]

if os.path.isfile(inputfile):
	print("infile",inputfile)

hmms = {}
with open(inputfile,"r") as fh:
	lines = []
	name  = ""
	for line in fh:
		if re.match("^HMMER3",line):
			if len(name) and len(lines):
				hmms[name] = lines
				lines = []
				name = ""
		elif re.match("^NAME",line):
			name = re.split("\s+",line)[1]
			name = name.strip()
		lines.append(line)
	hmms[name] = lines

for name in hmms:
	data = hmms[name]
	if not name.endswith(".hmm"):
		name += ".hmm"

	with open(name, "w") as fh:
		for l in data:		
			fh.write(l)
