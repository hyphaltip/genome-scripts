#!/bin/bash
module load stajichlab
module load cdbfasta
cdbfasta $1
cdbyank $1.cidx | sort > $1.names
cdbyank $1.cidx -o $1.reorder < $1.names
