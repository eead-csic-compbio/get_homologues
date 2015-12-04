#!/usr/bin/python

import sys
import os
import re

Usage="""
./compile_pangenome_output.py pangenome_matrix_t0__*.list.txt

This script will pull all the genes from a list produced by get_homologues script 
"parse_pangenome_matrix.pl" and amalgamate the .faa data into a single .faa. This can be (not by this script) converted into a .gbk file and used in Pathway-tools to identify pathways encoded by each inter-genome group.

Note: Be sure to remove the "compiled_output.txt" before re-running the script, since the script uses the append function.

"""

start = timeit.default_timer()

if len(sys.argv)<1:
	print Usage
else:
	List = sys.argv[1]

#Open file containing all of the file names and kmer values
InFile = open(List, 'r')
Gene_List = []

for Line in InFile:
        Line = Line.strip('\n\r')
        Gene_List.append(Line)

r = re.compile(r'([\d]+)_')

List = re.sub(r'pangenome_matrix_t0__', '', List)

for Line in Gene_List:
	match = r.match(Line)
	os.system(' '.join([
		"cat",
		match.group() + "*",
		">>",
		"output_"+List
	]))

print "All of the sequence data for the genes contained in the list file you provided have been concatenated in the file: output_" +List
stop = timeit.default_timer()

print stop - start
