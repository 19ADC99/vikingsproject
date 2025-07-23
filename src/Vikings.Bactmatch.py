#!/usr/bin/env python3.5

'''
___________________________________________________

I scan the output of check_eukbact_contigs.1.py and filter for putative
bacterial genes surrounded by eukaryotic genes.
___________________________________________________

Andrea Del Cortona
2020/06/15
'''



#==================================================================#
#   LOAD LIBRARIES                                                 #
#==================================================================#

import argparse
import functools
import os
import re
import shutil
import string
import sys
import itertools
import operator
from datetime import datetime



#==================================================================#
#   INPUT PARSER                                                   #
#==================================================================#

parser = argparse.ArgumentParser(description='''
I scan the output of check_eukbact_contigs.1.py and filter for putative
bacterial genes surrounded by eukaryotic genes.

The inputfile looks like:
scaffold70|size63194	Voss100003601-RA	['Eukaryota']
scaffold70|size63194	Voss100003602-RA	match
scaffold70|size63194	Voss100003603-RA	[]

	usage:
	python3.5 Vikings.Bactmatch.py --input > bacterialgenes.txt''')

parser.add_argument("--input",
	metavar ='INPUT',
	action = 'store',
	type = str,
	dest = 'INPUT',
	help = 'The output of check_eukbact_contigs.1.py script.',
	required = True)

args = parser.parse_args()



#==================================================================#
#   FUNCTIONS                                                      #
#==================================================================#

# read the input table and print the output
def main():

	# import input file
	with open(args.INPUT) as infile:

		INPUT_DB = []
		for line in infile:
			line = line.rstrip('\n')
			CHR, GENE, FEAT = line.split('\t')
			INPUT_DB.append([CHR, GENE, FEAT])
	
	# discover bacterial positions surrounded by eukaryotic genes
	for k in range(len(INPUT_DB)):
		if INPUT_DB[k][2] == "match":
			if INPUT_DB[k][0] == INPUT_DB[k-1][0] and INPUT_DB[k-1][2] == "['Eukaryota']":
				print('\t'.join(INPUT_DB[k]))	
			elif INPUT_DB[k][0] == INPUT_DB[k+1][0] and INPUT_DB[k+1][2] == "['Eukaryota']":
				print('\t'.join(INPUT_DB[k])) 


	
#==================================================================#
#   RUN                                                            #
#==================================================================#

# run script and give the running time
if __name__ == '__main__':
	t0 = datetime.now()
	main()
	dt = datetime.now() - t0
sys.stderr.write( "# Time elapsed: %s\n" % dt )
