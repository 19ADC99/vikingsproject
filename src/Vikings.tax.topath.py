#! /home/andrea/anaconda_ete/bin/	python3.6


'''
I retrieve the full taxonomic path based on the NCBI taxid.

Andrea Del Cortona
2019/02/07
'''


####################################################################
### LOAD LIBRARIES #################################################
####################################################################

import argparse
import functools
#import numpy
import os
#import pandas
import re
import shutil
import string
import sys
import csv
from ete3 import NCBITaxa
ncbi = NCBITaxa()


####################################################################
### INPUT PARSER ###################################################
####################################################################

parser = argparse.ArgumentParser(description='''I assign full taxomonic path to a the DIAMOND output

	usage:
	~/anaconda_ete/bin/python3.6 tax.topath.py --input DIAMOND.output.tab >> OUTPUT.tab''')

parser.add_argument("--input",
	metavar='INPUT',
	action = 'store',
	type = str,
	dest = 'INPUT',
	help = 'the DIAMOND.output.tab input file in format 6, with taxid as last (13th) column',
	required = True)

args = parser.parse_args()


####################################################################
### FUNCTIONS ######################################################
####################################################################

def get_desired_ranks(taxid):
	lineage = ncbi.get_lineage(taxid)
	names = ncbi.get_taxid_translator(lineage)
	return [names[taxid] for taxid in lineage]

def main(taxids):
	mytax = get_desired_ranks(taxids)
	return mytax


####################################################################
### RUN ############################################################
####################################################################

with open(args.INPUT, 'r') as input_file:

	for line in input_file:
		
		line = line.rstrip('\n')
		prot, AC, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, staxids = line.split('\t')
		main_taxid, *throwaway = staxids.split(";") 		

		if main_taxid != '':

			# generates the lineage
			taxids = main_taxid
			mytax = main(taxids)		
		
			# print output
			mytax.insert(0, AC)
			mytax.insert(0, prot)		
			print(*mytax, sep = '\t')	

		else:

			mytax = ["taxid not available"]
			mytax.insert(0, AC)
			mytax.insert(0, prot)		
			print(*mytax, sep = '\t')	




