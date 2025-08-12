#!/usr/bin/env python


'''
I remove null values for coverage.bed files, in order to obtain nice plots in R.

Andrea Del Cortona
2018/12/13
'''



####################################################################
### LOAD LIBRARIES #################################################
####################################################################

import argparse
import functools
import numpy
import os
import re
import shutil
import string
import sys


####################################################################
### INPUT PARSER ###################################################
####################################################################

parser = argparse.ArgumentParser(description='''I remove 0 values

	usage:
	python bed_remove0coverage.py --input CD-HIT-EST.fa > CD-HIT-EST.clstr''')

parser.add_argument("--input",
	metavar='INPUT',
	action = 'store',
	type = str,
	dest = 'INPUT',
	help = 'the input bed file',
	required = True)

args = parser.parse_args()


####################################################################
### RUN ############################################################
####################################################################

with open(args.INPUT, 'r') as infile:
	
	for line in infile:

		uno, due, tre, quattro, cinque = line.split()
		
		if cinque != "0.0000000" :

			print(line)	
