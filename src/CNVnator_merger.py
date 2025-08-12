#!/usr/bin/env python3.5

'''
___________________________________________________

I take two CNVnator files called with different bins (e.g.: 500 bp and 1000 bp) and check concordance, significance and overlap and I print a unified file.
___________________________________________________

Andrea Del Cortona
2019/10/28
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

parser = argparse.ArgumentParser(description='''I take two CNVnator files called with different bins (e.g.: 500 bp and 1000 bp) and check concordance, significance and overlap and I print a unified file.

The input files look like:
deletion	x1156_PM_chrIV:2001-17000	15000	0.00848363	1.06248e-11	0	1.22594e-11	0	1
deletion	x1156_PM_chrIV:44001-82000	38000	0.0501951	4.19401e-12	0	4.42701e-12	0	1

	usage:
	python3.5 CNVnator_merger.py --input_1 CNVnator_bin1 --input_2 CNVnator_bin2 --sample sample_name > CNV_merged.tab''')

parser.add_argument("--input_1",
	metavar ='INPUT_1',
	action = 'store',
	type = str,
	dest = 'INPUT_1',
	help = 'CNVnator output with bin1.',
	required = True)
	
parser.add_argument("--input_2",
	metavar ='INPUT_2',
	action = 'store',
	type = str,
	dest = 'INPUT_2',
	help = 'CNVnator output with bin2.',
	required = True)
	
parser.add_argument("--sample",
	metavar ='SAMPLE',
	action = 'store',
	type = str,
	dest = 'SAMPLE',
	help = 'The name of the sample.',
	required = True)

args = parser.parse_args()



#==================================================================#
#   FUNCTIONS                                                      #
#==================================================================#

# read the input table and print the output
def main():

	# load databases
	with open(args.INPUT_1) as infile:
 		INPUT_1_DB = read_table(infile)

	with open(args.INPUT_2) as infile:
 		INPUT_2_DB = read_table(infile)

	# check positions, chromosome by chromosome. Print overlaps
	for KEY in INPUT_1_DB:
		if KEY in INPUT_2_DB:
			POS_1 = []
			POS_2 = []
			for k in range(len(INPUT_1_DB[KEY])):
				POS_1.append(INPUT_1_DB[KEY][k][0:2])									# here I save the positions of each record
			for k in range(len(INPUT_2_DB[KEY])):
				POS_2.append(INPUT_2_DB[KEY][k][0:2])									# here I save the positions of each record
			
			# iterate through positions and find overlapping regions
			for k in range(len(POS_1)):													# window 1
				for j in range(len(POS_2)):												# window 2
					if int(POS_1[k][1]) < int(POS_2[j][0]):								# window 1 < window 2
						continue
					elif int(POS_1[k][0]) > int(POS_2[j][1]):							# window 1 > window 2
						continue
					else:																# the windows overlap
						if int(POS_1[k][0]) <= int(POS_2[j][0]):
							if int(POS_1[k][1]) <= int(POS_2[j][1]):
								virtual_printer(args.SAMPLE, KEY, POS_2[j][0], POS_1[k][1], INPUT_1_DB[KEY][k][2], INPUT_1_DB[KEY][k][3], INPUT_2_DB[KEY][j][3])
							else:
								virtual_printer(args.SAMPLE, KEY, POS_2[j][0], POS_2[j][1], INPUT_1_DB[KEY][k][2], INPUT_1_DB[KEY][k][3], INPUT_2_DB[KEY][j][3])
						else:
							if int(POS_1[k][1]) <= int(POS_2[j][1]):
								virtual_printer(args.SAMPLE, KEY, POS_1[k][0], POS_1[k][1], INPUT_1_DB[KEY][k][2], INPUT_1_DB[KEY][k][3], INPUT_2_DB[KEY][j][3])
							else:
								virtual_printer(args.SAMPLE, KEY, POS_1[k][0], POS_2[j][1], INPUT_1_DB[KEY][k][2], INPUT_1_DB[KEY][k][3], INPUT_2_DB[KEY][j][3])



# import the input files
def read_table(infile):			

	IN_DB = {}

	for line in infile:
		line = line.rstrip('\n')
		KIND, LOCUS, LENGTH, NORM_RD, EVAL1, EVAL2, EVAL3, EVAL4, QUAL = line.split('\t')
		CHR, COHORD = LOCUS.split(':')
		START, STOP = COHORD.split('-')
		
		if float(EVAL1) < 0.05:
			if CHR in IN_DB:
				IN_DB[CHR].append([START, STOP, KIND, NORM_RD, EVAL1])
			else:
				IN_DB[CHR] = [[START, STOP, KIND, NORM_RD, EVAL1]]
	
	return(IN_DB)
	
	

# print the output
def virtual_printer(SAMPLE, CHR, START, STOP, KIND, NORM_RD1, NORM_RD2):
	MEAN_RD = float((float(NORM_RD1) + float(NORM_RD2))/2)
	print('\t'.join([SAMPLE, CHR, START, STOP, KIND, str(MEAN_RD)]))
	return()
	
	
	
#==================================================================#
#   RUN                                                            #
#==================================================================#

# run script and give the running time
if __name__ == '__main__':
	t0 = datetime.now()
	main()
	dt = datetime.now() - t0
sys.stderr.write( "# Time elapsed: %s\n" % dt )
