#!/usr/bin/env python3.5

'''
___________________________________________________

This scripts takes a file with all the Vikings.CNVsmerged.all.tab (resulting from CNVnator_merger.py) concatenated
and evaluates the CNVs that overlaps among the samples.
___________________________________________________

Andrea Del Cortona
2020/06/01
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
import statistics
import sys
import itertools
import operator
from datetime import datetime



#==================================================================#
#   INPUT PARSER                                                   #
#==================================================================#

parser = argparse.ArgumentParser(description='''
This scripts takes a file with all the Vikings.CNVsmerged.all.tab (resulting from CNVnator_merger.py) concatenated
and evaluates the CNVs that overlaps among the samples.

Vikings.CNVsmerged.all.tab looks like:
42R31	I	1001	12000	duplication	3.0931949999999997
42R31	I	13501	24500	deletion	0.08476824999999999
42R31	I	24501	25000	duplication	0.749117
42R31	I	25001	45000	duplication	1.3506749999999998
42R31	I	46001	160000	duplication	1.35033
42R31	I	166001	181500	duplication	1.348675
42R31	I	208001	219000	deletion	0.12166195
42R31	XV	1	12000	deletion	0.0649642



	usage:
	python3.5 Vikings.overlapCNVs.py --allCNVs Vikings.CNVsmerged.all.tab --bed Vikings.CNVsmerged.all.bed > Vikings.CNVs_overlap.txt''')

parser.add_argument("--allCNVs",
	metavar ='ALLCNVS',
	action = 'store',
	type = str,
	dest = 'ALLCNVS',
	help = 'Vikings.CNVsmerged.all.tab.',
	required = True)

parser.add_argument("--bed",
	metavar ='BED',
	action = 'store',
	type = str,
	dest = 'BED',
	help = 'Vikings.CNVsmerged.all.bed.',
	required = True)

args = parser.parse_args()



#==================================================================#
#   FUNCTIONS                                                      #
#==================================================================#

# read the input table and print the output
def main():

	# import CNV file
	with open(args.ALLCNVS) as infile:

		SAMPLES_DB = {}
		SAMPLES_LIST = []
		for line in infile:
			line = line.rstrip('\n')
			SAMPLE, CHR, START, STOP, KIND, FOLD = line.split('\t')
			if SAMPLE in SAMPLES_DB:
				SAMPLES_DB[SAMPLE].append([CHR, START, STOP, KIND, FOLD])
			else:
				SAMPLES_DB[SAMPLE] = [[CHR, START, STOP, KIND, FOLD]]
				SAMPLES_LIST.append(SAMPLE)

	# import windows
	with open(args.BED) as infile:

		WIND_DB = {}
		for line in infile:
			line = line.rstrip('\n')
			CHR, START, STOP = line.split('\t')
			if CHR in WIND_DB:
				WIND_DB[CHR].append([START, STOP])
			else:
				WIND_DB[CHR] = [[START, STOP]]

	CHR_LIST = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "Mito"]
	
	# iterate through chromosomes and overlap windows
	OUTPUT = []
	j = 0
	for CHR in CHR_LIST:
		for POS in WIND_DB[CHR]:
			OUTPUT.append([CHR, POS[0], POS[1]])
			k = 3
			for SAMPLE in SAMPLES_LIST:
				OUTPUT[j].append([])
				#OUTPUT.append([])
				for ENTRY in SAMPLES_DB:
					# if sample is correct
					if ENTRY == SAMPLE:
						for LINE in SAMPLES_DB[ENTRY]:
							# correct chromosome?
							if LINE[0] == CHR:
								START = LINE[1]
								STOP = LINE[2]
								if int(START) >= int(POS[0]) and int(STOP) <= int(POS[1]):
									OUTPUT[j][k].append(str(LINE[4]))
				if OUTPUT[j][k] == []:
					OUTPUT[j][k] = ''
				else:
					OUTPUT[j][k] = ';'.join(OUTPUT[j][k])
				k += 1
			j += 1

	# collapse to the average
	for i in range(len(OUTPUT)):
		for k in range(3, len(OUTPUT[i])):
			if OUTPUT[i][k] != '':
				NUMS = OUTPUT[i][k].split(";")
				NUMS = [float(x) for x in NUMS]
				AVG = statistics.mean(NUMS)
				OUTPUT[i][k] = str(round(AVG, 3))

	# print output
	print('\t'.join(['', '', '', '\t'.join(SAMPLES_LIST)]))
	for k in OUTPUT:
		print('\t'.join(k))


	
#==================================================================#
#   RUN                                                            #
#==================================================================#

# run script and give the running time
if __name__ == '__main__':
	t0 = datetime.now()
	main()
	dt = datetime.now() - t0
sys.stderr.write( "# Time elapsed: %s\n" % dt )

