#!/usr/bin/env python3.5

'''
___________________________________________________

This scripts takes a vcf files containing a gene of interest in the yeast domestication process
and it outputs a matrix reporting interesting genotypes and coverages for a set of samples.
___________________________________________________

Andrea Del Cortona
2020/06/23
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

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
	description='''
	This scripts takes a vcf files containing a gene of interest in the yeast domestication process	
	and it outputs a matrix reporting interesting genotypes and coverages for a set of samples.

	vcf_sample_position.txt is a tab delimtied file, on the first column is the name of the sample
	with the order they appear on the vcf file, on the second column there is metannotation:
	e.g.:
	_____________________
	
	Laerdal2	Kveik
	Muri	Kveik
	SortdalEbbe1	Kveik
	Voss1	Kveik
	X1002	Beer1
	X1007	Beer2
	X1011	Mixed
	X1014	Beer1
	_____________________

	usage:
	python3.5 Vikings.domesticationVCFparser.py --input gene.vcf --samples vcf_sample_position.txt > gene.matrix.txt''')

parser.add_argument("--input",
	metavar ='INPUT',
	action = 'store',
	type = str,
	dest = 'INPUT',
	help = 'A gene.vcf file.',
	required = True)

parser.add_argument("--samples",
	metavar ='SAMPLES',
	action = 'store',
	type = str,
	dest = 'SAMPLES',
	help = 'List of samples names wiht the order they appear on the vcf file.',
	required = True)

args = parser.parse_args()



#==================================================================#
#   FUNCTIONS                                                      #
#==================================================================#

# read the input table and print the output
def main():

	# import VCF file
	with open(args.INPUT) as infile:

		GENE_DB = []
		for line in infile:
			line = line.rstrip('\n')
			CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *SAMPLES = line.split('\t')
			GENE_DB.append([[CHROM, POS],[REF, ALT, FILTER, INFO, FORMAT], SAMPLES])
			
	# import sample list (vcf_sample_position.txt)
	with open(args.SAMPLES) as infile:
		
		GROUP, GENE, *OTHER = args.INPUT.split('.')
		SAMPLE_LIST = []
		for line in infile:
			line = line.rstrip('\n')
			SAMPLE, ANNOT = line.split('\t')
			SAMPLE_LIST.append([GROUP, GENE, SAMPLE, ANNOT])
	
	# gene walker
	OUTPUT_MATRIX = []
	NO_COV_POS = []
	NO_ALT_POS = []
	NO_ALT_COV = ["REF" for x in range(len(SAMPLE_LIST))]
	ALT_POS = ["REF" for x in range(len(SAMPLE_LIST))]
	for POSIT in GENE_DB:
		
		# check format
		# GT ---> 1511, no coverage in all the samples
		if POSIT[1][4] == "GT":
			# reset no alternative
			if NO_ALT_POS != []:
				if "NO_COV" not in set(NO_ALT_COV):
					OUTPUT_MATRIX.append([NO_ALT_POS, "no alternative start", ["REF" for x in range(len(SAMPLE_LIST))]])
					#OUTPUT_MATRIX.append([POSIT[0], "no alternative stop", ["REF" for x in range(len(SAMPLE_LIST))]])
				else:
					OUTPUT_MATRIX.append([NO_ALT_POS, "no alternative start +", NO_ALT_COV])
					#OUTPUT_MATRIX.append([POSIT[0], "no alternative stop +", NO_ALT_COV])
				NO_ALT_POS = []
				NO_ALT_COV = ["REF" for x in range(len(SAMPLE_LIST))]
			# process position
			if NO_COV_POS == []:
				NO_COV_POS = POSIT[0]
				
		
		# GT:DP ---> 204302, no variants any sample OR no variants + no coverage
		elif POSIT[1][4] == "GT:DP":
			# reset no coverage all samples
			if NO_COV_POS != []:
				OUTPUT_MATRIX.append([NO_COV_POS, "no coverage start", ["NO_COV" for x in range(len(SAMPLE_LIST))]])
				#OUTPUT_MATRIX.append([POSIT[0], "no coverage stop", ["NO_COV" for x in range(len(SAMPLE_LIST))]])
				NO_COV_POS = []
			# process position
			if NO_ALT_POS == []:
				NO_ALT_POS = POSIT[0]
				for k in range(len(SAMPLE_LIST)):
					if POSIT[2][k] == "./.:.":
						NO_ALT_COV[k] = "NO_COV"
			else:
				NEW_ALT_COV = ["REF" for x in range(len(SAMPLE_LIST))]
				for k in range(len(GENE_DB[0][2])):
					if POSIT[2][k] == "./.:.":
						NEW_ALT_COV[k] = "NO_COV"
						if NEW_ALT_COV != NO_ALT_COV:
							OUTPUT_MATRIX.append([NO_ALT_POS, "no alternative start +", NO_ALT_COV])
							#OUTPUT_MATRIX.append([POSIT[0], "no alternative stop +", NO_ALT_COV])    # # # ## # I MAY WANT TO CHANGE THIS AFTERWARDS
							NO_ALT_POS = POSIT[0]
							NO_ALT_COV = NEW_ALT_COV		


		# GT:AD:DP ---> 1, ignore them, just one positions: [VIII:526893]
		elif POSIT[1][4] == "GT:AD:DP":
			continue
			
			
		# GT:AD:DP:GQ ---> 534, multiallelic variant + no coverage
		elif POSIT[1][4] == "GT:AD:DP:GQ":
			# reset no coverage all samples
			if NO_COV_POS != []:
				OUTPUT_MATRIX.append([NO_COV_POS, "no coverage start", ["NO_COV" for x in range(len(SAMPLE_LIST))]])
				#OUTPUT_MATRIX.append([POSIT[0], "no coverage stop", ["NO_COV" for x in range(len(SAMPLE_LIST))]])
				NO_COV_POS = []	
			# reset no alternative
			if NO_ALT_POS != []:
				if "NO_COV" not in set(NO_ALT_COV):
					OUTPUT_MATRIX.append([NO_ALT_POS, "no alternative start", ["REF" for x in range(len(SAMPLE_LIST))]])
					#OUTPUT_MATRIX.append([POSIT[0], "no alternative stop", ["REF" for x in range(len(SAMPLE_LIST))]])
				else:
					OUTPUT_MATRIX.append([NO_ALT_POS, "no alternative start +", NO_ALT_COV])
					#OUTPUT_MATRIX.append([POSIT[0], "no alternative stop +", NO_ALT_COV])
				NO_ALT_POS = []
				NO_ALT_COV = ["REF" for x in range(len(SAMPLE_LIST))]	
			# process position
			*ALT, = POSIT[1][1].split(',')
			NALT = len(ALT)
			for k in range(len(SAMPLE_LIST)):
				GT, AD, DP, GQ = POSIT[2][k].split(':')
				if POSIT[2][k] == "./.:.:.:.":
					ALT_POS[k] = "NO_COV"
				elif GT != "0/0":
					ALT1, ALT2 = GT.split('/')
					# ALT_POS[k] = "ALT" + str(ALT1) + "," + "ALT" + str(ALT2)
					ALT_POS[k] = GT
			OUTPUT_MATRIX.append([POSIT[0], "multiallelic locus", ALT_POS])
			ALT_POS = ["REF" for x in range(len(SAMPLE_LIST))]
			
			############ ALT_POS = "REF" for x in range(len(GENE_DB[0][2]))]
			
			### check impact of variatn


		# GT:AD:DP:GQ:PL ---> 33537, biallelic variant + no covaerage
		elif POSIT[1][4] == "GT:AD:DP:GQ:PL":
			# reset no coverage all samples
			if NO_COV_POS != []:
				OUTPUT_MATRIX.append([NO_COV_POS, "no coverage start", ["NO_COV" for x in range(len(SAMPLE_LIST))]])
				#OUTPUT_MATRIX.append([POSIT[0], "no coverage stop", ["NO_COV" for x in range(len(SAMPLE_LIST))]])
				NO_COV_POS = []	
			# reset no alternative
			if NO_ALT_POS != []:
				if "NO_COV" not in set(NO_ALT_COV):
					OUTPUT_MATRIX.append([NO_ALT_POS, "no alternative start", ["REF" for x in range(len(SAMPLE_LIST))]])
					#OUTPUT_MATRIX.append([POSIT[0], "no alternative stop", ["REF" for x in range(len(SAMPLE_LIST))]])
				else:
					OUTPUT_MATRIX.append([NO_ALT_POS, "no alternative start +", NO_ALT_COV])
					#OUTPUT_MATRIX.append([POSIT[0], "no alternative stop +", NO_ALT_COV])
				NO_ALT_POS = []
				NO_ALT_COV = ["REF" for x in range(len(SAMPLE_LIST))]	
			# process position
			for k in range(len(SAMPLE_LIST)):
				GT, AD, DP, GQ, PL = POSIT[2][k].split(':')
				if POSIT[2][k] == "./.:.:.:.:.":
					ALT_POS[k] = "NO_COV"
				elif GT != "0/0":
					ALT1, ALT2 = GT.split('/')
					# ALT_POS[k] = "ALT" + str(ALT1) + "," + "ALT" + str(ALT2)
					ALT_POS[k] = GT
			OUTPUT_MATRIX.append([POSIT[0], "biallelic locus", ALT_POS])
			ALT_POS = ["REF" for x in range(len(SAMPLE_LIST))]
			
			############ ALT_POS = "REF" for x in range(len(GENE_DB[0][2]))]
			
			### check impact of variatn

	'''
	#### OUTPUT
	## stupid output without ranges
	for ENTRY in OUTPUT_MATRIX:
		if set(ENTRY[2]) == REF:
			continue
		else:
			for k in range(len(SAMPLE_LIST)):
				#SAMPLE_LIST[k].append('\t'.join(ENTRY[0]))
				SAMPLE_LIST[k].append(ENTRY[2][k])
	
	for k in SAMPLE_LIST:
		print('\t'.join(k))
	#print(OUTPUT_MATRIX)
	'''
	#### OUTPUT
	## stupid output without ranges
	OUTPUT = []
	COUNT = 1
	for ENTRY in OUTPUT_MATRIX:
		if set(ENTRY[2]) == REF:
			continue
		else:
			for k in range(len(SAMPLE_LIST)):
				OUTPUT.append(['\t'.join(SAMPLE_LIST[k]), str(COUNT), ENTRY[2][k]])
			COUNT += 1
	
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
