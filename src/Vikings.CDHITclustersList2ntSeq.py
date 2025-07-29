#!/usr/bin/env python3.6

'''
I take a fasta file with nucleotide sequences and a folder with list of genes headers.
I return a fasta file for each of the file in the folder with the corresponding sequences.

Andrea Del Cortona
2021/01/27
'''



#------------------------------------------------------------------#
# LOAD LIBRARIES

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



#------------------------------------------------------------------#
# INPUT PARSER

parser = argparse.ArgumentParser(formatter_class = argparse.RawTextHelpFormatter,
	description = '''
	
	==================================================================
	I take a fasta file with nucleotide sequences and a folder with list of genes headers.
	I return a fasta file for each of the file in the folder with the corresponding sequences.
	
	__________________________________________________________________
	
	## all_plusKV.nt.fa example
	
	>14R30_14R3000006713-RA
	ATGACTATTTCTCATCATTTGCGTCATCTTCTAACACCGTATATGATAATATACTAG
	>14R30_14R3000001296-RA
	ATGGACAACTTACAGGTATCTGATATAGAAACTGCTTTACAATGCATATCGTCTACTGCA
	TCTCAAGATGATAAAAACAAAGCGCTTCAATTTTTAGAACAATTCCAAAGATCAACTGTT
	GCCTGGTCTATTTGCAATGAAATATTGTCTAAAGAAGATCCTACAAACGCTCTTCTAGAA



	__________________________________________________________________
	
	Usage:
	python3.6 Vikings.CDHITclustersList2ntSeq.py --fasta all_plusKV.nt.fa --indir ./01_filter_clusters/02_clusters_311-350_noDuplicates/ --outdir ./03_nt_clusters/
	
	==================================================================
	''',
	epilog = '''
	__________________________________________________________________
	         Andrea Del Cortona - andrea.delcortona@gmail.com
			               2021-01-27
	__________________________________________________________________
	''')

parser.add_argument("--fasta", metavar ='FASTA', action = 'store',
	type = str, dest = 'FASTA', required = True,
	help = 'The all_plusKV.nt.fa sequence file with all the nucleotidic sequences.')

parser.add_argument("--indir", metavar ='INDIR', action = 'store',
	type = str, dest = 'INDIR', required = True,
	help = 'The directory with all the deduplicated list of genes for each cluset.')

parser.add_argument("--outdir", metavar ='OUTDIR', action = 'store',
	type = str, dest = 'OUTDIR', required = True,
	help = 'The output directory where to write all the output fasta files.')
	
args = parser.parse_args()



#------------------------------------------------------------------#
# FUNCTIONS                                                      

# main function
def main():
	'''
	I read the fasta file.
	I iterate through the input directory.
	I restrieve fasta sequences for each gene list file in the input directory.
	I write selected fasta sequences to output directory.
	'''

	# read input datasets
	fastaDB = read_fasta(args.FASTA)
	
	# read all list files in input dir
	process_list(args.INDIR, fastaDB)



def read_fasta(fasta):
	'''
	I read an input fasta file.
	I return a dictionary of sequences.
	'''
	
	# declare output dictionary
	fastaDB = {}

	# import clstr file
	with open(fasta) as infile:
		lastSeq = ''
		for line in infile:
			line = line.rstrip('\n')
			# is it an header?
			if line[0] == ">":
				header = line[1:]
				lastSeq = header
				fastaDB[header] = ''
			# is sequence
			else:
				fastaDB[lastSeq] = fastaDB[lastSeq] + line

	return(fastaDB)


	
def process_list(indir, fastaDB):
	'''
	I iterate through the input directory.
	I restrieve fasta sequences for each gene list file in the input directory.
	I write selected fasta sequences to output directory..
	'''
	
	# iterate through the input directory
	for filename in os.listdir(indir):
	
		# read list of gene headers
		with open('%s/%s' % (indir, filename)) as infile:
			GeneList = []
			for line in infile:
				line = line.rstrip('\n')
				GeneList.append(line)
				
		# retrieve fasta seq
		seqDB = {}
		for Gene in GeneList:
			if Gene in fastaDB:
				seqDB[Gene] = ''.join(fastaDB[Gene])
			else:
				sys.stderr.write( "\n# WARNING: the sequence %s is missing from %s \n" % (Gene, filename))
		
		# write fasta to file	
		outfile = filename + '.fa'
		with open('%s/%s' % (args.OUTDIR, outfile), 'w') as oF:
			for seq in seqDB:
				oF.write('>' + seq + '\n')
				oF.write(seqDB[seq] + '\n')
		oF.close()
	


#------------------------------------------------------------------#
# RUN                                                            

# run script and give the running time
if __name__ == '__main__':
	t0 = datetime.now()
	main()
	dt = datetime.now() - t0
sys.stderr.write( "# Time elapsed: %s\n" % dt )
