#!/usr/bin/env python3.6

'''
I take a list of genes generated with Vikings.filterCDHITclusters.nt.py.
I identify S288C gene name.
I look into a folder of fasta files with Hittigner genes.
I identify Skud & Spar corresponding headers and add them to the gene list.

Andrea Del Cortona
2024/11/26
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
	I take a list of genes generated with Vikings.filterCDHITclusters.nt.py.
	I identify S288C gene name.
	I look into a folder of fasta files with Hittigner genes.
	I identify Skud & Spar corresponding headers and add them to the gene list.
	
	__________________________________________________________________
	
	Usage:
	python3.6 Vikings.addHittingerID.py --CDlists 01_nt_clusters_lst_select --Hit 99_Hittinger/geneList/ --outdir 02_nt_clusters_lst_select_Hittinger
	
	==================================================================
	''',
	epilog = '''
	__________________________________________________________________
	         Andrea Del Cortona - andrea.delcortona@gmail.com
			               2021-03-13
	__________________________________________________________________
	''')

parser.add_argument("--CDlists", metavar ='CDlists', action = 'store',
	type = str, dest = 'CDlists', required = True,
	help = 'Folder with gene ID list generated with Vikings.filterCDHITclusters.nt.py.')

parser.add_argument("--Genes", metavar ='Genes', action = 'store',
	type = str, dest = 'Genes', required = True,
	help = 'Genes to filter.')
	
args = parser.parse_args()



#------------------------------------------------------------------#
# FUNCTIONS                                                      

# main function
def main():
	'''
	I take a list of genes generated with Vikings.filterCDHITclusters.nt.py.
	I identify S288C gene name.
	I look into a folder of fasta files with Hittigner genes.
	I identify Skud & Spar corresponding headers and add them to the gene list.
	'''


	sequences_DB = {}
	for file in os.listdir(args.CDlists):
		sequences_DB[file] = []
		with open('%s/%s' % (args.CDlists, file)) as infile:
			for line in infile:
				line = line.strip('\n')
				sequences_DB[file].append(line)

	geneList = []
	with open(args.Genes) as infile:
		for line in infile:
			line = line.strip('\n')
			geneList.append(line)

	for Gene in geneList:
		print(Gene)
		for cluster in sequences_DB:
			if "S288C_"+Gene in sequences_DB[cluster]:
				# print output file
				outfile = "./01_nt_clusters_lst_select/" + cluster
				oF = open(outfile, 'w+')

				# print output
				for Sequence in sequences_DB[cluster]:
					oF.write(Sequence + '\n')

				oF.close()



#------------------------------------------------------------------#
# RUN                                                            

# run script and give the running time
if __name__ == '__main__':
	t0 = datetime.now()
	main()
	dt = datetime.now() - t0
sys.stderr.write( "# Time elapsed: %s\n" % dt )
