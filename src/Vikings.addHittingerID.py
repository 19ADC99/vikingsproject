#!/usr/bin/env python3.6

'''
I take a list of genes generated with Vikings.filterCDHITclusters.nt.py.
I identify S288C gene name.
I look into a folder of fasta files with Hittigner genes.
I identify Skud & Spar corresponding headers and add them to the gene list.

Andrea Del Cortona
2021/03/13
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

parser.add_argument("--Hit", metavar ='Hit', action = 'store',
	type = str, dest = 'Hit', required = True,
	help = 'The folder with Hittinger gene IDs.')

parser.add_argument("--outdir", metavar ='OUTDIR', action = 'store',
	type = str, dest = 'OUTDIR', required = True,
	help = 'The output directory where to write all the output files.')
	
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

	# loop input files
	for file in os.listdir(args.CDlists):
		outDB = []
		with open('%s/%s' % (args.CDlists, file)) as infile:
			for line in infile:
				line = line.strip('\n')
				outDB.append(line)
				# get Scer gene ID
				if line[:5] == "S288C":
					Prefix, Scer_gene = line.split('_') 
				
			# identify corresponding Hittinger Skud & Spar IDs	
			regex = re.compile('OG[0-9].*_' + Scer_gene)
			for file1 in os.listdir(args.Hit):
				if regex.match(file1):
					with open(args.Hit + "/" + file1) as Hitfile:
						for line in Hitfile:
							line = line.strip('\n')
							if line[:4] == "Skud" or line[0:3] == "Spar":
								outDB.append(line)
				
		# print output file
		outfile = args.OUTDIR + "/" + file + ".Hit"
		oF = open(outfile, 'w')
		
		# print output
		for Entry in outDB:
			oF.write(Entry + '\n')
	
		oF.close()



#------------------------------------------------------------------#
# RUN                                                            

# run script and give the running time
if __name__ == '__main__':
	t0 = datetime.now()
	main()
	dt = datetime.now() - t0
sys.stderr.write( "# Time elapsed: %s\n" % dt )
