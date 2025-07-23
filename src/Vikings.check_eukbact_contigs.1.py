
#! /usr/bin/env python3.6

'''
Based on sequence similarity search and taxonomic binning, I detect Horizontal Gene Transfer in genes,
by determining if a putatively bacterial gene is surrounded by eukaryotic genes.

Andrea Del Cortona
2019/06/14
'''


#==================================================================#
#   LOAD LIBRARIES                                                 #
#==================================================================#

import argparse
import os
import sys
import re
import statistics
from datetime import datetime


#==================================================================#
#   INPUT PARSER                                                   #
#==================================================================#

parser = argparse.ArgumentParser(description='''I count the number of good snps in your vcf samples

usage:

python3.6 check_eukbact_contigs.py --input KEVEIKS_names.lst''')

parser.add_argument("--input",
	metavar = 'INPUT',
	action = 'store',
	type = str,
	dest = 'INPUT',
	help = "a list of KVEIKS names",
	required = "TRUE")

args = parser.parse_args()


#==================================================================#
#   FUNCTIONS                                                      #
#==================================================================#

def main():

	os.chdir("/media/DISK2-3TB/03_KVEIK/07_HGT_test")

	# get samples names
	SAMPLE_LIST = []
	with open(args.INPUT) as infile:
		for line in infile:
			line = line.rstrip('\n')
			SAMPLE_LIST.append(line)

	# open corresponding files and run the test for each sample
	for sample in SAMPLE_LIST:
		BACT_F = "{}.SPAdes.redundans.all.maker.proteins.fasta.diamond.tax.bact.lst".format(sample)
		TAX_F = "{}.SPAdes.redundans.all.maker.proteins.fasta.diamond.tax.all".format(sample)
		ANNOT_F = "{}.SPAdes.redundans.all.CDSonly.gff".format(sample)

		# identify contigs of interest
		BACT_LIST = []
		with open(BACT_F, "r") as bactfile:
			for line in bactfile:
				line = line.rstrip('\n')
				BACT_LIST.append(line.split("|", 1)[0])
						
		# get contig of interest
		GENE_LIST = []
		CONTIG_LIST = []
		GENOME_CONT = {}
		with open(ANNOT_F, "r") as annotfile:
			for line in annotfile:
				line = line.rstrip('\n')
				POS1, POS2, POS3, POS4, POS5, POS6, POS7, POS8, POS9 = line.split('\t')
				ADD_GENE = re.sub(r"ID=", "", POS9)
				ADD_GENE = re.sub(r"\:cds\;.*", "", ADD_GENE)
				if POS1 in GENOME_CONT:
					GENOME_CONT[POS1].append(ADD_GENE.split("|", 1)[0])
				else:
					GENOME_CONT[POS1] = [ADD_GENE.split("|", 1)[0]]
				
		for BACT in BACT_LIST:
				for k in GENOME_CONT:
					for y in GENOME_CONT[k]:
						z = re.search(BACT, y)
						if z:
							CONTIG_LIST.append(k)
							GENE_LIST.append([k, y])				

		WHOLE_TAX = []
		with open(TAX_F, "r") as taxfile:
			for line in taxfile:
				line = line.rstrip('\n')
				GENE, *TAX = line.split('\t')				
				WHOLE_TAX.append([GENE.split("|", 1)[0], TAX])		
		
		# prepare the output 
		OUTPUT = "{}.SPAdes.redundans.all.maker.proteins.fasta.diamond.tax.HGT.table".format(sample) 
		with open(OUTPUT, "w+") as outfile:
			for CONTIG in CONTIG_LIST:
				for GENE in GENOME_CONT[CONTIG]:
					if GENE in BACT_LIST:
						outfile.write(CONTIG + '\t' + GENE + '\t' + "match" + '\n')
						continue
					else:
						for INFO in WHOLE_TAX:
							if GENE == INFO[0]:
								outfile.write(CONTIG + '\t' + GENE + '\t' + str(INFO[1]) + '\n')
								continue



#==================================================================#
#   RUN                                                            #
#==================================================================#

# run script and give the running time
if __name__ == '__main__':
	t0 = datetime.now()
	main()
	dt = datetime.now() - t0
sys.stderr.write( "# Time elapsed: %s\n" % dt )
