#!/usr/bin/env python3.6

"""
I filter Orthofinder Orthogroups tabular output.
I identify kveiks specific Orthogroups (if any).

Andrea Del Cortona
2023/03/20
"""



#------------------------------------------------------------------#
# LOAD LIBRARIES

import argparse
import sys
from datetime import datetime



#------------------------------------------------------------------#
# INPUT PARSER

parser = argparse.ArgumentParser(formatter_class = argparse.RawTextHelpFormatter,
	description = """
	==================================================================
	I filter Orthofinder Orthogroups tabular output.
	I identify kveiks specific Orthogroups (if any).
	__________________________________________________________________
	Usage:
	python3.8 Vikings.GeneFamiliesFilter.py --in Vikings.Orthogroups.GeneCount.tsv
	==================================================================
	""",
	epilog = """
	__________________________________________________________________
			 Andrea Del Cortona - andrea.delcortona@gmail.com
						   2023-03-20
	__________________________________________________________________
	""")

parser.add_argument("--in", metavar = "INFILE", action = "store",
	type = str, dest = "INFILE" , required = True,
	help = "Input Vikings.Orthogroups.GeneCount.tsv file.")

args = parser.parse_args()



#------------------------------------------------------------------#
# FUNCTIONS

# main function
def main():
	"""
	I filter Orthofinder Orthogroups tabular output.
	I identify kveiks specific Orthogroups (if any).
	"""

	# declare list of kveiks
	kveiks_DB = ["41R10", "21R38", "9R40", "17P5", "SortdalEbbe1", "3R11", "21P1", "41R15", "Hornindal1", "Hornindal2",
	    "1R16", "2R23", "8R19", "Muri", "k7R25", "38R16", "44R32", "19R18", "44R7", "6R15", "Laerdal2", "7R7", "14R6",
		"14R30", "27R17", "28P1", "28P6", "28R21", "28R33", "28R8", "42R20", "42R31", "45P5", "45R11", "46R12", "46R37",
		"16R23", "16R37", "39R20", "40R14", "40R1", "40R20", "Granvin1", "Voss1"]
	
	# declare Orthologs database
	OrthoDB = {}
	Kveiks_Orthogroups = {}

	# import all Orthogroups
	with open(args.INFILE) as infile:
		for line in infile:
			line = line.rstrip("\n")

			# get sample list
			if line[0:10] == "Orthogroup":
				sample_list = line.split("\t")[1:-1]
				continue
			
			NAME, *SAMPLES = line.split("\t")[:-1]
			OrthoDB[NAME] = SAMPLES

	# find kveiks unique orthogroups
	for ORTHOGROUP in OrthoDB:
		kveiks_count = 0
		others_count = 0
		for k in range(len(OrthoDB[ORTHOGROUP])):
			if int(OrthoDB[ORTHOGROUP][k]) > 0:
				if sample_list[k] in kveiks_DB:
					kveiks_count += 1
				else:
					others_count += 1

		if kveiks_count >= 4 and others_count <= 40:
			Kveiks_Orthogroups[ORTHOGROUP] = OrthoDB[ORTHOGROUP]

	# print output
	print("\t".join(["Orthogroups", "\t".join(sample_list)]))
	for ORTHOGROUP in Kveiks_Orthogroups:
		print("\t".join([ORTHOGROUP, "\t".join(Kveiks_Orthogroups[ORTHOGROUP])]))



#------------------------------------------------------------------#
# RUN

# run script and give the running time
if __name__ == '__main__':
	t0 = datetime.now()
	main()
	dt = datetime.now() - t0
	sys.stderr.write("# Time elapsed: %s\n" % dt)
