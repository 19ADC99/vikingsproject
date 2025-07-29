#!/usr/bin/env python3.6

'''
I take the output of whatshap trio phasing, and I print a summary of the phasing.

Andrea Del Cortona
2025/03/17
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
    I take the output of whatshap trio phasing, and I print a summary of the phasing.

    __________________________________________________________________
    
    Usage:
    python3.6 Vikings.culture_phasing_summary.py --input phased.vcf --name run_ID
    
    ==================================================================
    ''',
    epilog = '''
    __________________________________________________________________
             Andrea Del Cortona - andrea.delcortona@gmail.com
                           2025-03-17
    __________________________________________________________________
    ''')

parser.add_argument("--input", metavar ='INPUT', action = 'store',
    type = str, dest = 'INPUT', required = True,
    help = 'Input phased trio.')

parser.add_argument("--name", metavar ='NAME', action = 'store',
    type = str, dest = 'NAME', required = True,
    help = 'Name of the samples analysed.')
    
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

    # import vcf file
    with open(args.INPUT) as infile:
        
        # counters
        current_chr = "I"
        CHILD_counter = 0
        N1_parents_counter = 0
        All_parents_counter = 0
        Same_GT_counter = 0

        output_DB = {
            "I": [],
            "II": [],
            "III": [],
            "IV": [],
            "V": [],
            "VI": [],
            "VII": [],
            "VIII": [],
            "IX": [],
            "X": [],
            "XI": [],
            "XII": [],
            "XIII": [],
            "XIV": [],
            "XV": [],
            "XVI": []
        }

        for line in infile:
            
            line = line.rstrip('\n')
            
            # skip header
            if line[0] != "#":
                
                CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, CHILD, MOTHER, FATHER = line.split('\t')
                GT_CHILD, GT_MOTHER, GT_FATHER = CHILD.split(':')[0], MOTHER.split(':')[0], FATHER.split(':')[0]
                
                # check if the trio is phased
                if "|" in GT_CHILD:

                    # get chromosome
                    if current_chr != CHROM:
                        # add to output
                        output_DB[current_chr] = [str(CHILD_counter), str(N1_parents_counter), str(All_parents_counter), str(Same_GT_counter)]
                        # reset counters
                        current_chr = CHROM
                        CHILD_counter = 0
                        N1_parents_counter = 0
                        All_parents_counter = 0
                        Same_GT_counter = 0

                    # get counts
                    CHILD_counter += 1
                    if GT_CHILD == GT_MOTHER and GT_CHILD == GT_FATHER:
                        Same_GT_counter += 1
                        All_parents_counter += 1
                    elif "|" in GT_MOTHER and "|" in GT_FATHER:
                        All_parents_counter += 1
                    elif "|" in GT_MOTHER and "|" not in GT_FATHER or "|" not in GT_MOTHER and "|" in GT_FATHER:
                        N1_parents_counter += 1

            # add last chromosome
            output_DB[current_chr] = [str(CHILD_counter), str(N1_parents_counter), str(All_parents_counter), str(Same_GT_counter)]
                    
    # print output
    for chrom in output_DB:
        print("\t".join([args.NAME, chrom, "\t".join(output_DB[chrom])]))


    
#------------------------------------------------------------------#
# RUN                                                            

# run script and give the running time
if __name__ == '__main__':
    t0 = datetime.now()
    main()
    dt = datetime.now() - t0
sys.stderr.write( "# Time elapsed: %s\n" % dt )
