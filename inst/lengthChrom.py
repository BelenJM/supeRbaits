#!/usr/bin/env python

# Input file: genome file (FASTA)
# What this script does:
# -Opens the genome file 
# -Count the length of each sequence
# -Writes output
# Output: 
# -column 1, name of contig/chromosome,
# -column 2, length 


import sys
import pandas as pd
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqIO.FastaIO import SimpleFastaParser



# genome file
input_file1 = sys.argv[1]

# define output to return


# main function
def length_genome(file1): 
	with open(input_file1) as genomeFile, open ("temp_folder_for_supeRbaits/genome_size.txt", "w") as output:
		# initializing dataframe
		df = pd.DataFrame(columns=['ChromName','Length'])
		# initializing index
		i = 0
		for values in SimpleFastaParser(genomeFile):
			
			# Names of chromosomes can contain different info
			# we are interested in the first element of the sequence
			# which is assumed to contain the name of the contig/chrom
			# get name of the chromosome
			name_chr = values[0].split(' ') # we first split the name by spaces
			#print(name_chr)
			name_chr = name_chr[0]
			# get the sequence, stored at value[1] of fasta file in generator
			seq = values[1][:]
			#print seq
			#print(len(seq))
			# print to the output file
			output.write("%s\t%s\n" % (name_chr, len(seq)))

			# store into the table to return into R 
			df.loc[i] = (name_chr, len(seq)+1)

			# Index 
			i += 1

		return(df)

print(length_genome(input_file1))