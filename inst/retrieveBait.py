#!/usr/bin/env python

# Opens the txt file, and copies the file into a common txt file
# Opens the fasta file and counts the number of lines


import sys
import pandas as pd
#import Bio
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqIO.FastaIO import SimpleFastaParser


# Input from the console
input_file1 = sys.argv[1] # genome file
input_file2 = sys.argv[2] # bait BED file


def counting_nt(sequence, a_count, t_count, g_count, c_count, at_count, gc_count, unknown):
	'''
	Function to count the number of nucleotides (A, T, G, C, GC, AT, UNK)
	from a given sequence '''

	nuc_str = list(sequence.strip())
	for n in nuc_str:
		if n == 'G' or n == 'g':
			g_count += 1
		elif n == 'C' or n == 'c':
			c_count += 1
		elif n == 'A' or n == 'a':
			a_count += 1
		elif n == 'T' or n == 't':
			t_count += 1
		else:
			unknown += 1  # usually represented by an 'N' in the fasta file.
	# count the number of AT
	at_count = a_count + t_count
	# count the number of GC
	gc_count = g_count + c_count

	return(a_count, t_count, g_count, c_count, at_count, gc_count, unknown)
    	

def processing_baitfile(file1, file2): 

	''' 
	Function to take the given bait file generated in the R script,
	find the chromosome name in the genome file,
	extract the sequence, count nucleotides,
	generate output file + dataframe in R
	'''

	# initializing the dataframe to return to R
	df = pd.DataFrame(columns=['Bait_no', 'ChromName','Start_bp', 'End_bp', 'Sequence_bait', 'Size_bait', 'Number_A', 'Number_T', 'Number_G', 'Number_C', 'Number_UNK', 'Number_AT', 'Number_GC'])
	
	with open(input_file1) as genomeFile, open ("output.txt", "w") as output:
		# prepare the output file with header
		output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("Bait_no", "Chrom", "Start_bp","End_bp","Sequence_bait", "Length_bait", "No_A", "No_T", "No_G", "No_C", "No_Unk", "No_AT", "No_GC"))

		with open(input_file2)	as bait_coord:
		
			# we find the chromosome of the input file
			# if name is "bed_CHROMOSOMENAME", remove "bed_"
			chromosome_name = str(input_file2.split('_')[1])
			#print(chromosome_name)
			# remove the extension (e.g. ".txt")
			chromosome_name2 = str(chromosome_name.split('.')[0])
			#print chromosome_name2
			# for testing: 
			# chrom_name_bait = "CM003279"

			# loop over the FASTA file using SimpleFastaParser from BIO
			for values in SimpleFastaParser(genomeFile):
			# Names of chromosomes can contain different info
			# We are interested in the first element of the sequence
			# which is assumed to contain a name
				name_chr = values[0].split(' ')[0]# we first split the name by spaces
#				print(name_chr)
#				print("Name of chr is, ", name_chr, " and name of the other chrom is ", chromosome_name2)
				
				if str(name_chr) == str(chromosome_name2):
					#print("hey! it is the same name")
					seq = values[1][:]
					#print(seq)
					break # we stop looking in the chromosome file as soon as the chrom names match
			#print("now here")
		
		# we extract the sequences from the chromosome and measure the stats
			for i, line in enumerate(bait_coord):
				# skip header
				if i == 0: # change if we do not want header
					pass
				else:
					chromosome, start_bp, end_bp = line.strip('\n').split('\t')
					start_bp = int(start_bp)
					end_bp = int(end_bp)
					# extract the sequence from the genome file
					seq_bait = seq[start_bp:end_bp]
					# count the length of the bait sequence
					length_bait = end_bp - start_bp +1
					
					# count the nucleotides
					a_count = 0
					t_count = 0
					g_count = 0
					c_count = 0
					at_count = 0
					gc_count = 0
					unknown = 0
					a, t, g, c, unk, at, gc = counting_nt(seq_bait, a_count, t_count, g_count, c_count, unknown, at_count, gc_count)
					
					# OUTPUT:
					output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (i, chromosome, start_bp, end_bp, seq_bait, length_bait, a, t, g, c, unk, at, gc))
					df.loc[i] = (i, chromosome, start_bp, end_bp, seq_bait, length_bait, a, t, g, c, unk, at, gc)

			return(df)


processing_baitfile(input_file1, input_file2)
