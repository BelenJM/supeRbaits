#!/usr/bin/env python

# Opens the txt file, and copies the file into a common txt file
# Opens the fasta file and counts the number of lines


import sys
#import Bio
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqIO.FastaIO import SimpleFastaParser


input_file1 = sys.argv[1]
input_file2 = sys.argv[2]
#output_file1 = sys.argv[2]
# output_file2 = sys.argv[4]


# Assuming that the name of the input file is coord_chrom1,
# being the second element in the name the chromosome name
# we extract the chromosome name
#chromosome_name = str(sys.argv[2].split('_')[1])
#print(chromosome_name)
#chromosome_name2 = str(chromosome_name.split('.')[0])

#print(chromosome_name2)

def processing_baitfile(file1, file2): 
	with open(input_file1) as genomeFile, open ("output.txt", "w") as output:
		with open(input_file2)	as bait_coord:
		
			# we find the chromosome of the input file
			#chromosome_name = str(input_file2.split('_')[1])
			#print(chromosome_name)
			#chromosome_name2 = str(chromosome_name.split('.')[0])
			chromosome_name2 = "CM003279"

			for values in SimpleFastaParser(genomeFile):
				print(values)
				print(len(values[0]))
			# Names of chromosomes can contain different info
			# we are interested in the first element of the sequence
			# which is assumed to contain a name
				name_chr = values[0].split(' ') # we first split the name by spaces
				print(name_chr)
				if len(name_chr) > 1:
					print("length is more than 1")
					name_chr = name_chr[0]
				else:
					name_chr = name_chr 
				
				print(name_chr)
				
				if name_chr == chromosome_name2:
					seq = values[1] [0:-1]
					print(seq)
					StopIteration
			print("now here")
		
		# we extract the sequences from the chromosome and measure the stats
			for i, line in enumerate(bait_coord):
				# skip header
				if i == 0: # change if we do not want header
					pass
				else:
					chromosome, start_bp, end_bp = line.strip('\n').split('\t')
					start_bp = int(start_bp)
					end_bp = int(end_bp)
					seq_bait = seq[start_bp:end_bp]
					length_bait = end_bp - start_bp +1
					output.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (i, chromosome, start_bp, end_bp, seq_bait, length_bait))
			return(start_bp)


processing_baitfile(input_file1, input_file2)
# def getGC(seq):
#     return(str(int(round((sum([1.0 for nucl in seq if nucl in ['G', 'C']]) / len(seq)) * 100))))

# def processFasta(fas):
#     print("identifier\tGCcontent\tlength")
#     for record in SeqIO.parse(fas, "fasta"):
#         print(getGC(str(record.seq)), str(len(record.seq)))

#processFasta(input_file1)
# seq = ''
# seq_extr = ''
# for line in chrom['CM003279.1']:
# 	line = str(line)
# 	seq += line.rstrip().replace('\n', '')
# 	seq_extr = seq[1:5]
# #print(seq)
# #print(seq)
# print(len(seq))
# print(seq_extr)
# with open("exo.txt", "r") as f, open ("test.txt", "w") as output :
#     seq =""
#     for line in f :
#         if line.startswith(">"):
#             if seq != "":
#                 output.write(seq)
#                 seq=""
#             seq += '\n'+line     
#         else:
#             seq += line.strip().replace('\n', '')
#     else:
#         #end of lines 
#          if seq != "":
#                 output.write(seq)
#                 seq=""