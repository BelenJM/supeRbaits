#!/usr/bin/env python3

from sys import stdout
from argparse import ArgumentParser

from pandas import DataFrame
from Bio.SeqIO.FastaIO import SimpleFastaParser

header = ["ChromName", "Length"]

def row_builder(index, name, seq):
    """Build a row of the data frame.

    The row consists of the name of the chromosome,
    and the length of its sequence.

    Keyword arguments:
    index -- the index of the chromosome (starting at 0)
    name  -- the name of the chromosome possibly with some noise
    seq   -- the sequence of the chromosome
    """
    get_chrom_name = lambda x: x.split(" ", 1)[0]

    return [get_chrom_name(name), len(seq)]

def get_lengths(genome_file_path, df_header, df_row_builder):
    """Construct a data frame with columns according to df_header and rows according to df_row_builder.

    The data frame is written to output_file_path.

    Keyword arguments:
    genome_file_path -- the path to the genome file
    df_header        -- the data frame header as a list of column names
    df_row_builder   -- a function that builds a row for the data frame given
                        the index, name and sequence of a chromosome in the genome file
    """
    df = DataFrame(columns = df_header)

    with open(genome_file_path, "r") as genome_file_handler:
        for index, (name, seq) in enumerate(SimpleFastaParser(genome_file_handler)):
            df.loc[index] = df_row_builder(index, name, seq)
                        
    return df

def main():
    parser = ArgumentParser(description="Extract chromosome lengths from genome in FASTA format")
    parser.add_argument("fasta_file_path", type=str, help="Path to the FASTA file")
    args = parser.parse_args()

    get_lengths(args.fasta_file_path, header, row_builder).to_csv(stdout, sep="\t", index=None, header=False)

if __name__ == "__main__":
    main()
