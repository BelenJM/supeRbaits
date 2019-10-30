#!/usr/bin/env python3

from pathlib import PurePath
from collections import Counter
from argparse import ArgumentParser

from pandas import DataFrame
from Bio.Alphabet import generic_dna
from Bio.SeqIO.FastaIO import SimpleFastaParser

header = [
    "Bait_no",
    "ChromName",
    "Start_bp",
    "End_bp",
    "Sequence_bait",
    "Size_bait",
    "Number_A",
    "Number_T",
    "Number_G",
    "Number_C",
    "Number_UNK",
    "Number_AT",
    "Number_GC"
]

def get_nuc_count(seq):
    """Count the number of A, T, G, C, AT, GC and UNK (unknown nuclueotides) in seq.
    
    Keyword arguments:
    seq -- a dna sequence
    """
    seq = seq.upper()
    count = Counter(seq)

    return {
        "A" : count["A"],
        "T" : count["T"],
        "G" : count["G"],
        "C" : count["C"],
        "AT" : count["A"] + count["T"],
        "GC" : count["G"] + count["C"],
        "UNK" : sum(count.values()) - (count["A"] + count["T"] + count["C"] + count["G"])
    }

def get_chrom_seq(chrom_name, genome_file_path):
    """Find chrom_name in genome_file_path and return its sequence.

    Keyword arguments:
    chrom_name       -- the name of the chromosome
    genome_file_path -- the path to the genome file

    Raises:
    RuntimeError -- if chrom_name does not appear on genome_file_path
    """
    get_chrom_name = lambda x: x.split(" ", 1)[0]

    with open(genome_file_path, "r") as genome_file_handler:
        for name, seq in SimpleFastaParser(genome_file_handler):
            if get_chrom_name(name) == chrom_name:
                return seq
            
    raise RuntimeError(f"{bait_chrom_name} was not found in {genome_file_path}")

def bait_parser(bait_file_handler, header=False):
    """Parse a bait file.

    Keyword arguments:
    bait_file_handler -- file handler of the bait file
    header            -- should yield header (default=False)
    """
    for index, line in enumerate(bait_file_handler):
        if not(header) and index == 0:
            continue
        yield map(int, line.strip().split()[0:2])

def row_builder(index, bait_chrom_name, seq, seq_start, seq_end):
    """Build a row of the data frame.

    The row consists of the index of the chromosome, its name,
    the start of the bait sequence, the end of the bait sequence,
    the bait sequence, the length of the bait sequence,
    the number of A nucleotides,
    the number of T nucleotides,
    the number of G nucleotides,
    the number of C nucleotides,
    the number of UNK (unknown) nucleotides,
    the total number of A and T nucleotides,
    the total number of G and C nucleotides.

    Keyword arguments:
    index           -- the index of the chromosome (starting at 0)
    bait_chrom_name -- the name of the bait chromosome
    seq             -- the sequence of the chromosome
    seq_start       -- the start of the bait sequence
    seq_end         -- the end of the bait sequence
    """
    bait_seq = seq[seq_start:seq_end]
    nuc_count = get_nuc_count(bait_seq)
        
    return [
        index + 1,
        bait_chrom_name,
        seq_start,
        seq_end,
        bait_seq,
        seq_end - seq_start + 1,
        nuc_count["A"],
        nuc_count["T"],
        nuc_count["G"],
        nuc_count["C"],
        nuc_count["UNK"],
        nuc_count["AT"],
        nuc_count["GC"]
    ]
        
def process_bait_file(genome_file_path, bait_file_path, df_header, df_row_builder, output_file_path = "output.txt"):
    """Construct a data frame with columns df_header and rows according to df_row_builder.

    The data frame is written to output_file_path.

    Keyword arguments:
    genome_file_path -- the path to the genome file
    bait_file_path   -- the path to the genome file
    df_header        -- the data frame header as a list of column names
    df_row_builder   -- a function that builds a row for the data frame given
                        the index, name and sequence of a chromosome in the genome file
    output_file_path -- the path to the output file including the file name and extension
    """
    get_bait_chrom_name = lambda x: PurePath(x).stem
    bait_chrom_name = get_bait_chrom_name(bait_file_path)
    
    try:
        seq = get_chrom_seq(bait_chrom_name, genome_file_path)
    except RuntimeError as e:
        print(f"Failed to get chromosome sequence: {str(e)}")
        exit()

    df = DataFrame(columns = df_header)
    
    with open(output_file_path, "w") as output_file_handler, open(bait_file_path, "r") as bait_file_handler:
        for index, (seq_start, seq_end) in enumerate(bait_parser(bait_file_handler)):
            df.loc[index] = row_builder(index, bait_chrom_name, seq, seq_start, seq_end)

        df.to_csv(output_file_handler, sep="\t", index=None)
                
    return df

def main():
    parser = ArgumentParser(description="Extract info from the genome FASTA file according to the provided chromosome bait file.")
    parser.add_argument("fasta_file_path", type=str, help="Path to the FASTA file.")
    parser.add_argument("bait_file_path", type=str, help="Path to the bait file.")
    parser.add_argument("output_file_path", type=str, nargs="?", help="Path to the output file (including file name).", default="output.txt")
    
    args = parser.parse_args()

    process_bait_file(args.fasta_file_path, args.bait_file_path, header, row_builder, args.output_file_path)

if __name__ == "__main__":
    main()
