#!/usr/bin/python3

import os, sys, argparse, subprocess, gzip, time
from Bio import SeqIO
from datetime import datetime
from pysam import FastxFile

# ARGUMENTS
# Create an ArgumentParser object
parser = argparse.ArgumentParser(description = 'Blast FASTQ.gz sequences remotely against the NCBI nucleotide database.')

# Required arguments
parser.add_argument('--sample', '-s', type = str, required=True, help = 'FASTQ.gz sample to blast.')

# Parse arguments to a dictionary
args = vars(parser.parse_args())

# Check for required arguments
if args['sample'] is None:
    print("Oops, somthing went wrong...")
    sys.exit('Argument missing: --sample argument is required.')
    
# DEFINE FUNCTIONS
def print_date():
    """Print the current date and time to stderr."""
    sys.stderr.write('                    \n')
    sys.stderr.write('********************\n')
    sys.stderr.write('{}\n'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    sys.stderr.write('********************\n\n')
    return

# Open FASTQ(.gz) and convert to fasta
def fastq_to_fasta(fastq_file):
    """Convert FASTQ.gz into basic FASTA format and blastn each sequence."""
    name = (os.path.splitext((os.path.splitext(fastq_file))[0]))[0]
    with gzip.open(fastq_file, 'rt') as fastq:
        for record in SeqIO.parse(fastq, "fastq"):
            fasta_file = open(name + ".temp.fasta", 'a')
            fasta_file.write(">"+record.id)
            fasta_file.write("\n")
            fasta_file.write(str(record.seq))
            fasta_file.write("\n")
            fasta_file.close()
            blast_file = open(name + ".blastn.out", 'a')
            subprocess.run(["blastn", "-db", "nt", "-query", name + ".temp.fasta", "-remote",  "-outfmt", "6"])
            blast_file.close()

# ANALYSIS

fastq_to_fasta(args['sample'])
