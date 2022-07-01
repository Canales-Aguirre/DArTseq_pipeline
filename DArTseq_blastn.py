#!/usr/bin/python3

import os, sys, argparse, subprocess, gzip
from Bio import SeqIO
from datetime import datetime

# ARGUMENTS
# Create an ArgumentParser object
parser = argparse.ArgumentParser(description = 'Blast FASTQ.gz sequences remotely against the NCBI nucleotide database.')

# Required arguments
parser.add_argument('--sample', '-s', type = str, required=True, help = 'FASTQ.gz sample to blast.')

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
    if '.gz' in fastq_file:
	    print 'Detected a .gz file.'
        name = (os.path.splitext((os.path.splitext(fastq_file))[0]))[0]
	    fastq = gzip.open(fastq_file, 'rb')
    else:
        name = (os.path.splitext(fastq_file))[0]
	    fastq = open(fastq_file, 'r')
    subprocess.run(["sed", "-n", "'1~4s/^@/>/p;2~4p'", fastq, ">", name + ".fasta"])
    
# Blast remotely against NCBI database
def DArTseq_blast(fastq_file)
    """Blast newly converted fasta file remotely against NCBI nucleotide database. Each sequence is blasted seperately and merged at the end."""
    if '.gz' in fastq_file:
        name = (os.path.splitext((os.path.splitext(fastq_file))[0]))[0] 
        fasta_file = name + ".fasta"
    else:
        name = (os.path.splitext(fastq_file))[0] + ".fasta"
        fasta_file = name + ".fasta"
    for record in SeqIO.parse(name, "fasta"):
        subprocess.run(["blastn", "-db", "nt", "-query", fasta_file, "-num_alignments", "10", "-evalue", "0.01", "-word_size", "50", "-remote", ">", name + ".blast.out"])

# ANALYSIS

fastq_to_fasta(args['sample'])
DArTseq_blast(args['sample'])
