#!/usr/bin/python3

import os, sys, argparse, subprocess, gzip
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

# Blast remotely against NCBI database
def DArT_blast(file)
    """Blast query FASTQ.gz file remotely against NCBI nucleotide database. Each sequence is blasted seperately and merged at the end."""
    with gzip.open(file,'r') as DArT_file:
        for line in DArT_file
