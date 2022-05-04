#/bin/python3

import os, argparse, subprocess

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description = 'Process raw fastq reads into VCF files.')

parser.add_argument('--extension', '--e', type = str, help = 'Extension of files')
parser.add_argument('--reference', '--r', type = str, help = 'Reference fasta file')

# Parse arguments to a dictionary
args = vars(parser.parse_args())

for file in os.listdir():
    if file.endswith(args['extension']):
        print(file)
        print(len(file))
        subprocess.run(["bwa", "mem", args['reference'], file])
