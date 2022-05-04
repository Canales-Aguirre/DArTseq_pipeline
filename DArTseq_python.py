#/bin/python3

import os, sys, argparse, subprocess

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description = 'Process raw fastq reads into VCF files.')

parser.add_argument('--extension', '--e', type = str, help = 'Extension of files')
parser.add_argument('--reference', '--r', type = str, help = 'Reference fasta file')

# Parse arguments to a dictionary
args = vars(parser.parse_args())

# DEFINE FUNCTIONS
# Mapping with BWA MEM
def mapping(fastq):
    if args['extension'] == None:
        sys.exit('--extension argument is missing')
    elif args['reference'] == None:
        sys.exit('--reference argument is missing')
    else:
        subprocess.run(["bwa", "mem", args['reference'], fastq, ">", "${fastq%.fastq.gz}.sam"])
    
# ANALYSIS        
# Run function mapping()
for file in os.listdir():
    if file.endswith(str(args['extension'])):
        print(file)
        mapping(file)
      
#/bin/python3

import os, sys, argparse, subprocess

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description = 'Process raw fastq reads into VCF files.')

parser.add_argument('--extension', '--e', type = str, help = 'Extension of files')
parser.add_argument('--reference', '--r', type = str, help = 'Reference fasta file')

# Parse arguments to a dictionary
args = vars(parser.parse_args())

# DEFINE FUNCTIONS
# Mapping with BWA MEM
def mapping(qz):
    fastq_gz = os.path.splitext(qz)
    fastq = fastq_gz[0]
    name_fastq = os.path.splitext(fastq)
    name = name_fastq[0]
    print(name)
    if args['extension'] == None:
        sys.exit('--extension argument is missing')
    elif args['reference'] == None:
        sys.exit('--reference argument is missing')
    else:
        subprocess.run(["bwa", "mem", args['reference'], gz, "-o", basename[0]])
    
# ANALYSIS        
# Run function mapping()
for file in os.listdir():
    if file.endswith(str(args['extension'])):
        print(file)
        mapping(file)
