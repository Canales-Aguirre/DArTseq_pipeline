#!usr/bin/python3

#===============================================================================
# Based on Identify_barcode.py by Yves BAWIN December 2018
#===============================================================================
#
# Sander dE BACKER 2022
# Meise Botanic Garden
#
# Adapted to process FASTQ files regardless of bgzipped-state.
# List all files in directory and finds the barcode of FASTQ.gz and FASTQ files.
#
#===============================================================================
# Import modules
#===============================================================================

import os, argparse, subprocess, sys
from datetime import datetime

#===============================================================================
# Parse arguments
#===============================================================================

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description = 'Identify the most plausible barcode in FASTQ read data.')

# Mandatory arguments
parser.add_argument('-e', '--enzyme',type = str,
                    help = 'Restriction enzyme used in the GBS library preparation (included: PstI, ApeKI, MspI, MseI). \
                    If the restriction enzyme is not included, add it in the script at line 59.')

# Optional arguments
parser.add_argument('-o', '--output_directory',type = str,default = '.',
                    help = 'Output directory (default = current directory).')
                    
parser.add_argument('--symbol',type = str,default = '@',
                    help = 'Symbol at the beginning of each identifyer line in the fastq file (default = @).')
                    
# Parse arguments to a dictionary
args = vars(parser.parse_args())

#===============================================================================
# Functions
#===============================================================================
def print_date ():
    """Print the current date and time to stderr."""
    sys.stderr.write('********************\n')
    sys.stderr.write('{}\n'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + "  " + name + suffix))
    sys.stderr.write('********************\n\n')
    return


def IdentifyBC(file, enzyme = args['enzyme'], o = args['output_directory'], symbol = args['symbol']):
    """Trace the barcode in the reads of a fastQ file and create a TableOfBarcodes.fasta."""
    #Create and open the input file.
    file = open('{}{}'.format(name, suffix))
    
    #Create a TableOfBarcodes.fasta.
    TableOfBarcodes = open('{}/TableOfBarcodes.fasta'.format(o), 'a+')
    
    #Extract the cutsite of the restriction enzyme from the enzymes dictionary.
    Enzymes = {'PstI': 'TGCAG', 'ApeKI':'CWGC', 'MspI':'CGG', 'MseI':'TTAN'}
    cutsite = Enzymes[enzyme]
    
    #Create a dictionary with all potential barcodes and the number of times this barcode was found.
    barcodes = dict()
    
    #Iterate over all reads in a file and check the potential barcode sequence in the read.
    for line in file:
        if line.startswith(symbol):
            n_line = next(file)
            pos_cutsite = n_line[:25].find(cutsite)
            if pos_cutsite != -1:
                barcode = n_line[:pos_cutsite]
                if barcode not in barcodes:
                    barcodes[barcode] = 1
                else:
                    barcodes[barcode] += 1
    barcode = ''
    NumberOfOccurrences = 0
    for BC, count in barcodes.items():
        if count > 500:
            assert NumberOfOccurrences != count, 'Two barcodes occur an equal number of times: {} and {}'.format(BC, barcode)
        if count > NumberOfOccurrences:
            barcode = BC
            NumberOfOccurrences = count
    
    #Print the info to the TableOfBarcodes.
    print('>{}\n{}'.format(name, barcode), file = TableOfBarcodes)
    return


#===============================================================================
# Script
#===============================================================================

if __name__ == '__main__':
    for file in os.listdir():
        if file.endswith(".FASTQ.gz"):
            name = os.path.splitext(os.path.splitext(file)[0])[0]
            suffix = os.path.splitext(os.path.splitext(file)[0])[1]
            print_date()
            
            subprocess.run(["gunzip", file])
    
            IdentifyBC(file)
            
            subprocess.run(["bgzip", name + suffix])
            
        elif file.endswith(".FASTQ"):
            name = os.path.splitext(file)[0]
            suffix = os.path.splitext(file)[1]
            print_date()
            
            IdentifyBC(file)
            