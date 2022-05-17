#!usr/bin/python3

import os, sys, argparse, subprocess
from datetime import datetime

# ARGUMENTS
# Create an ArgumentParser object
parser = argparse.ArgumentParser(description = 'Process raw fastq reads into VCF files.')

# Required arguments for entire pipeline
parser.add_argument('--reference', '--r', type = str, help = 'Reference fasta file.')

# Required arguments for some steps
parser.add_argument('--gbprocess', '--gbp', type = str, help = 'Path to the configuration file (.ini) for GBprocesS. Documentation available at https://gbprocess.readthedocs.io/en/latest/.')
parser.add_argument('--picard', '--p', type = str, help = 'Path to picard.jar.')
parser.add_argument('--gatk4', '--g', type = str, help = 'Path to GATK4 installation')

# Optional arguments
parser.add_argument('--threads', '--t', default = 1, type = int, help = 'number of threads to be used.')
parser.add_argument('--notrimming', action = 'store_true', help = 'skip the GBprocesS trimming step.')
parser.add_argument('--nomapping', action = 'store_true', help = 'skip the BWA MEM mapping step.')
parser.add_argument('--noaddreadgroup', action = 'store_true', help = 'skip the AddOrReplaceReadGroups (picard) step.')
parser.add_argument('--novariantcalling', action = 'store_true', help = 'skip the FreeBayes variant calling step.')
parser.add_argument('--novariantfilter', action = 'store_true', help = 'skip the GATK4 variant filtration step.')

# Parse arguments to a dictionary
args = vars(parser.parse_args())

# Check for required arguments
if args['reference'] is None:
    print("Oops, somthing went wrong...")
    sys.exit('Argument missing: --reference argument is required.')
elif not args['notrimming'] and args['gbprocess'] is None:
    print("Oops, somthing went wrong...")
    sys.exit('Argument missing: --gbprocess argument is required for the GBprocesS trimming step. Documentation available at https://gbprocess.readthedocs.io/en/latest/. Otherwise use the --notrimming argument.')
elif not args['noaddreadgroup'] and args['picard'] is None:
    print("Oops, somthing went wrong...")
    sys.exit('Argument missing: --picard argument is required for the AddOrReplaceReadGroup step. Otherwise use the --noaddreadgroup argument.')
elif not args['novariantfilter'] and args['gatk4'] is None:
    print("Oops, somthing went wrong...")
    sys.exit('Argument missing: --gatk4 argument is required for the variant filtration step. Otherwise use the --novariantfilter argument.')

# DEFINE FUNCTIONS
def print_date():
    """Print the current date and time to stderr."""
    sys.stderr.write('********************\n')
    sys.stderr.write('{}\n'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    sys.stderr.write('********************\n\n')
    return
    
# Trimming with GBprocesS
def trimming():
    """Trimming of the raw fastq files with GBprocesS (cutadapt)."""
    subprocess.run(["gbprocess", "-c", args['gbprocess'], "--debug"])

# Mapping with BWA MEM
def mapping(gz):
    """Mapping of raw fastq.gz file on reference sequence."""
    name = (os.path.splitext((os.path.splitext(gz))[0]))[0]
    ref = os.path.basename(args['reference'])
    print("Mapping " + gz + " on " + ref + "... Writing " + (name + ".sam"))
    subprocess.run(["bwa", "mem", "-t", str(args['threads']), args['reference'], gz, "-o", name + ".sam"])

# Reorganising SAM with AddOrReplaceReadGroups (picard)
def add_rg(sam):
    """Reorganising SAM file with AddOrReplaceReadGroups (picard) and writing BAM file."""
    name = (os.path.splitext(sam))[0]
    print("Reorganising " + sam + "... Writing " + name + ".recalibrated.bam")
    subprocess.run(["java", "-jar", args['picard'], "AddOrReplaceReadGroups", "-I", sam, "-O", name + ".recalibrated.bam", "--CREATE_INDEX", "-SO", "coordinate", "-RGLB", name, "-RGPL", "illumina", "-RGPU", "run", "-RGSM", name, "-RGID", name])

# Variant calling with FreeBayes
def varcal_freebayes(bam):
    """Variant calling on the BAM file with FreeBayes."""
    name = (os.path.splitext(bam))[0]
    print("Variant calling on " + bam + "... Writing VCF file")
    subprocess.run(["freebayes", "-f", args['reference'], "--genotype-qualities", "--strict-vcf", "--vcf", name + ".freebayes.vcf", bam])
    
# SNP filtering with GATK4 SelectVariants
def filter_vcf(fb_vcf):
    """Filter SNPs from the FreeBayes VCF with GATK4 SelectVariants."""
    name = (os.path.splitext(fb_vcf))[0]
    print("Filtering SNPs from " + fb_vcf + " and writing new VCF file...")
    subprocess.run([args['gatk4'], "SelectVariants", "-R", args['reference'], "-V", fb_vcf, "--select-type-to-include", "SNP", "-O", name + ".filtered.vcf"])

# ANALYSIS  
if not args['notrimming']:
    print_date()
    trimming()
      
for file in os.listdir():
    if file.endswith(".GBprocesS.FASTQ.gz"):
        if not args['nomapping']:
            print_date()
            mapping(file)
       
for file in os.listdir():
    if file.endswith(".sam"):
        if not args['noaddreadgroup']:
            print_date()
            add_rg(file)
        
for file in os.listdir():
    if file.endswith(".recalibrated.bam"):
        if not args['novariantcalling']:
            print_date()
            varcal_freebayes(file)
            
for file in os.listdir():
    if file.endswith(".recalibrated.freebayes.vcf"):
        if not args['novariantfilter']:
            print_date()
            filter_vcf(file)
