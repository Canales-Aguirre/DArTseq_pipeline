#!usr/bin/python3

#===============================================================================
# FilterVCF.py
#===============================================================================

# Yves BAWIN March 2020
#
#===============================================================================
# Import modules
#===============================================================================

import os, sys, argparse
from datetime import datetime

#===============================================================================
# Parse arguments
#===============================================================================

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description = 'Filter a VCF file for genotype .')

'''
Mandatory argument.
'''
parser.add_argument('--vcf',
                    type = str,
                    help = 'Name of the input VCF file.')

'''
Input data options.
'''
parser.add_argument('-i', '--input_directory',
                    type = str,
                    default = '.',
                    help = 'Input directory containing the input VCF file (default = current directory).')

'''
Analysis options.
'''
parser.add_argument('-a', '--AD',
                    default = 0,
                    type = int,
                    help = 'Minimum allele depth. Genotype calls with a lower allele depth are converted into missing data (default = 0).')
parser.add_argument('-d', '--DP',
                    default = 0,
                    type = int,
                    help = 'Minimum total read depth (sum of allele depths). Genotype calls with a lower total read depth are converted into missing data (default = 0).')
parser.add_argument('-q', '--GQ',
                    default = 0,
                    type = int,
                    help = 'Minimum genotype quality. Genotype calls with a lower genotype quality are converted into missing data (default = 0).')
parser.add_argument('-m', '--MAF',
                    default = 0.0,
                    type = float,
                    help = 'Minimal minor allele frequency. SNPs with a lower minimal minor allele frequency are discarded (default = 0).')

parser.add_argument('-mh', '--maximum_proportion_heterozygous',
                    type = float,
                    default = 1.0,
                    help = 'Maximum proportion of heterozygous genotype calls (0/1) in a SNP. SNPs with a higher proportion of heterozygous genotype calls are discarded (default = 0).')
parser.add_argument('-c', '--completeness',
                    type = float,
                    default = 0.0,
                    help = 'Completeness level of SNPs. SNPs with a lower completeness level are discarded (default = 0).')

'''
Output data options.
'''
parser.add_argument('-o', '--output_directory',
                    type = str,
                    default = '.',
                    help = 'Output directory (default = current directory).')

# Parse arguments to a dictionary
args = vars(parser.parse_args())

#===============================================================================
# Functions
#===============================================================================
def print_date ():
    '''
    Print the current date and time to stderr.
    '''
    sys.stderr.write('-------------------\n')
    sys.stderr.write('{}\n'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    sys.stderr.write('-------------------\n\n')
    return


def filterGT(vcf = args['vcf'], in_dir = args['input_directory'], out_dir = args['output_directory'], a = args['AD'], d = args['DP'], q = args['GQ'], m = args['MAF']):
    '''
    Convert genotype calls with an allele depth, read depth, genotype quality, and/or minor allele frequency lower than the predefined minima into missing data.
    One output file is created: *.AD{}_DP{}_GQ{}_MAF{}.vcf
    '''
    #Extract the name of the vcf file.
    name_vcf = os.path.splitext(vcf)[0]
    
    #Open the vcf file and create a new vcf file.
    vcf = open('{}/{}'.format(in_dir, vcf))
    vcf_new = open('{}/{}.AD{}_DP{}_GQ{}_MAF{}.vcf'.format(out_dir, name_vcf, a, d, q, m), 'w+')
    
    #Iterate over lines in the vcf file and evaluate the genotype calls.
    tags = ['GT', 'AD', 'DP', 'GQ']
    total, retained = 0, 0
    for line in vcf:
        if line.startswith('#CHROM'):
            print(line, end = '', file = vcf_new)
        elif not line.startswith('##'):
            total += 1
            line = line.rstrip().split('\t')
            ref_count, alt_count = 0, 0
            #Find index of GT, AD, DP and GQ in format field
            criteria = {tag:line[8].split(':').index(tag) for tag in tags}
            line[8] = 'GT'
            
            #Iterate over genotype calls.
            for i, gt in enumerate(line[9:]):
                gt = gt.split(':')
                if len(gt) > 3:
                    #Check if the total read depth and the genotype quality are not lower than the predefined minima.
                    if int(gt[criteria['DP']]) >= d and int(gt[criteria['GQ']]) >= q:
                        #Check if the allele depths are not lower than the predefined minima:
                        AD_ref, AD_alt = gt[criteria['AD']].split(',')
                        if ('0' not in gt[criteria['GT']] or int(AD_ref) >= a) and ('1' not in gt[criteria['GT']] or int(AD_alt) >= a):
                            ref_count += gt[criteria['GT']].split('/').count('0')
                            alt_count += gt[criteria['GT']].split('/').count('1')
                            line[i + 9] = gt[criteria['GT']]
                        else:
                            line[i + 9] = './.'
                    else:
                        line[i + 9] = './.'
                        
            #Only retain sites that are still polymorphic after filtering and with a minor allele count not lower than the predefined minimum.
            if ref_count > 0 and alt_count > 0 and min(ref_count, alt_count)/(ref_count + alt_count) >= m:
                print('\t'.join(e for e in line), file = vcf_new)
                retained += 1
    sys.stderr.write('\t After filtering the genotype calls, kept {} out of a possible {} Sites\n'.format(retained, total))
    return


def filterSNP(vcf = args['vcf'], in_dir = args['input_directory'], out_dir = args['output_directory'], a = args['AD'], d = args['DP'], q = args['GQ'], m = args['MAF'], max_prop_het = args['maximum_proportion_heterozygous'], c = args['completeness']):
    '''
    Exclude SNPs with a higher proportion of heterozygous genotype calls than the predefined maximum and with a higher completeness level than the predefined minimum.
    One output file is created: *(.AD{}_DP{}_GQ{}_MAF{})_mh{}_c{}.vcf
    '''
    #Extract the name of the vcf file.
    name_vcf = os.path.splitext(vcf)[0]
    if a > 0 or d > 0 or q > 0 or m > 0:
        name_vcf += '.AD{}_DP{}_GQ{}_MAF{}'.format(a, d, q, m)
        
    #Open the vcf file and create a new vcf file.
    vcf = open('{}/{}.vcf'.format(out_dir if any(x > 0 for x in [a, d, q, m]) else in_dir, name_vcf))
    vcf_new = open('{}/{}_mh{}_c{}.vcf'.format(out_dir, name_vcf, max_prop_het, c), 'w+')
    
    #Extract the header and print it to the new vcf file.
    total, retained = 0, 0
    for line in vcf:
        if line.startswith('#CHROM'):
            print(line, end = '', file = vcf_new)
        elif not line.startswith('##'):
            total += 1
            line = line.rstrip().split('\t')
            line[8] = 'GT'
            
            #Only retain polymorphic SNP lines after filtering
            ref, het, alt, missing = 0, 0, 0, 0
            for i, genotype in enumerate(line[9:]):
                genotype = genotype.split(':')[0].split('/')
                if all(x == '0' for x in genotype):
                    ref += 1
                elif all(x == '.' for x in genotype):
                    missing += 1
                elif all(x == '1' for x in genotype):
                    alt += 1
                else:
                    het += 1
                line[9 + i] = '/'.join(e for e in genotype)
            if ref + missing < len(line[9:]) and alt + missing < len(line[9:]) and het / (het + ref + alt) <= max_prop_het and (missing / (missing + ref + alt + het)) <= (1 - c):
                print('\t'.join(e for e in line), file = vcf_new)
                retained += 1
    sys.stderr.write('\t After filtering the SNPs, kept {} out of a possible {} Sites\n'.format(retained, total))
    return


#===============================================================================
# Script
#===============================================================================

if __name__ == '__main__':
    print_date()
    if any(float(x) > 0 for x in [args['AD'], args['DP'], args['GQ'], args['MAF']]):
        sys.stderr.write('* Filtering genotype calls in {} for AD >= {}, DP >= {}, GQ >= {} and MAF >= {} ...\n'.format(args['vcf'], args['AD'], args['DP'], args['GQ'], args['MAF']))
        filterGT()
    if any(float(x) > 0 for x in [args['maximum_proportion_heterozygous'], args['completeness']]):
        sys.stderr.write('* Filtering SNPs in {} for mh >= {} and c >= {} ...\n'.format(args['vcf'], args['maximum_proportion_heterozygous'], args['completeness']))
        filterSNP()
    sys.stderr.write("* Finished! \n\n")
    print_date()
