#!/usr/bin/env python

#===============================================================================
# SMAPapp_Select.py
#===============================================================================

#Yves BAWIN October 2021
#Python script to filter a SMAP haplotypes table.
#
#===============================================================================
# Import modules
#===============================================================================

import os, argparse
from datetime import datetime

#===============================================================================
# Parse arguments
#===============================================================================

#Create an ArgumentParser object
parser = argparse.ArgumentParser(description = 'Filter output file of SMAP haplotype')

'''Mandatory arguments.'''
parser.add_argument('-t', '--table',
                    type = str,
                    help = 'Name of the haplotypes table retrieved from SMAP haplotype-sites or SMAP haplotype-windows in the input directory.')

'''Input data options.'''
parser.add_argument('-i', '--input_directory',
                    default = '.',
                    type = str,
                    help = 'Input directory (default = current directory).')
parser.add_argument('-n', '--samples',
                    type = str,
                    default = None,
                    help = 'Name of a tab-delimited text file in the input directory defining the subset and order of the (new) sample names in the matrix: first column = old names, second column (optional) = new names (default = no sample list, the order of samples in the matrix equals their order in the haplotypes table).')
parser.add_argument('-l', '--loci',
                    type = str,
                    default = None,
                    help = 'Name of a tab-delimited text file in the input directory containing a subset of locus IDs formatted as in the haplotypes table in a one-column list (default = no list provided).')

'''Analysis options.'''
parser.add_argument('--keep_non_polymorphic_loci',
                    dest = 'non_poly', 
                    action = 'store_true',
                    help = 'Keep non-polymorphic loci in the filtered output file (default = non-polymorphic loci are removed).')

'''Output data options.'''
parser.add_argument('-o', '--output_directory',
                    default = '.',
                    type = str,
                    help = 'Output directory (default = current directory).')
parser.add_argument('-s', '--suffix',
                    type = str,
                    default = 'selection',
                    help = 'Suffix of the output haplotypes table (default = selection).')

#Parse arguments to a dictionary
args = vars(parser.parse_args())

#===============================================================================
# Functions
#===============================================================================
def print_date():
    print('-------------------')
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    print('-------------------\n')
    return

def txt2list(file, original, option = 'Sample', i = args['input_directory']):
    new_list = list()
    for line in open('{}/{}'.format(i, file)):
        line = line.rstrip().rstrip('\t')
        if line:
            line = line.split('\t')
            assert line[0] in original, '{} {} is not present in the haplotypes table. Please check if the {} IDs in your file specified with the --{} option match exactly with the IDs in the haplotypes table.'.format(option, line[0], option.lower(), option.lower(), 'samples' if option == 'Sample' else 'loci')
            if option == 'Sample':
                new_list.append((original.index(line[0]), line[1]))
            else:
                new_list.append(line[0])
    return new_list

def import_table(table = args['table'], i = args['input_directory'], names_list = args['samples']):
    locus_dict = dict()
    table = open('{}/{}'.format(i, table))
    samples = table.readline().rstrip().split('\t')
    indexes = range(2, len(samples))
    if names_list:
        print('\t- Selecting samples in the haplotypes table based on sample ID ...')
        sample_list = txt2list(names_list, samples)
        samples = samples[:2] + [x[1] for x in sample_list]
        indexes = [x[0] for x in sample_list]
    
    for line in table:
        line = line.rstrip().split('\t')
        calls = [line[i] for i in indexes]
        if line[0] not in locus_dict:
            locus_dict[line[0]] = {line[1]:calls}
        else:
            locus_dict[line[0]][line[1]] = calls
    return locus_dict, samples

def polymorphic(locus_dict):
    filter_dict = dict()
    for locus in locus_dict:
        poly = 0
        haplotypes = dict()
        for haplotype in locus_dict[locus]:
            h = [float(x) for x in locus_dict[locus][haplotype] if x != '<NA>']
            if len(set(h)) > 1:
                poly += 1
            if sum(h) > 0:
                haplotypes[haplotype] = locus_dict[locus][haplotype]
        if poly > 0:
            filter_dict[locus] = haplotypes
    return filter_dict

def write_output(locus_dict, samples, o = args['output_directory'], t = args['table'], s = args['suffix']):
    out_file = open('{}/{}_{}.tsv'.format(o, t[:-4], s), 'w+')
    print('\t'.join(samples), file = out_file)
    for locus, haplotypes in locus_dict.items():
        for haplotype in haplotypes.keys():
            print('{}\t{}\t{}'.format(locus, haplotype, '\t'.join(locus_dict[locus][haplotype])), file = out_file)
    return

#===============================================================================
# Script
#===============================================================================

if __name__ == '__main__':
    print_date()
    #Import haplotypes table with discrete calls. Create a subset and/or rename samples if names list is specified.
    print('* Importing haplotypes table ...')
    locus_dict, samples = import_table()
    
    if args['loci']:
        print('\t- Selecting loci in the haplotypes table based on locus ID ... ')
        loci = txt2list(args['loci'], locus_dict.keys(), option = 'Locus')
        locus_dict = {locus:locus_dict[locus] for locus in loci}
        
    #Remove non-polymorphic loci in the haplotypes table.
    if not args['non_poly']:
        print('* Removing non-polymorphic loci in the haplotypes table ...')
        locus_dict = polymorphic(locus_dict)
    
    print('* Creating output file ...')
    write_output(locus_dict, samples)
    print('* Finished!\n')
    print_date()
