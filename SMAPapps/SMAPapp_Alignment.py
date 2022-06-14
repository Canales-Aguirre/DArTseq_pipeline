#!usr/bin/python3

#===============================================================================
# Description
#===============================================================================

# Yves BAWIN March-May 2019

#===============================================================================
# Import modules
#===============================================================================

from datetime import datetime
from natsort import natsorted
from pybedtools import BedTool
from Bio import SeqIO
import os, sys, argparse, shelve, multiprocessing, random
import pysam # make sure this is version 1.11
version = pysam.__version__
if not int(version.split('.')[1]) > 11 or (int(version.split('.')[1]) == 11 and int(version.split('.')[2]) >= 1):
    sys.exit('Please make sure Pysam version 0.11.1 (or later) is installed (current version = {})\nRun sudo pip3 install pysam==0.11.1 to fix this\n'.format(version))

#===============================================================================
# Parse arguments
#===============================================================================

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description = 'Create nucleotide alignments based on a combination of GBS and amplicon sequencing data.')

#Define arguments.
'''
Mandatory argument describing the input data.
'''
parser.add_argument('--bed',
                    type = str,
                    help = 'Comma separated string of BED files containing loci for alignment construction.')
parser.add_argument('--vcf',
                    type = str,
                    default = None,
                    help = 'Comma separated string of VCF files containing variant positions.')
parser.add_argument('--names',
                    type = str,
                    help = 'Tab delimited text file containing a table with 3-5 columns: the sample name, the corresponding GBS bam file name, the sequencing type (SE for single-end, PE for paired-end), the amplicon sequencing bam file name, and the sequencing type.')
parser.add_argument('-ref', '--reference',
                    type = str,
                    help = 'Reference genome in FastA format that was used to delineate genomic regions and for SNP calling.')

'''
Input data options.
'''
parser.add_argument('--sample_type',
                    type = str,
                    default = None,
                    help = 'Comma separated list specifying the origin of read data in the same order as the BED directories: MAS (multiplex amplicon sequencing), GBS (genotyping-by-sequencing), or other.')
parser.add_argument('--bed_dir',
                    type = str,
                    default = '.',
                    help = 'Comma-separated list of directories containing BED file (default = current directory).')
parser.add_argument('--vcf_dir',
                    type = str,
                    default = '.',
                    help = 'Comma-separated list of directories containing VCF files (default = current directory).')
parser.add_argument('--off_targets', 
                    type = str,
                    default = None,
                    help = 'Tab delimited text file containing the off-target regions and their primers (default = no file provided).')
parser.add_argument('--ploidy',
                    type = int,
                    default = 2,
                    help = 'Ploidy level of samples (default = 2).')

'''
Analysis options.
'''
parser.add_argument('--polymorphic',
                    dest = 'polymorphic_clusters_only',
                    action = 'store_true',
                    help = 'Ignore non-polymorphic clusters (default).')
parser.add_argument('--minimum_overlap',
                    type = float,
                    default = 0.6,
                    help = 'Minimum proportion of samples with data in an alignment (default = include all samples).')
parser.add_argument('--partial',
                    dest = 'ignore_partial',
                    action = 'store_false',
                    help = 'Exclude reads in multiplex amplicon sequencing data that do not completely align to the locus (default = partial reads not included).')
parser.add_argument('--no_deletions',
                    dest = 'ignore_deletions',
                    action = 'store_true',
                    help = 'Set all data on a site with a deletion in at least one haplotype to missing (default = include deletions).')
parser.add_argument('-c', '--sample_completeness',
                    type = float,
                    default = 0,
                    help = 'Minimum proportion of samples with data in an alignment (default = include all samples).')
parser.add_argument( '-sc', '--sequence_completeness',
                    type = float,
                    default = 0,
                    help = 'Minimum sequence length in proportion to the total region length. Sequence data of samples with shorter sequences is discarded (default = no filtering for sequence completeness).')
parser.add_argument('--ploidy_filtered',
                    dest = 'filter_ploidy',
                    action = 'store_true',
                    help = 'Exclude regions with more haplotypes than the defined ploidy level (default = no filtering on ploidy level).')
parser.add_argument('-r', '--min_nr_reads',
                    default = 10,
                    type = int,
                    help = 'Minimum read depth of a haplotype in a sample. Haplotypes with a lower read depth in a sample are ignored (default = 10).')
parser.add_argument('-f', '--min_frequency',
                    default = 0,
                    type = float,
                    help = 'Minimum frequency of a haplotype in a sample. Haplotypes with a lower frequency in a sample are ignored (default = 0).')
parser.add_argument('-q', '--mapping_quality',
                    default = 20,
                    type = int,
                    help = 'Minimum mapping quality score of reads. Reads with a lower mapping quality score are ignored (default = 20).')
parser.add_argument('-p', '--processes',
                    default = 4,
                    type = int,
                    help = 'Number of processes run in parallel (default = 4)')
'''
Output options.
'''
parser.add_argument('--out_dir',
                    default = '.',
                    type = str,
                    help = 'Output directory (default = current directory).')
parser.add_argument('--sequence_type',
                    choices = ['variant', 'all_sites'],
                    default = 'variant',
                    help = 'Construct sequences only containing variant sites (variant, default) or variant and invariant sites (all_sites) (default = variant).')
parser.add_argument('--alignment_type',
                    choices = ['separated', 'concatenated'],
                    default = 'separated',
                    help = 'Construct alignments per region (separated, default) or with all loci concatenated to each other (concatenated) (default = separated).')
parser.add_argument('--min_nr_samples',
                    default = 2,
                    type = int,
                    help = 'Minimum number of samples with a consensus sequence of a region. Regions with less samples with a consensus sequence are ignored (default = 2).')
parser.add_argument('--format',
                    choices = ['FastA', 'Phylip', 'Nexus'],
                    default = 'FastA',
                    help = 'Format of the alignment output file (Default = FastA).')
parser.add_argument('-o', '--out',
                    default = 'Concatenated_alignment',
                    type = str,
                    help = 'Name of the output alignment if the alignment type is concatenated (default = Alignment).')
parser.add_argument('--delete_intermediate_files',
                    dest = 'del_intermediate',
                    action = 'store_true',
                    help = 'Delete the intermediate shelve files and vcf file that were created by the script (default = intermediate files not deleted).')


#Parse arguments to a dictionary
args = vars(parser.parse_args())

#===============================================================================
# Functions
#===============================================================================
def print_date():
    '''
    Print the current date and time to stderr.
    '''
    print('-------------------')
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    print('-------------------\n')
    return


def delineate_regions(out_dir = args['out_dir'], poly_only = args['polymorphic_clusters_only']):
    '''
    Define a multi-level dictionary with region ID as key and values containing the variants.
    '''
    
    #Initiate a multi-level dictionary and open the vcf file with the SNP calls.
    regions = dict()
    
    #Initiate bedtool objects for the BED and VCF file(s) and calculate intersect between bed and vcf.
    for type in bed_dict.keys():
        vcf = open('{}/vcf_alignment.vcf'.format(out_dir))
        A = BedTool(bed_dict[type])
        B = BedTool(vcf)
        C = A.intersect(B, loj = True, nonamecheck = True)
        vcf.close()
        
        #Add loci the regions dictionary.
        for i in C:
            if poly_only and i[10] != '.':
                regions = create_region(i, regions, type)
            
    print('\tCreating consensus sequences for {} regions.\n'.format(len(regions)))
    return regions


def create_region(i, regions, type, min_overlap = args['minimum_overlap']):
    '''
    Convert the BED and VCF input into a region dictionary.
    '''
    chrom, start, stop = i[0], int(i[1]), int(i[2])
    locus = '{}-{}:{}-{}/{}'.format(type.split(':')[0], type.split(':')[1], start, stop, i[5].split(',')[0] if type.split(':')[0] == 'GBS' else '')
    
    #Check if the locus overlaps with another region in the regions dictionary.
    overlapping = list()
    for r in regions.keys():
        if chrom == r.split(':')[0]:
            rstart = int(r.split(':')[1].split('-')[0])
            rstop = int(r.split(':')[1].split('-')[1])
            if rstart < start < rstop or start < rstart < stop:
                overlapping.append(r)
    
    #Redefine the start and stop position of the combined locus.
    if any(regions[r]['type'] == 'MAS' for r in overlapping):
        overlapping_pos = natsorted([r for r in overlapping if regions[r]['type'] == 'MAS'])
        start, stop = int(overlapping_pos[0].split(':')[1].split('-')[0]), int(overlapping_pos[0].split(':')[1].split('-')[1])
    else:
        overlapping_pos = overlapping
        for r in overlapping_pos:
            if int(i[1]) > int(r.split(':')[1].split('-')[0]):
                start = int(r.split(':')[1].split('-')[0])
            if stop < int(r.split(':')[1].split('-')[1]):
                stop = int(r.split(':')[1].split('-')[1])
                
    #Define the region.
    region = '{}:{}-{}'.format(chrom, start, stop)
    
    #Add the locus to dictionary, define the locus type, and initiate the variant dictionary.
    if region not in regions:
        regions[region] = dict()
        regions[region]['variants'] = {}
        regions[region]['type'] = type.split(':')[0]
        regions[region]['loci'] = {locus}
    else:
        regions[region]['loci'].add(locus)
    
    #Add variant position and alleles.
    if i[10] != '.' and int(i[11]) - 1 not in regions[region]['variants']:
        regions[region]['variants'][int(i[11]) - 1] = {'ref' : i[13].upper(), 'alt' : i[14].upper()}
    return regions


def extract_consensus(sample, ploidy = args['ploidy'], no_del = args['ignore_deletions'], sc = args['sequence_completeness'], seq_type = args['sequence_type'], dir = args['out_dir'], ploidy_filtered = args['filter_ploidy'], min_nr_reads = args['min_nr_reads'], min_freq = args['min_frequency']):
    '''
    Create a consensus sequence for each region in the regions dictionary.
    '''
    #Check if the bam file was already parsed.
    sample_shelve = dir + '/' + sample + '.shelve'
    
    #If not, create a new shelve file and extract the consensus sequences from the bam file.
    if not os.path.isfile(sample_shelve):
        #Initiate a new shelve file to store the consensus sequences.
        shelf = shelve.open(sample_shelve, 'c')
        
        #Open the BAM files.
        GBS_samfile, MAS_samfile, other_samfile = None, None, None
        for x in samples_dict[sample]:
            if x.startswith('GBS'):
                GBS_samfile = pysam.AlignmentFile(samples_dict[sample][x], 'rb')
                GBS_type = x.split(':')[1]
            elif x.startswith('MAS'):
                MAS_samfile = pysam.AlignmentFile(samples_dict[sample][x], 'rb')
            else:
                other_samfile = pysam.AlignmentFile(samples_dict[sample][x], 'rb')
                
        #Define the abbreviations for heterozygous sites (following IUPAC convention).
        IUPAC = {'AC': 'M', 'CT': 'Y', 'CU': 'Y', 'AT': 'W', 'AU': 'W', 'AG': 'R', 'CG': 'S', 'GT': 'K', 'GU': 'K'}
        
        #For each bam file, get all haplotypes present in the read data.
        for region in natsorted(regions):
            #Define chrom, start, and stop.
            chrom, start, stop = region.split(':')[0], int(region.split(':')[1].split('-')[0]), int(region.split('-')[1])
            MAS_only, GBS_only = False, False
            
            #Define a dictionary to store all haplotypes.
            region_hts = dict()
            
            #Extract all positions in a region.
            if seq_type == 'variant':
                positions = sorted([int(p) for p in list(regions[region]['variants'].keys())])
            else:
                positions = [int(p) for p in range(start, stop)]
                
            for locus in natsorted(regions[region]['loci']):
                #Define the sequencing type, locus borders, and strand direction.
                locus_type = locus.split(':')[0]
                loc_p = [int(locus.split(':')[1].split('-')[0]), int(locus.split(':')[1].split('-')[1].split('/')[0])]
                loc_d = locus.split('/')[1]
                
                #Create a read iter containing all reads overlapping with the region sorted by the first aligned base.
                if MAS_samfile:
                    read_iter = MAS_samfile.fetch(chrom, loc_p[0], loc_p[1])
                    type = 'MAS'
                    
                    #Check amplicon sequencing reads and call haplotypes in these samples.
                    region_hts = check_reads(read_iter, region_hts, sample, type, region, locus, loc_p, positions)
                    if any(count >= min_nr_reads for count in region_hts.values()):
                        MAS_only = True
                    else:
                        region_hts = dict()
                        
                #Only check GBS reads if amplicon sequencing data are not available.
                if not MAS_only and GBS_samfile:
                    read_iter = GBS_samfile.fetch(chrom, loc_p[0], loc_p[1])
                    type = 'GBS'
                    GBS_only = True
                    #Check GBS reads and call haplotypes in these samples.
                    region_hts = check_reads(read_iter, region_hts, sample, type, region, locus, loc_p, positions, GBS_type)
                    
                #Only check other reads if GBS and amplicon sequencing data are not available.
                if not (MAS_only or GBS_only) and other_samfile:
                    read_iter = other_samfile.fetch(chrom, loc_p[0], loc_p[1])
                    type = 'other'
                    
                    #Check reads and call haplotypes in these samples.
                    region_hts = check_reads(read_iter, region_hts, sample, type, region, locus, loc_p, positions)
            
            if len(region_hts) > 0 and type != 'other':
                #Remove haplotypes with a read count lower than the predefined minimum read count and frequency.
                region_hts = {ht : region_hts[ht] for ht in region_hts if region_hts[ht] >= min_nr_reads and region_hts[ht] / sum(region_hts.values()) >= min_freq}
                
                #Remove regions with more haplotypes than the ploidy level (optional))
                if ploidy_filtered and len(region_hts) > ploidy:
                    continue
            
            #Convert positions in non stack-based sequencing data to missing data if the read count is lower than the predefined minimum.
            elif len(region_hts) > 0:
                for i, _ in enumerate(list(region_hts.keys())[0]):
                    region_hts_new = dict()
                    count = sum([region_hts[ht] for ht in region_hts if ht[i] not in ('+', '.')])
                    for ht in region_hts.keys():
                        if count < min_nr_reads:
                            ht_new = ht[:i] + '+' + ht[i + 1:]
                        else:
                            ht_new = ht
                        if ht_new not in region_hts_new:
                            region_hts_new[ht_new] = region_hts[ht]
                        else:
                            region_hts_new[ht_new] += region_hts[ht]
                    region_hts = region_hts_new
                    
            #Create a consensus haplotype.
            if len(region_hts) > 0:
                consensus = list(list(region_hts.keys())[0])
                
                for ht in region_hts.keys():
                    for i, site in enumerate(ht):
                        if site != consensus[i]:
                            if consensus[i] == '.':
                                if type == 'other':
                                    consensus[i] = site
                                    
                            #Delete site from the alignment if it was not included in any haplotype.
                            elif site == '.':
                                if type != 'other':
                                    consensus[i] = '.'
                                    
                            #Ignore all data on sites with a deletion if the no-deletions option is used.
                            elif not no_del or consensus[i] != '-':
                                if no_del and site == '-':
                                    consensus[i] == '-'
                                    
                                elif site not in ('-', '+'):
                                    #Heterozygous positions
                                    if consensus[i] in ('0', '1'):
                                        consensus[i] = '2'
                                        
                                    #Homozygous positions
                                    else:
                                        consensus[i] = site
                #Convert the consensus sequence into a nucleotide sequence.
                seq = ''
                for i, site in enumerate(consensus):
                    if site == '-':
                        seq += '-'
                    elif site in ('.', '+'):
                        seq += '?'
                    elif site == '0':
                        seq += ref[chrom][positions[i]].upper()
                    elif site == '1':
                        seq += regions[region]['variants'][positions[i]]['alt']
                    else:
                        seq += IUPAC[''.join(e for e in sorted([ref[chrom][positions[i]].upper(), regions[region]['variants'][positions[i]]['alt']]))]
                        
                #Add region to shelf if the number of nucleotides in the sequence equals or exceeds the sequence completeness.
                if seq.count('?') / len(seq) <= 1 - sc:
                    shelf['{}:{}-{}'.format(chrom, start, stop)] = seq
                    
        #Close the bam file.
        if GBS_samfile:
            GBS_samfile.close()
        if MAS_samfile:
            MAS_samfile.close()
        shelf.close()
    return sample_shelve


def check_reads(read_iter, region_hts, sample, type, region, locus, loc_p, positions, GBS_type = None, ignore_partial = args['ignore_partial'], quality_threshold = args['mapping_quality']):
    '''
    Check if reads fulfill the predefined requirements.
    '''
    #Keep track of the reads we already saw before.
    reads = dict()
    for read in read_iter:
        #Define the start of the read as the difference between reference position and alignment start.
        read_start = read.reference_start - read.query_alignment_start
        read_end = read_start + read.query_alignment_end
        
        #Skip read if the mapping quality is below the minimum.
        if int(read.mapping_quality) < quality_threshold:
            continue
            
        #Skip amplicon sequencing reads if they only partially overlap to the locus.
        elif type == 'MAS' and ignore_partial and (read.reference_start > loc_p[0] or read_end < loc_p[1]):
            continue
            
        #Skip single-end GBS reads if they are mapped into the opposite strand direction than the locus.
        elif type == 'GBS' and GBS_type == 'SE':
            if read.is_reverse:
                read_d = '-'
            else:
                read_d = '+'
            if read_d != locus.split('/')[1]:
                continue
                
        #Get CIGAR string of the read.
        cigar = read.cigartuples
        
        #Define the haplotype.
        signature = (read_start, tuple(cigar), read.query_alignment_sequence)
        
        if signature in reads:
            ht_seq = reads[signature]
        else:
            #Get a list of tuples of aligned read and reference positions (read pos, ref pos, ref seq).
            aln_pairs = read.get_aligned_pairs(with_seq = True)
            
            #Convert positions from genomic coordinates to indexes in aln_pairs.
            indexes = [get_index_from_cigar(pos, cigar, read_start) for pos in positions]
            
            #Get the corresponding aligned pairs, we put (None, index, None) if the index is out of bounds (this becomes a '.' in haplotype).
            pos_aln_pairs = [aln_pairs[i] if 0 <= i < len(aln_pairs) else (None, positions[indexes.index(i)], None) for i in indexes]
            
            #Generate haplotype string.
            ht_seq = call_haplotype(pos_aln_pairs, region, loc_p)
            
        reads[signature] = ht_seq
        
        if ht_seq != '':
            if ht_seq not in region_hts:
                region_hts[ht_seq] = 1
            else:
                region_hts[ht_seq] += 1
    return region_hts


def get_index_from_cigar(position, cigar, start):
    '''
    Transform the genomic coordinate "position" into the corresponding index in the list of aligned pairs.
    '''
    #Start genomic coordinate has index 0.
    result = position - start
    
    #Keep track of where we are in the alignment.
    i = 0
    
    #Iterate over the cigar operations.
    for operation, n in cigar:
        if operation == 1: #I (insert)
            if result >= i:
                result += n
            i += n 
        elif operation in {0,2,3,4}: #M, D, N, S
            i += n
        elif operation == 5: #H (hard clip)
            pass
        else:
            #Other operations than I, M, D, N, S, H exist, but they do not appear in the datasets
            raise NotImplementedError("CIGAR Operation " + str(operation) + " not implemented in CIGAR" + str(cigar))
    return result


def call_haplotype(pos_aln_pairs, region, loc_p):
    '''
    Transform the aligned pairs into a haplotype string.
    '''
    ht_seq = ''
    for read_pos, ref_pos, ref_seq in pos_aln_pairs:
        if ref_seq is None:
            #Soft-clipping or SMAP-related indel.
            if ref_pos is None or loc_p[0] <= int(ref_pos) < loc_p[1]:
                ht_seq += '.'
                
            #position in other part of combined locus.
            else:
                ht_seq += '+'
                
        #Deletions
        elif read_pos is None:
            ht_seq += '-'
            
        #Substitutions
        elif ref_pos in regions[region]['variants'] and ref_seq.islower():
            ht_seq += '1'
            
        #Equal to reference
        else:
            ht_seq += '0'
    return ht_seq


def create_alignment(sequences, seq_type = args['sequence_type'], poly_only = args['polymorphic_clusters_only'], c = args['sample_completeness'], sc = args['sequence_completeness'], f = args['format'], min_nr_samples = args['min_nr_samples'], dir = args['out_dir']):
    '''
    Combine the consensus sequences of all samples into an alignment dictionary.
    '''
    #Rearrange the variable sequences to create a dictionary with three levels ({region : {sample_name: sequence}}).
    seq_dict = dict()
    for sample in sequences.keys(): 
        with shelve.open('{}'.format(sequences[sample]), 'r') as shelf:
            for region, sequence in shelf.items():
                if region not in seq_dict:
                    seq_dict[region] = dict()
                seq_dict[region][sample] = list(sequence)
                
    #Filter the alignment for completeness, missing sites, and invariant sites (if variant alignment).
    alignment_dict = dict()
    for region, alignment in seq_dict.items():
        non_poly = 0
        
        #Check if the alignment contains data for the minimum required number of samples and if the proportion of samples with data equals or exceeds the completeness value.
        if len(alignment) >= min_nr_samples and len(alignment) / len(samples_dict) >= c:
            seq_length = len(list(alignment.values())[0])
            
            #Delete missing and invariant sites (if variant alignments are preferred).
            for i in range(seq_length):
                missing, sites = 0, set()
                for sequence in alignment.values():
                    if sequence[i] in ('?', '-'):
                        missing += 1
                        
                    elif sequence[i] not in 'MYWRSK':
                        sites.add(sequence[i])
                        
                if missing == len(alignment) or (seq_type == 'variant' and len(sites) < 2):
                    for sample, sequence in alignment.items():
                        seq_dict[region][sample][i] = ''
                        
                #Count the number of non-polymorphic sites in the alignment.
                elif poly_only and len(sites) < 2:
                    non_poly += 1
                    
            #Convert sequences into strings.
            if (not poly_only or non_poly < seq_length) and len(list(filter(None, list(alignment.values())[0]))) > 0:
                for sample, sequence in alignment.items():
                    if f != 'Nexus':
                        sequence = ['-' if n == '?' else n for n in sequence]
                    sequence = ''.join(e for e in sequence)
                    if not all(x in ('-', '?') for x in sequence):
                        if region not in alignment_dict:
                            alignment_dict[region] = dict()
                        alignment_dict[region][sample] = sequence
    return alignment_dict


def export(alignment, region = None, f = args['format'], out_file = args['out'], c = args['sample_completeness'], dir = args['out_dir']):
    '''
    Export the alignment dictionary into the preferred alignment format.
    '''
    
    #Create an output file.
    if region != None:
        region = [region.split(':')[0], region.split(':')[1].split('-')]
        out = open('{}/{}_{}_{}_c{}.{}'.format(dir, region[0], region[1][0], region[1][1], c, 'nex' if f == 'Nexus' else 'fa' if f == 'FastA' else 'phy'), 'w+')
    else:
        out = open('{}/{}_c{}.{}'.format(dir, out_file, c, 'nex' if f == 'Nexus' else 'fa' if f == 'FastA' else 'phy'), 'w+')
        
    #Calculate the number of samples and the length of the alignment.
    NumberOfSamples, NumberOfNucleotides = len(alignment), len(alignment[next(iter(alignment))])
    
    #Print the Phylip header to the output file if the preferred alignment format is Phylip
    if f == 'Phylip':
        print('{} {}'.format(NumberOfSamples, NumberOfNucleotides), file = out)
    
    #Print the Nexus header with the sample names to the output file if the preferred alignment format is Nexus.
    elif f == 'Nexus':
        print('#NEXUS\nBegin taxa;\n\tDimensions ntax={};\n\tformat datatype=dna missing=? gap=-;\n\ttaxlabels'.format(NumberOfSamples), file = out)
        for name, sequence in alignment.items():
            print('\t{}'.format(name), file = out)
        print(';\nend;\n\nBegin characters;\n\tDimensions nchar={};\n\tformat datatype=dna missing=? gap=-;\n\tmatrix'.format(NumberOfNucleotides), file = out)
    
    #Print the alignment
    samples_list = list(alignment.items())
    random.shuffle(samples_list)
    alignment = dict(samples_list)
    for name in alignment.keys():
        if f == 'FastA': #Default
            print('>{}\n{}'.format(name, alignment[name]), file = out)
        elif f == 'Phylip':
            print('{}\t{}'.format(name, alignment[name]), file = out)
        else:
            print('\t{}\t{}'.format(name, alignment[name]), file = out)
    
    if f == 'Nexus':
        print(';\nend;', file = out)
    return


def delete_intermediate(dir = args['out_dir']):
    for file in os.listdir(dir):
    #or file == 'vcf_alignment.vcf'
        if file.endswith('.shelve.bak') or file.endswith('.shelve.dat') or file.endswith('.shelve.dir'):
            os.remove('{}/{}'.format(dir, file))
    return


#===============================================================================
# Script
#===============================================================================
if __name__ == '__main__':
    
    print_date()
    
    #Extract bam file names from sample list.
    print('* Checking input files ...')
    samples_dict = dict()
    sample_list = open(args['names'])
    for s in sample_list:
        s = s.rstrip().split('\t')
        samples_dict[s[0]] = dict()
        for i in range(1, len(s), 2):
            samples_dict[s[0]][s[i + 1]] = s[i]
    
    #Convert info about the BED and VCF files into a dictionary.
    b_dir = args['bed_dir'].split(',')
    b_files = args['bed'].split(',')
    if args['sample_type']:
        s_types = args['sample_type'].split(',')
        for s in s_types:
            s = s.strip(' ').split(':')
            assert s[0] in ('GBS', 'MAS', 'other'), '{} is not an accepted library type. Please use GBS (genotyping-by-sequencing), MAS (multiplex amplicon sequencing), or other (for not stack-based read data).'.format(s[0])
            assert s[1] in ('PE', 'SE'), '{} is not an accepted sequencing type. Please use SE (single-end sequencing) or PE (paired-end sequencing).'.format(s[1])
        while len(b_dir) < len(b_files):
            b_dir.append(b_dir[-1])
        while len(s_types) < len(b_files):
            s_types.append(s_types[-1])
    else:
        s_types = ['MAS:PE'] * len(b_files)
        
    #Merge VCF files.
    v_dir = args['vcf_dir'].split(',')
    v_files = args['vcf'].split(',')
    while len(v_dir) < len(v_files):
        v_dir.append(v_dir[-1])
        
    #Open each vcf file and extract the SNP calls.
    SNP_dict = dict()
    for i, v in enumerate(v_files):
        file = open('{}/{}'.format(v_dir[i].strip(' ').rstrip('/'), v.strip(' ')))
        for line in file:
            if not line.startswith('#'):
                line = line.rstrip().split('\t')
                if '{}:{}'.format(line[0], line[1]) not in SNP_dict:
                    SNP_dict['{}:{}'.format(line[0], line[1])] = [line[3], line[4]]
                #Remove sites that are multi-allelic across the different vcf files.
                elif SNP_dict['{}:{}'.format(line[0], line[1])] != [line[3], line[4]]:
                    SNP_dict['{}:{}'.format(line[0], line[1])] = 'NA'
        file.close()
        
    vcf_new = open('{}/vcf_alignment.vcf'.format(args['out_dir'].rstrip('/')), 'w+')
    print('##fileformat=VCFv4.1', file = vcf_new)
    print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT', file = vcf_new)
    for SNP in natsorted(SNP_dict.keys()):
        if SNP_dict[SNP] != 'NA':
            print('{}\t{}\t.\t{}\t{}\t.\t.\t.\t.'.format(SNP.split(':')[0], SNP.split(':')[1], SNP_dict[SNP][0], SNP_dict[SNP][1]), file = vcf_new)
            
    SNP_dict.clear()
    vcf_new.close()
    
    #Create a dictionary with the bed file names as keys and the type of bed file as values.
    bed_dict = dict()
    for i, b in enumerate(b_files):
        bed_dict[s_types[i].strip(' ')] = '{}/{}'.format(b_dir[i].strip(' ').rstrip('/'), b.strip(' '))
        
    #Create a dictionary with all information about the regions.
    print('* Creating a dictionary containing all loci and their variant positions ...')
    regions = delineate_regions()
    
    #Convert the reference genome into a dictionary if a reference genome is specified.
    ref = open(args['reference'])
    ref = SeqIO.to_dict(SeqIO.parse(ref, "fasta"))
    
    #Generate alignments for each bam file.
    print('* Generating consensus sequences for each region ...')
    with multiprocessing.Pool(args['processes']) as p:
        bam_shelves = p.map(extract_consensus, list(samples_dict.keys()))
        sequences = dict(zip(list(samples_dict.keys()), bam_shelves))
        
    #Rearrange and filter the sequencing data into a dictionary with three levels ({region : {sample_name: sequence}}).
    print('* Creating {}nucleotide alignment{} ...'.format('a concatenated ' if args['alignment_type'] == 'concatenated' else '', '' if args['alignment_type'] == 'concatenated' else 's for each region'))
    alignment_dict = create_alignment(sequences)
    
    #Print the alignments to the output file.
    if args['alignment_type'] == 'concatenated':
        conc_dict = dict()
        if args['format'] == 'Nexus':
            missing = '?'
        else:
            missing = '-'
            
    print('* Exporting {}nucleotide alignment{} {}of {} regions ...'.format('the concatenated ' if args['alignment_type'] == 'concatenated' else '', '' if args['alignment_type'] == 'concatenated' else 's', 'consisting ' if args['alignment_type'] == 'concatenated' else '', len(alignment_dict)))
    for region in natsorted(alignment_dict):
        if args['alignment_type'] == 'concatenated':
            for sample in samples_dict.keys():
                if sample not in alignment_dict[region]:
                    sequence = missing * len(alignment_dict[region][list(alignment_dict[region].keys())[0]])
                else:
                    sequence = alignment_dict[region][sample]
                    
                if sample not in conc_dict:
                    conc_dict[sample] = sequence
                else:
                    conc_dict[sample] += sequence
            export(conc_dict)
        else:
            export(alignment_dict[region], region)
    if args['del_intermediate']:
        delete_intermediate()
    print('* Finished!\n')
    
    print_date()
