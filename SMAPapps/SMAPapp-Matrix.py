#!usr/bin/python3

#===============================================================================
# SMAPapp_Matrix.py
#===============================================================================

# Yves BAWIN March-May 2019

#===============================================================================
# Import modules
#===============================================================================

from datetime import datetime
from natsort import natsorted
import os, argparse, math, random
import multiprocessing as mp
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use('Agg')

#===============================================================================
# Parse arguments
#===============================================================================

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description = 'Convert the haplotype table from SMAP haplotype-sites or SMAP haplotype-windows into a genetic similarity/distance matrix and/or a locus information matrix.')

#Mandatory argument.
'''
Mandatory argument specifying the name of the haplotypes table (haplotypes_R*_F*_*.txt).
'''
parser.add_argument('-t', '--table',
					type = str,
					help = 'Name of the haplotypes table retrieved from SMAP haplotype-sites or SMAP haplotype-windows in the input directory.')

#Optional arguments.
'''
Input data options.
'''
parser.add_argument('-i', '--input_directory',
					default = '.',
					type = str,
					help = 'Input directory containing the haplotypes table, the --samples text file, and/or the --loci text file (default = current directory).')
parser.add_argument('-n', '--samples',
					type = str,
					default = None,
					help = 'Name of a tab-delimited text file in the input directory defining the order of the (new) sample IDs in the matrix: first column = old IDs, second column (optional) = new IDs (default = no list provided, the order of sample IDs in the matrix equals their order in the haplotypes table).')
parser.add_argument('-l', '--loci',
					type = str,
					default = None,
					help = 'Name of a tab-delimited text file in the input directory containing a one-column list of locus IDs formatted as in the haplotypes table (default = no list provided).')

'''
Analysis options.
'''
parser.add_argument('-lc', '--locus_completeness',
					default = 0,
					type = float,
					help = 'Minimum proportion of samples with haplotype data in a locus. Loci with less data are removed (default = all loci are included).')
parser.add_argument('-sc', '--sample_completeness',
					default = 0,
					type = int,
					help = 'Minimum number of loci with haplotype data in a sample. Samples with less data are removed (default = all samples are included).')
parser.add_argument('--include_non_shared_loci',
					dest = 'non_shared_loci',
					action = 'store_true',
					help = 'Loci with data in only one out of two samples in each comparison are included in genetic similarity and locus information calculations (default = only loci with data in both samples of each comparison are included in calculations).')
parser.add_argument('-s', '--similarity_coefficient',
					choices = ['Jaccard', 'Sorensen-Dice', 'Ochiai'],
					default = 'Jaccard',
					help = 'Coefficient used to express pairwise genetic similarity between samples (default = Jaccard, other options are: Sorensen-Dice and Ochiai).')
parser.add_argument('--distance',
					dest = 'calc_distance', 
					action = 'store_true',
					help = 'Convert genetic similarity estimates into genetic distances (default = no conversion to distances).')
parser.add_argument('--distance_method',
					choices = ['Inversed', 'Euclidean'],
					default = 'Inversed',
					help = 'Method used for genetic distance calculations (default = Inversed, other option is: Euclidean).')
parser.add_argument('-lic', '--locus_information_criterion',
					choices = ['Shared', 'Unique'],
					default = 'Shared',
					help = 'Create a matrix showing the number of loci with shared or unique haplotypes in each comparison (default = Shared (matrix showing the number of loci with shared haplotypes), other option is: Unique (matrix showing the number of loci with unique haplotypes).')
parser.add_argument('--partial',
					dest = 'partially_informative',
					action = 'store_true',
					help = 'Include loci in locus information matrix with at least one shared/unique haplotype (default = only loci with only shared/unique haplotypes are included).')
parser.add_argument('-b', '--bootstrap',
					default = 0,
					type = int,
					help = 'The number of bootstrap replicates of the genetic similarity/distance matrix (default = no bootstrap replicates).')
parser.add_argument('-p', '--processes',
					default = 4,
					type = int,
					help = 'Number of processes used by the script (default = 4).')

'''
Output data options.
'''
parser.add_argument('-o', '--output_directory',
					default = '.',
					type = str,
					help = 'Output directory (default = current directory).')
parser.add_argument('-u', '--suffix',
					default = None,
					type = str,
					help = 'Suffix added to all output file names (default = no suffix added).')
parser.add_argument('--print_sample_information',
					choices = ['Matrix', 'Plot', 'All'], 
					default = 'Matrix',
					help = 'Print the similarity/distance matrix and the number of (shared) loci to the output directory as matrixes (Matrix option, .csv file) and/or as plots (Plot option, file type specified by --plot_format option) (default = Matrix, other options are: Plot and All (Matrix and Plot).')
parser.add_argument('--print_locus_information',
					choices = [None, 'Matrix', 'Plot', 'List', 'All'],
					default = None,
					help = 'Print locus information to the output directory as a matrix (Matrix option, .csv file) and/or as a plot (Plot option, file type specified by --plot_format option) showing the number of loci with shared or unique haplotypes in each comparison. The locus information per locus can also be printed as a tab-delimited list (List option, .txt file) showing the number of comparisons wherein the locus was informative and the sample IDs with unique haplotypes in all comparisons for each locus (default = locus information is not printed, other options are: Matrix, Plot, List, and All (Matrix, Plot, and List).')
parser.add_argument('--matrix_format',
					choices = ['Phylip', 'Nexus'],
					default = 'Phylip',
					help = 'Format of the similarity/distance matrix (default = Phylip, other option is Nexus).')
parser.add_argument('--plot_format',
					choices = ['pdf', 'png', 'svg', 'jpg', 'jpeg', 'tif', 'tiff'],
					default = 'pdf',
					help = 'File format of plots (default = pdf, other options are png, svg, jpg, jpeg, tif, and tiff).')
parser.add_argument('--mask',
					choices = [None, 'Upper', 'Lower'], 
					default = None,
					help = 'Mask values on the main diagonal of each matrix and above (upper) or below (lower) the main diagonal (default = None (no masking), other options are: upper (mask upper half) and lower (mask lower half).')
parser.add_argument('--annotate_matrix_plots',
					dest = 'add_values', 
					action = 'store_true',
					help = 'Annotate the matrix plots with values (default = no annotations).')
parser.add_argument('--plot_line_curves',
					dest = 'curves',
					action = 'store_true',
					help = 'Plot line curves showing the genetic correspondence between two samples following a cumulative number of loci (default = curves are not plotted).')
parser.add_argument('--list_line_curves',
					default = None,
					type = str,
					help = 'Tab-delimited text file containing a list of sample comparisons (first column = ID of the first sample, second column = ID of the second sample) of which line curves need to be plotted (default = no list specified, all comparisons are plotted if the --plot_line_curves option is specified).')
parser.add_argument('--locus_interval',
					type = int,
					default = 10,
					help = 'Interval in the number of loci between two consecutive points in the line curves (default = 10).')
parser.add_argument('--no_matrix_plot_labels',
					dest = 'disable_plot_labels',
					action = 'store_true',
					help = 'Do not plot labels on the axes of matrix plots (default = labels are plotted on the axes of matrix plots).')
parser.add_argument('-f', '--font',
					choices = ['Times New Roman', 'Arial', 'Calibri', 'Consolas', 'Verdana', 'Helvetica', 'Comic Sans MS'],
					default = 'Times New Roman',
					help = 'Font used in all plots (default = Times_New_Roman, other options are: Arial, Calibri, Consolas, Verdana, Helvetica, and Comic_Sans_MS).')
parser.add_argument('--title_fontsize',
					type = float,
					default = 12,
					help = 'Title font size in points (default = 12).')
parser.add_argument('--label_fontsize',
					type = float,
					default = 12,
					help = 'Label font size in points (default = 12).')
parser.add_argument('--tick_fontsize',
					type = float,
					default = 8,
					help = 'Tick font size in points,(default = 8).')
parser.add_argument('--legend_fontsize',
					type = float,
					default = 10,
					help = 'Legend font size in points (default = 10).')
parser.add_argument('--legend_position',
					default = [1, 1],
					nargs = 2,
					type = float,
					help = 'Pair of coordinates defining the x (first number) and y (second number) position of the legend in the line curve plots (default = 1,1, i.e. position the legend in the top right corner of the plot).')
parser.add_argument('-r', '--plot_resolution',
					type = int,
					default = 300,
					help = 'Plot resolution in dots per inch (dpi) (default = 300).')
parser.add_argument('--colour_map',
					choices = ['viridis', 'plasma', 'inferno', 'magma', 'cividis', 'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds', 'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn', 'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink','spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia','hot', 'afmhot', 'gist_heat', 'copper', 'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic', 'twilight', 'twilight_shifted', 'hsv'],
					default = "bone",
					help = 'Colour palette used in the matrix plots (default = bone). The (Perceptually Uniform) All sequential colormaps listed on the site of MatPlotLib (https://matplotlib.org/stable/tutorials/colors/colormaps.html) are accepted.')

# Parse arguments to a dictionary
args = vars(parser.parse_args())

#===============================================================================
# Functions
#===============================================================================
def print_date ():
	'''
	Print the current date and time to stderr.
	'''
	print('-------------------')
	print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
	print('-------------------\n')
	return

def txt2dict(file, original, option = 'Sample', i = args['input_directory'], n = args['samples'], sc = args['sample_completeness']):
	'''
	Convert a tab-delimited text file into a dictionary.
	'''
	#Create an empty dictionary and set the different options, types, and lists.
	new = dict()
	options = ['Sample', 'Locus', 'Plot']
	types = ['Sample', 'Locus', 'Sample']
	lists = ['--samples', '--loci', '--list_line_curves']
	opt = options.index(option)
	
	#Iterate over the lines in the tab-delimited text file.
	for index, line in enumerate(open('{}/{}'.format(i.rstrip('/'), file))):
		#Remove trailing whitespaces and tabs. Check if any characters are left.
		line = line.rstrip().rstrip('\t')
		if line:
			line = line.split('\t')
			
			#Check if the first element in the line (i.e. original sample ID, locus ID, or reference sample in --list_line_curves) is present in the haplotypes table.
			assert line[0] in original, '{} {} on line {} in the {} list is not listed in the {}. Please check if the IDs in the {} list exactly match with those in the {}.'.format(types[opt], line[0], index + 1, lists[opt], 'haplotypes table' if not n or option == 'Sample' else '--samples list', lists[opt], 'haplotypes table' if not n or option == 'Sample' else '--samples list')
			
			#In case of --list_line_curves, check if the line contains at least two elements and if the second element (i.e. second sample ID) is present in the haplotypes table.
			if option == 'Plot':
				assert len(line) > 1, 'Line {} only contains one sample ID. Please make sure that every line in the --list_line_curves list contains two sample IDs'.format(index + 1) 
				assert line[1] in original, 'Sample {} on line {} in the --list_line_curves list is not listed in the {}. Please check if the sample IDs in the {} exactly match those in the --list_line_curves list{}.'.format(line[1], index + 1, '--samples list' if n else 'haplotypes table', '--samples list' if n else 'haplotypes table', 'and if the sample was not deleted due to a sample completeness lower than {}'.format(sc) if sc > 0 else '')
				line[1] = [line[1]]
			
			#Check if every element in the first column of the --samples list is listed only once.
			elif option == 'Sample':
				assert line[0] not in new, 'Sample {} on line {} in the --samples list is listed twice. Please make sure that each sample is only mentioned once in this list.'.format(line[0], index + 1)
			
			#Add all elements in the line to the dictionary.
			if len(line) == 1:
				line.append(None)
			if line[0] not in new:
				new[line[0]] = line[1]
			elif option == 'Plot':
				new[line[0]] += line[1]
	return new

def create_relatedness_matrix(o = args['output_directory'], pc = args['processes'], non_shared = args['non_shared_loci'], mask = args['mask'], print_sample_info = args['print_sample_information'], print_locus_info = args['print_locus_information'], lic = args['locus_information_criterion'], partial = args['partially_informative'], u = args['suffix']):
	'''
	Create a matrix based on the haplotypes table.
	'''
	#Create similarity matrixes for the estimates, the number of (shared) loci, and the locus information criterion.
	values = pd.DataFrame(index = samples, columns = samples)
	loci = pd.DataFrame(index = samples, columns = samples)
	locus_info = pd.DataFrame(index = samples, columns = samples)
	
	#Calculate genetic correspondence between each sample pair.
	with mp.Pool(pc) as p:
		results = p.map(calc_haplotype_similarity, samples)
	
	locus_dict = dict()
	#Add the calculated values to the dataframes.
	for sample_information, locus_information, plot_information in results:
		s1 = list(sample_information.keys())[0]
		for s2 in sample_information[s1]:
			if mask != 'Upper':
				values.loc[[s1], [s2]] = sample_information[s1][s2]['value']
				loci.loc[[s1], [s2]] = sample_information[s1][s2]['count']
				locus_info.loc[[s1], [s2]] = sample_information[s1][s2]['informative_count']
			if mask != 'Lower':
				values.loc[[s2], [s1]] = sample_information[s1][s2]['value']
				loci.loc[[s2], [s1]] = sample_information[s1][s2]['count']
				locus_info.loc[[s2], [s1]] = sample_information[s1][s2]['informative_count']
		
		#Structure the locus information (if preferred).
		if lic:
			for locus in locus_information:
				if locus not in locus_dict:
					locus_dict[locus] = {'count' : locus_information[locus]['count'], 'informative_count' : locus_information[locus]['informative_count'], 'samples' : {s : 0 for s in samples}}
				else:
					locus_dict[locus]['count'] += locus_information[locus]['count']
					locus_dict[locus]['informative_count'] += locus_information[locus]['informative_count']
				for sample_1, sample_2 in locus_information[locus]['comparisons']:
					locus_dict[locus]['samples'][sample_1] += 1
					locus_dict[locus]['samples'][sample_2] += 1
		
		#Add the lists with the number of loci and the genetic correspondence values of each line plot to the plot_dict.
		if len(plot_information) > 0:
			p1 = list(plot_information.keys())[0]
			for p2 in plot_information[p1]:
				plot_dict[p1][p2]['values'] = plot_information[p1][p2]['values']
				plot_dict[p1][p2]['loci'] = plot_information[p1][p2]['loci']
	
	#Write matrices to output files.
	if print_sample_info in ['Matrix', 'All']:
		print_matrix(values)
		loci.to_csv('{}/NumberOf{}Loci{}.txt'.format(o.rstrip('/'), 'Shared' if not non_shared else '', '_{}'.format(u) if u else ''), sep = '\t')
	if print_locus_info in ['Matrix', 'All']:
		locus_info.to_csv('{}/NumberOf{}{}Loci{}.txt'.format(o.rstrip('/'), 'Partially' if partial else 'Completely', lic, '_{}'.format(u) if u else ''), sep = '\t')
	return values, loci, locus_info, locus_dict, plot_dict

def calc_haplotype_similarity(s1, bs = 0, count_dict = None, s = args['similarity_coefficient'], lic = args['locus_information_criterion'], partial = args['partially_informative'], plot_int = args['locus_interval'], non_shared = args['non_shared_loci']):
	'''
	Calculate the number of loci with shared/unique haplotypes in two samples.
	'''
	#Define the first sample in the list that must be compared to s1.
	first_sample = samples.index(s1)
	
	#Create a dictionary to store the information of each comparison: estimate of genetic correspondence (value), number of assessed loci in comparison (count), and list of informative loci (informative loci).
	sample_information = {s1 : {s2 : {'value': None, 'count' : 0, 'informative_count' : 0} for s2 in samples[first_sample:]}}
	
	#Create a dictionary to store the information of each locus: number of times the locus was assessed (count), the number of times the locus was informative (informative_count), and a list with all comparisons wherein the locus was informative (comparisons).
	locus_information = {locus : {'count' : 0, 'informative_count' : 0, 'comparisons' : list()} for locus in loci}
	
	#Create a dictionary to store intermediate calculations for sample pairs in plot_dict.
	plot_information = dict()
	
	if not count_dict:
		count_dict = {locus : 1 for locus in loci}
	
	#Calculate estimates between s1 and all specified s2's.
	for s2 in samples[first_sample:]:
		#Count the number of haplotypes on sample level (all loci included): shared (a), unique for s1 (b), and unique for s2 (c).
		a, b, c, = 0, 0, 0
		
		#Check if the sample pair is in plot_dict.
		p1 = None
		if s1 in plot_dict and s2 in plot_dict[s1]:
			p1, p2 = s1, s2
		elif s2 in plot_dict and s1 in plot_dict[s2]:
			p1, p2 = s2, s1
		if p1:
			if p1 not in plot_information:
				plot_information[p1] = {p2 : {'values' : [0], 'loci' : [0]}}
			else:
				plot_information[p1][p2] = {'values' : [0], 'loci' : [0]}
		
		#Iterate over all loci in the locus_dict and infer locus similarity.
		for locus in count_dict:
			#Count the number of haplotypes on locus level within one comparison: shared (x), unique for s1 (y), and unique for s2 (z).
			x, y, z = 0, 0, 0
			
			#Extract the haplotype calls of both samples from the haplotypes table.
			s1_calls = [float(x) for x in t[locus][s1]]
			s2_calls = [float(x) for x in t[locus][s2]]
			
			#Count the number of shared and unique haplotypes in the locus if both samples have data.
			if not any(math.isnan(x) for x in s1_calls) and not any(math.isnan(x) for x in s2_calls):
				for s1_call, s2_call in zip(s1_calls, s2_calls):
					if s1_call > 0 and s2_call > 0:
						x += count_dict[locus]
					elif s1_call > 0:
						y += count_dict[locus]
					elif s2_call > 0:
						z += count_dict[locus]
			
			#Count the number of unique haplotypes in one sample if the other has not data and non-shared loci must be taken into account.
			elif non_shared:
				y += count_dict[locus] * len([call for call in s1_calls if call > 0])
				z += count_dict[locus] * len([call for call in s2_calls if call > 0])
			
			#Add a, b, and c to the sample-level counts x, y, and z.
			a += x
			b += y
			c += z
			
			#Check whether the locus was assessed.
			if x > 0 or y > 0 or z > 0:
				sample_information[s1][s2]['count'] += 1
				if bs == 0 and s1 != s2:
					locus_information[locus]['count'] += 1
					
				#Check if the locus contained shared haplotypes if a locus information matrix based on shared haplotypes must be created.
				if lic == 'Shared' and x > 0:
					if partial or (y == 0 and z == 0):
						sample_information[s1][s2]['informative_count'] += 1
						if bs == 0 and s1 != s2:
							locus_information[locus]['informative_count'] += 1
					
				#Check if the locus contained unique haplotypes if a locus information matrix based on unique haplotypes must be created
				elif lic == 'Unique' and (y > 0 or z > 0):
					if partial or x == 0:
						sample_information[s1][s2]['informative_count'] += 1
						if bs == 0 and s1 != s2:
							locus_information[locus]['informative_count'] += 1
							locus_information[locus]['comparisons'].append((s1, s2))
			
			#Add intermediate calculations to the plot_dictionary.
			if p1:
				if sample_information[s1][s2]['count'] % plot_int == 0 and plot_information[p1][p2]['loci'][-1] != sample_information[s1][s2]['count']:
					d = calc_relatedness(a, b, c)
					plot_information[p1][p2]['values'].append(d)
					plot_information[p1][p2]['loci'].append(sample_information[s1][s2]['count'])
				
		#Calculate genetic similarity/distance.
		sample_information[s1][s2]['value'] = calc_relatedness(a, b, c)
		if a == 0 and b == 0 and c == 0:
			sample_information[s1][s2]['count'] = '<NA>'
			sample_information[s1][s2]['informative_count'] = '<NA>'
			
		#Add the final calculated values to the plot_dict.
		if p1:
			if sample_information[s1][s2]['count'] != '<NA>' and sample_information[s1][s2]['count'] % plot_int != 0:
				plot_information[p1][p2]['values'].append(sample_information[s1][s2]['value'])
				plot_information[p1][p2]['loci'].append(sample_information[s1][s2]['count'])
	return (sample_information, locus_information, plot_information) if bs == 0 else sample_information

def calc_relatedness(a, b, c, s = args['similarity_coefficient'], distance = args['calc_distance'], d_method = args['distance_method']):
	'''
	Calculate pairwise genetic similarity/distance.
	'''
	#Calculate genetic similarity.
	if a > 0 or b > 0 or c > 0:
		if s == 'Jaccard':
			d = a / (a + b + c)
		elif s == 'Sorensen-Dice':
			d = (2 * a) / ((2 * a) + b + c)
		elif s == 'Ochiai':
			d = a / math.sqrt((a + b) * (a + c))
	else:
		d = '<NA>'
	
	#Optional: calculate genetic distance between two samples.
	if distance and d != '<NA>':
		if d_method == 'Euclidean':
			d = math.sqrt((1 - d)**2)
			print('{}\n'.format(d))
		else:
			d = 1 - d
	return d

def create_bootstrap_reps(bs, mask = args['mask']):
	'''
	Create bootstrap replicates of the similarity/distance matrix.
	'''
	print(' * Creating bootstrap replicate {} ...'.format(bs))
	
	#Create a dataframe to store the bootstrap estimates.
	values = pd.DataFrame(index = samples, columns = samples)
	
	#Initiate a dictionary to count the occurrence of each locus in the bootstrap replicate.
	count_dict = dict()
	
	#Randomly sample loci in the locus_dict.
	for _ in range(len(loci)):
		locus = random.choice(loci)
		if locus not in count_dict:
			count_dict[locus] = 1
		else:
			count_dict[locus] += 1
	
	#Calculate genetic correspondance between each pair of samples.
	for sample in samples:
		sample_information = calc_haplotype_similarity(sample, bs, count_dict)
		
		#Reorder sample_information into the values dataframe
		for s1, s2_dict in sample_information.items():
			for s2 in s2_dict:
				if mask != 'Upper':
					values.loc[[s1], [s2]] = sample_information[s1][s2]['value']
				if mask != 'Lower':
					values.loc[[s2], [s1]] = sample_information[s1][s2]['value']
	
	#Write matrix to output file.
	print_matrix(values, bs = bs)
	return

def set_font_style(f = args['font'], legend_fsize = args['legend_fontsize'], label_fsize = args['label_fontsize'], title_fsize = args['title_fontsize'], tick_fsize = args['tick_fontsize']):
	'''
	Set font settings for the plot title, axes titles, axes ticks, and the legend.
	'''
	font_set = dict()
	font_set['title'] = {'fontname' : f, 'fontsize' : title_fsize}
	font_set['label'] = {'fontname' : f, 'fontsize' : label_fsize}
	font_set['ticks'] = {'fontname' : f, 'fontsize' : tick_fsize}
	font_set['legend'] = {'family' : f, 'size' : legend_fsize}
	return font_set

def plot_line_curves(s1, s = args['similarity_coefficient'], d = args['calc_distance'], d_method = args['distance_method'], o = args['output_directory'], int = args['locus_interval'], plot_format =args['plot_format'], dpi = args['plot_resolution'], non_shared = args['non_shared_loci'], u = args['suffix'], lxy = args['legend_position']):
	'''
	Plot genetic relatedness values in function of the cumulative number of (shared) loci.
	'''
	#Create a matplotlib figure.
	if any(len(plot_dict[s1][x]['values']) > 1 for x in plot_dict[s1]):
		ax = plt.figure().gca()
		#Plot a line for each comparison in plot_dict.
		for s2 in natsorted(plot_dict[s1]):
			plt.plot(plot_dict[s1][s2]['loci'], plot_dict[s1][s2]['values'], linewidth = 1.0, c = colours[p_samples_list.index(s2)], linestyle = l_styles[p_samples_list.index(s2)%len(l_styles)])
		
		#Set plot title, legend, axes titles, axes limits, and axes ticks.
		plt.title('{} {} - Number of {}loci\n{}'.format(s, '{} distance'.format('inversed' if d_method == 'Inversed' else 'Euclidean') if d else 'similarity', 'shared ' if not non_shared else '', s1), **font_set['title'])
		plt.legend(natsorted(plot_dict[s1]), prop = font_set['legend'], bbox_to_anchor = lxy)
		plt.xlabel('Number of {}loci'.format('shared ' if not non_shared else ''), **font_set['label'])
		plt.ylabel('{} {}'.format(s, 'similarity' if not d else '{} distance'.format('Inversed' if d_method == 'Inversed' else 'Euclidean')), **font_set['label'])
		plt.xticks(ax.get_xticks(), **font_set['ticks'])
		plt.yticks(ax.get_yticks(), **font_set['ticks'])
		ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
		ax.set_xlim(0)
		ax.set_ylim(0, min([1, max([max(plot_dict[s1][x]['values']) for x in plot_dict[s1]]) + 0.1]))
		
		#Set tide layout.
		plt.tight_layout()
		
		#Save and close plot.
		plt.savefig('{}/{}_{}_{}_I{}{}.{}'.format(o.rstrip('/'), s, '{}_distance'.format('Inversed' if d_method == 'Inversed' else 'Euclidean') if d else 'similarity', s1, int, '_{}'.format(u) if u else '', plot_format), format = plot_format, dpi = dpi)
		plt.close()
	return

def plot_heatmap(df, option = 'values', o = args['output_directory'], dpi = args['plot_resolution'], f = args['font'], s = args['similarity_coefficient'], d = args['calc_distance'], d_method = args['distance_method'], non_shared = args['non_shared_loci'], partial = args['partially_informative'], lic = args['locus_information_criterion'], annotate = args['add_values'], masked = args['mask'], plot_format = args['plot_format'], no_labels = args['disable_plot_labels'], u = args['suffix']):
	#Set the type of the heatmap
	if option == 'values':
		if d:
			type = '{}_{}_genetic_distance_heatmap'.format(s, d_method)
		else:
			type = '{}_genetic_similarity_heatmap'.format(s)
	else:
		type = 'NumberOf{}Loci_heatmap'.format('{}'.format('Shared' if not non_shared else '') if option == 'loci' else '{}{}'.format('Partially' if partial else 'Completely', lic))
	
	#Change values in dataframe to numeric and extract values from dataframe.
	df = df.apply(pd.to_numeric, errors='coerce')
	
	#Define the minimum and maximum value of the heatmap.
	min_value, max_value = df.min().min(), df.max().max()
	
	#Create a matplotlib figure.
	fig, ax = plt.subplots(figsize=(min(len(samples) + 10, 100), min(len(samples) + 10, 100)))
	
	#Change the font of the heatmap and add a title to the heatmap.
	plt.rcParams["font.family"] = f
	plt.title(type, **font_set['title'])
	
	#Mask NaN values and define colour map. Reverse the colour code if genetic distances were calculated.
	mask = df.isnull()
	cmap = args['colour_map']
	if d and option == 'values':
		cmap += '_r'
	
	#Create heatmap with seaborn. Position the colour bar legend on the left of the figure if the lower half is masked.
	if masked != 'Lower':
		ax = sns.heatmap(df, linewidths = 0.5, cmap = cmap, vmin = min_value, vmax = max_value, mask = mask, xticklabels=True, yticklabels=True)
	else:
		ax = sns.heatmap(df, linewidths = 0.5, cmap = cmap, vmin = min_value, vmax = max_value, mask = mask, cbar_kws = dict(use_gridspec=False,location="left"), xticklabels=True, yticklabels=True)
	
	#Change label settings of the legend.
	cax = plt.gcf().axes[-1]
	cax.tick_params(labelsize = font_set['legend']['size'])
	if option != 'values':
		cax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
	else:
		cax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
	
	#Change font settings of x and y labels.
	if no_labels:
		ax.axes.get_xaxis().set_visible(False)
		ax.axes.get_yaxis().set_visible(False)
	ax.set_xticklabels(samples, fontdict = font_set['label'], rotation = 'vertical')
	ax.set_yticklabels(samples, fontdict = font_set['label'], rotation = 'horizontal')
	ax.set_facecolor('w')
	
	#Set horizontal axis to the top and vertical axis to the right if the lower half of the matrix is masked.
	if masked == 'Lower':
		ax.tick_params(axis='both', which='major', labelbottom = False, bottom = False, top = True, labeltop = True, labelleft = False, left = False, labelright = True, right = True)
		
	#Add values to the output file (optional).
	if annotate:
		#Adjust the text colour (white or black) in the heatmap depending on the chosen cmap.
		darker = ['Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds', 'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn', 'binary', 'gist_yarg']
		lighter = ['gist_gray', 'gray', 'bone', 'pink', 'hot', 'afmhot', 'gist_heat', 'inferno', 'magma', 'plasma']
		if cmap in darker or (cmap.endswith('_r') and cmap[:-2] in lighter):
			text_colours = ['k', 'w']
		elif cmap in lighter or (cmap.endswith('_r') and cmap[:-2] in darker):
			text_colours = ['w', 'k']
		else:
			text_colours = ['w', 'w']
		
		#Print the values to the heatmap.
		for y in range(len(samples)):
			for x in range(len(samples)):
				if not masked or (masked == 'Upper' and x <= y) or (masked == 'Lower' and x >= y):
					if option != 'values':
						NumberOfDigits = 0
					else:
						NumberOfDigits = 2
					if float('%.2f' % df.iat[y,x]) < min_value + (max_value - min_value) / 2:
						plt.text(x + 0.5, y + 0.5, '%.{}f'.format(NumberOfDigits) % df.iat[y, x], horizontalalignment='center', verticalalignment='center', color=text_colours[0], fontdict = font_set['ticks'])
					else:
						plt.text(x + 0.5, y + 0.5, '%.{}f'.format(NumberOfDigits) % df.iat[y, x], horizontalalignment='center', verticalalignment='center', color=text_colours[1], fontdict = font_set['ticks'])
	
	#Set tide layout if lower half of matrix is not masked.
	if masked != 'Lower':
		plt.tight_layout()
	
	#Save and close heatmap.
	fig.savefig('{}/{}{}{}.{}'.format(o.rstrip('/'), type, '_annot' if args['add_values'] else '', '_{}'.format(u) if u else '', plot_format), format = plot_format, dpi = dpi)
	plt.close()
	return

def print_matrix(df, bs = 0, s = args['similarity_coefficient'], f = args['matrix_format'], d = args['calc_distance'], d_method = args['distance_method'], o = args['output_directory'], mask = args['mask'], u = args['suffix']):
	'''
	Export similarity/distance matrix as Nexus or Phylip text file.
	'''
	if f == 'Nexus':
		#Open a new text file to print the nexus matrix and transpose the data frame if the upper half is masked.
		file = open('{}/{}_{}{}{}.nex'.format(o.rstrip('/'),'Distances' if d else 'Similarities', s, '_{}'.format(u) if u else '', '_bootstrap_{}'.format(bs) if bs > 0 else ''), 'w+')
		df = df.transpose()
		
		#Print sample information to the output file.
		print('#nexus\nBEGIN Taxa;\nDIMENSIONS ntax={};\nTAXLABELS'.format(len(samples)), file = file)
		for index, sample in enumerate(samples):
			print("[{}] '{}'".format(index + 1, sample), file = file)
		
		#Print matrix information to the output file.
		print(';\nEND; [Taxa]\nBEGIN Distances;\nDIMENSIONS ntax={};\nFORMAT labels=left diagonal triangle={};\nMATRIX'.format(len(samples), 'lower' if mask == 'Lower' else 'upper' if mask == 'Upper' else 'both'), file = file)
		for index, sample in enumerate(samples):
			print("[{}] '{}' {}".format(index + 1, sample, ' '.join(df[sample].astype(str)).lstrip()), file = file)
		print(';\nEND; [Distances]', file = file)
	else:
		df.to_csv('{}/{}{}_{}{}{}.dist'.format(o.rstrip('/'), d_method if d else '', 'Distances' if d else 'Similarities', s, '_{}'.format(u) if u else '', '_bootstrap_{}'.format(bs) if bs > 0 else ''), sep = '\t')
	return

def print_locus_info(o = args['output_directory'], t = args['table'], partial = args['partially_informative'], lic = args['locus_information_criterion'], u = args['suffix']):
	#Open a new text file to print the locus information and print a header to the output file.
	file = open('{}/{}{}Loci{}.txt'.format(o.rstrip('/'), 'Partially' if partial else 'Completely', lic, '_{}'.format(u) if u else ''), 'w+')
	print('Locus_ID\tNumberOfComparisons\tNumberOfComparisonsWith{}Haplotypes\tProportionOfComparisonsWith{}Haplotypes{}'.format(lic, lic, '\tSamplesWithUniqueHaplotypes' if lic == 'Unique' else ''), file = file)
	
	#Iterate over all loci in the locus_dict and print the info to the output file.
	for locus in natsorted(locus_dict):
		unique = natsorted([s for s, c in locus_dict[locus]['samples'].items() if c == len(samples) - 1])
		print('{}\t{}\t{}\t{}\t{}'.format(locus, locus_dict[locus]['count'], locus_dict[locus]['informative_count'], round((locus_dict[locus]['informative_count']/locus_dict[locus]['count']),2) if locus_dict[locus]['count'] > 0 else 0, ', '.join(unique)), file = file)
	return

#===============================================================================
# Script
#===============================================================================
if __name__ == '__main__':
	print_date()
	
	#Create the output directory if the directory does not exist.
	if not os.path.isdir(args['output_directory']):
		os.mkdir(args['output_directory'])

	#Import the haplotypes table (row_names = loci, column_names = samples) and extract the original sample and locus IDs.
	print(' * Loading input files ...')
	t = pd.read_csv('{}/{}'.format(args['input_directory'].rstrip('/'), args['table']), sep = '\t', dtype = str, index_col = ['Locus', 'Haplotypes'])
	samples = t.columns.tolist()
	loci = list(dict.fromkeys(t.index.get_level_values('Locus').tolist()))
	
	#Create a list with the selected (original, new) sample and locus IDs.
	if args['samples']:
		IDs = txt2dict(args['samples'], samples, 'Sample')
		
		#Pool samples if specified in the samples file. Rename samples in the haplotypes table.
		for sample, ID in IDs.items():
			if ID:
				pool = [s for s in IDs if IDs[s] == ID]
				if len(pool) > 1 and ID not in t:
					#Pool samples with the same sample ID.
					print('\t- The samples {} and {} were combined into one sample with ID {}.'.format(', '.join(pool[:-1]), pool[-1], ID))
					t[ID] = t.loc[:,pool].astype(float).sum(axis = 1, numeric_only = True, min_count = 1)
					t.drop(pool, axis = 1, inplace = True)
				else:
					#Rename samples if a new ID is specified.
					t.rename(columns={sample:ID}, inplace = True)
		
		#Reorder the columns in the dataframe based on the order of the samples in IDs.
		order = list()
		for sample, ID in IDs.items():
			if ID:
				if ID not in order:
					order.append(ID)
			else:
				order.append(sample)
		t = t[order]
		samples = t.columns.tolist()
	
	#Create a selection of loci (if specified). Remove loci that are not part of the selection and keep track of the loci that were deleted.
	if args['loci']:
		locusIDs = txt2dict(args['loci'], loci, 'Locus')
		deleted = [x for x in loci if x not in locusIDs]
		if len(deleted) > 0:
			print('\n\t- {} loc{} removed from the haplotypes table because {} not present in the --loci list:\n\n\t\t - {}\n'.format(len(deleted), 'i were' if len(deleted) > 1 else 'us was', 'it was' if len(deleted) == 1 else 'they were', '\n\t\t - '.join(deleted)))
			t.drop(deleted, level = 0, inplace = True)
		loci = list(dict.fromkeys(t.index.get_level_values('Locus').tolist()))
		assert len(loci) > 0, 'No loci were retained after the selection of loci based on the --locus_list. Please disable the --loci option or list minimum one locus ID in the --loci list.'
	
	#Remove loci with a higher percentage of missing data than allowed by the locus completeness.
	t.dropna(thresh = int(args['locus_completeness'] * len(samples)), inplace = True)
	deleted = [x for x in loci if x not in list(dict.fromkeys(t.index.get_level_values('Locus').tolist()))]
	
	#Keep track of the loci that were deleted.
	if len(deleted) > 0:
		print('\n \t- {} loc{} ignored due to a completeness lower than {}:\n\t\t- {}\n'.format(len(deleted), 'i were' if len(deleted) > 1 else 'us was', args['locus_completeness'], '\n\t\t- '.join(deleted)))
		loci = list(dict.fromkeys(t.index.get_level_values('Locus').tolist()))
		assert len(loci) > 0, 'No loci were retained after the removal of loci with a data completeness of {}. Please adjust the locus completeness.'.format(args['locus_completeness'])
	
	#Remove samples with more missing data than allowed by the sample completeness.
	samples_with_data = t.astype(float).groupby(level = 0, dropna = False).sum(numeric_only = True, min_count = 1).dropna(axis = 1, thresh = args['sample_completeness']).columns.tolist()
	deleted = [x for x in samples if x not in samples_with_data]
	t.drop(deleted, axis = 1, inplace = True)
	samples = t.columns.tolist()
	
	#Keep track of samples that were deleted.
	if len(deleted) > 0:
		print('\n \t- {} sample{} ignored because they had data in less than {} loc{}:\n\n\t\t- {}\n'.format(' {}'.format(len(deleted)) if len(deleted) > 1 else '', 's were' if len(deleted) > 1 else ' was', args['sample_completeness'], 'us' if args['sample_completeness'] == 1 else 'i', '\n\t\t- '.join(deleted)))
	assert len(samples) > 0, 'No samples were retained after the removal of samples with data in less than {} loc{}. Please adjust the sample completeness.'.format(args['sample_completeness'], 'us' if args['sample_completeness'] == 1 else 'i')
	
	#Create a list with all sample pairs that must be included in the line curve(s).
	plot_dict = dict()
	if args['curves']:
		if args['list_line_curves']:
			p_samples = txt2dict(args['list_line_curves'], samples, 'Plot')
			p_samples_list = set()
			for p1, p2s in p_samples.items():
				p_samples_list.add(p1)
				for p2 in p2s:
					p_samples_list.add(p2)
					if p1 not in plot_dict:
						plot_dict[p1] = {p2 : {'values': [0], 'loci' : [0]}}
					else:
						plot_dict[p1][p2] = {'values': [0], 'loci' : [0]}
		else:
			plot_dict = {s1 : {s2 : {'values': [0], 'loci' : [0]} for s2 in samples[samples.index(s1) + 1: ]} for s1 in samples[:-1]}
			p_samples_list = set(samples)
	
	#Set font style.
	font_set = set_font_style()
	
	#Randomly shuffle the list of loci.
	#random.shuffle(loci)
	t = t.groupby(level = 0).apply(lambda t: t.xs(t.name).to_dict(orient = 'list')).to_dict()
	
	#Create a distance matrix based on haplotypes.
	print(' * Creating matrices for {} loc{} in {} sample{} ...'.format(len(loci), 'i' if len(loci) > 1 else 'us', len(samples), '' if len(samples) == 1 else 's'))
	values_df, loci_df, locus_info_df, locus_dict, plot_dict = create_relatedness_matrix()
	
	#Create line curves.
	if args['curves']:
		print(' * Creating line curves ...')
		#Define colours and line styles for line curves.
		colours = sns.color_palette('husl', n_colors=len(p_samples_list))
		l_styles = ['solid', 'dashed', 'dashdot', 'dotted']
		p_samples_list = natsorted(list(p_samples_list))
		
		#Create a plot for each samples in the first column of the line_curves list.
		with mp.Pool(args['processes']) as p:
			p.map(plot_line_curves, plot_dict.keys())
	
	#Print locus information.
	if args['print_locus_information'] in ['List', 'All']:
		print(' * Printing locus information as a list to the output directory ...')
		print_locus_info()
	
	#Create heat maps for the genetic relatedness values and the number of common loci.
	if args['print_sample_information'] in ['Plot', 'All']:
		print(' * Creating heatmaps showing the {} {} and the number of {}loci between each pair of samples ...'.format(args['similarity_coefficient'] if not args['calc_distance'] else args['distance_method'], 'similarity' if not args['calc_distance'] else 'distance', 'shared ' if not args['non_shared_loci'] else ''))
		plot_heatmap(values_df)
		plot_heatmap(loci_df, 'loci')
	if args['print_locus_information'] in ['Plot', 'All']:
		plot_heatmap(locus_info_df, 'locus_info')
		
	#Create a number of bootstrap replicates.
	if args['bootstrap'] > 0:
		with mp.Pool(args['processes']) as p:
			p.map(create_bootstrap_reps, range(1, args['bootstrap'] + 1))
			
	print(' * Finished!\n')
	print_date()
