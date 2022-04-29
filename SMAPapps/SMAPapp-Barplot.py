#!/usr/bin/env python

#===============================================================================
# SMAPapp_Barplot.py
#===============================================================================

#Yves BAWIN October 2021
#Python script to create bar plots based on a SMAP haplotypes table.
#
#===============================================================================
# Import modules
#===============================================================================

import os, argparse
from datetime import datetime
import pandas as pd
import numpy as np
import multiprocessing as mp
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use('Agg')

#===============================================================================
# Parse arguments
#===============================================================================

#Create an ArgumentParser object
parser = argparse.ArgumentParser(description = 'Create bar plots showing the number of loci with x haplotypes or the similarity of samples with reference samples.')

'''Mandatory arguments.'''
parser.add_argument('-t', '--table',
                    type = str,
                    help = 'Name of the haplotypes table retrieved from SMAP haplotype-sites or SMAP haplotype-windows in the input directory.')

'''Input data options.'''
parser.add_argument('-i', '--input_directory',
                    default = '.',
                    type = str,
                    help = 'Input directory containing the haplotypes table, the --samples text file, the --loci text file, --reference_samples text file, and/or the --reference_samples_table file (default = current directory).')
parser.add_argument('-n', '--samples',
                    type = str,
                    default = None,
                    help = 'Name of a tab-delimited text file in the input directory defining the order of the (new) sample names in the barplot: first column = old names, second column (optional) = new names (default = no sample list, the order of samples in the bar plot equals their order in the haplotypes table).')
parser.add_argument('-l', '--loci',
                    type = str,
                    default = None,
                    help = 'Name of a tab-delimited text file in the input directory containing a one-column list of locus IDs formatted as in the haplotypes table (default = no list provided).')
parser.add_argument('-r', '--reference_samples',
                    type = str,
                    default = None,
                    help = 'Name of a tab-delimited text file in the input directory listing the (new) IDs of samples used as references in the bar plot: first column = sample name, second column (optional): colour ID (default = no list with reference samples IDs is not provided).')
parser.add_argument('-a', '--additional_table',
                    type = str,
                    default = None,
                    help = 'Name of an additional haplotypes table retrieved from SMAP haplotype-sites or SMAP haplotype-windows in the input directory. (default = no additional haplotypes table provided).')

'''Analysis options.'''
parser.add_argument('--type',
                    choices = ['haplotype_count', 'reference_similarity'],
                    default = 'haplotype_count',
                    help = 'Type of the categories in the bar plot. Two options are allowed: haplotype_count (i.e. the number of loci with x haplotypes, with x between one and the ploidy level plus one) and reference_similarity (the number of loci with x haplotypes equal to one of the specified reference samples (default = haplotype_count).')
parser.add_argument('--ploidy',
                    type = int,
                    default = 2,
                    help = 'Integer defining the (highest) ploidy level of the samples in the haplotypes table (default = 2, diploid).')
parser.add_argument('--value_type',
                    choices = ['count', 'mean', 'median', 'percentage', 'minimum', 'minimum_without_zero', 'maximum'],
                    default = 'count',
                    help = 'Type of values displayed in the bar plots. Seven options are allowed: count, mean, median, percentage, minimum, minimum_without_zero, and maximum (default = count).')
parser.add_argument('--unique_haplotypes',
                    dest = 'unique',
                    action = 'store_true',
                    help = 'Only use haplotypes that are unique to one of the reference samples to calculate the reference similarity (default = all haplotypes (unique and non-unique) are used to calculate the reference similarity).')
parser.add_argument('-p', '--processes',
                    default = 4,
                    type = int,
                    help = 'Number of processes used by the script (default = 4).')

'''Output data options.'''
parser.add_argument('-o', '--output_directory',
                    default = '.',
                    type = str,
                    help = 'Output directory (default = current directory).')
parser.add_argument('-s', '--suffix',
                    default = None,
                    type = str,
                    help = 'Suffix added to all output file names (default = no suffix added).')
parser.add_argument('-b', '--bar_width',
                    default = 0.3,
                    type = float,
                    help = 'Float number defining the width of each bar in the bar plot (default = 0.3).')
parser.add_argument('--plot_type',
                    choices = ['unstacked', 'stacked'],
                    default = 'unstacked',
                    help = 'Bar plot type. Two options are allowed: unstacked and stacked (default = unstacked).')
parser.add_argument('--plot_direction',
                    choices = ['vertical', 'horizontal'],
                    default = 'vertical',
                    help = 'Direction of the bar plot. Two options are allowed: vertical and horizontal (default = vertical).')
parser.add_argument('--plot_reference',
                    dest = 'plot_refs',
                    action = 'store_true',
                    help = 'Include reference samples in the barplot (default = reference samples are not included).')
parser.add_argument('--plot_format',
                    choices = ['pdf', 'png', 'svg', 'jpg', 'jpeg', 'tif', 'tiff'],
                    default = "pdf",
                    help = 'File format of plots (default = pdf).')
parser.add_argument('-f', '--font',
                    choices = ['Times New Roman', 'Arial', 'Calibri', 'Consolas', 'Verdana', 'Helvetica', 'Comic Sans MS'],
                    default = 'Times New Roman',
                    help = 'Font used in all plots (default = Times_New_Roman, other options are: Arial, Calibri, Consolas, Verdana, Helvetica, and Comic_Sans_MS).')
parser.add_argument('--title_fontsize',
                    type = int,
                    default = 12,
                    help = 'Title font size in points (default = 12).')
parser.add_argument('--label_fontsize',
                    type = int,
                    default = 12,
                    help = 'Label font size in points (default = 12).')
parser.add_argument('--tick_fontsize',
                    type = int,
                    default = 8,
                    help = 'Tick font size in points,(default = 8).')
parser.add_argument('--legend_fontsize',
                    type = int,
                    default = 10,
                    help = 'Legend font size in points (default = 10).')
parser.add_argument('--legend_position',
                    default = [1, 1],
                    nargs = 2,
                    type = float,
                    help = 'Pair of coordinates defining the x (first number) and y (second number) position of the legend (default = 1,1, i.e. position the legend in the top right corner of the plot).')
parser.add_argument('--plot_resolution',
                    type = int,
                    default = 300,
                    help = 'Plot resolution in dots per inch (dpi) (default = 300).')
parser.add_argument('--colours',
                    type = str,
                    default = '#1F77B4, #FF7F0E, #2CA02C, #D62728, #9467BD, #8C564B, #E377C2, #7F7F7F, #BCBD22, #17BECF',
                    help = 'Comma-separated list of colours used in the barplot (default = category10 color palette, more info on: https://github.com/d3/d3-3.x-api-reference/blob/master/Ordinal-Scales.md#category10).)')
parser.add_argument('--print_sample_information',
                    dest = 'print_sample_info',
                    action = 'store_true',
                    help = 'Print a table with the haplotype count or the similarity to the reference samples (columns) for each sample (rows) to a tab-delimited text file (default = sample information is not printed).')
parser.add_argument('--print_locus_information',
                    dest = 'print_locus_info',
                    action = 'store_true',
                    help = 'Print a table with the haplotype count or the similarity to the reference samples (columns) for each locus in each sample (rows) to a tab-delimited text file (default = locus information is not printed).')
parser.add_argument('--plot_heatmap',
                    dest = 'plot_heatmap_ref_samples',
                    action = 'store_true',
                    help = 'Plot a heatmap displaying the haplotype count or reference similarity on all loci (rows) per sample (columns) (default = no heatmap is plotted).')
parser.add_argument('--no_heatmap_labels',
                    dest = 'disable_plot_labels',
                    action = 'store_true',
                    help = 'Do not plot labels on heatmap axes (default = labels are plotted on the axes of heatmap).')
parser.add_argument('--annotate_heatmap',
                    dest = 'add_values', 
                    action = 'store_true',
                    help = 'Annotate heatmap plots with haplotype counts or reference similarities (default = no heatmap annotations).')
parser.add_argument('--colour_map',
                    choices = ['viridis', 'plasma', 'inferno', 'magma', 'cividis', 'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds', 'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn', 'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink','spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia','hot', 'afmhot', 'gist_heat', 'copper', 'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic', 'twilight', 'twilight_shifted', 'hsv'],
                    default = "bone",
                    help = 'Colour palette used in the heatmap plots (default = bone). The (Perceptually Uniform) All sequential colormaps listed on the site of MatPlotLib (https://matplotlib.org/stable/tutorials/colors/colormaps.html) are accepted.')

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

def import_table(table, i = args['input_directory']):
    return pd.read_csv('{}/{}'.format(args['input_directory'].rstrip('/'), table), sep = '\t', dtype = str, index_col = ['Locus', 'Haplotypes'])

def txt2dict(original, file = args['samples'], opt = 'Sample', i = args['input_directory'], n = args['samples']):
    new = dict()
    lists = {'Sample' : '--samples', 'Locus' : '--loci', 'Ref_sample' : '--reference_samples'}
    for index, line in enumerate(open('{}/{}'.format(i.rstrip('/'), file))):
        line = line.rstrip().rstrip('\t')
        if line:
            line = line.split('\t')
            assert line[0] in original, '{} {} on line {} in the {} list is not listed in the {}. Please check if the IDs in the {} list exactly match with those in the {}.'.format(opt if opt != 'Ref_sample' else 'Sample', line[0], index + 1, lists[opt], 'haplotypes table' if not n or opt == 'Sample' else '--samples list', lists[opt], 'haplotypes table' if not n or opt == 'Sample' else '--samples list')
            if opt != 'Locus':
                assert line[0] not in new, 'Sample {} on line {} in the {} list is listed twice. Please make sure that each sample is only mentioned once in this list.'.format(line[0], index + 1, lists[opt])
                if opt == 'Ref_sample' and len(line) > 1:
                    assert line[1] not in new.values(), 'Sample {} on line {} in the --reference_samples list has the same colour as sample {}. Please specify a unique colour for each reference sample.'.format(line[0], index + 1, [k for k, v in new if v == line[1]][0])
            if len(line) == 1:
                line.append(None)
            new[line[0]] = line[1]
    return new

def rename_samples(t):
    IDs = txt2dict(samples)
    for sample in samples:
        #Remove samples from the haplotypes table that are not in the sampleIDs list.
        if sample not in IDs:
            t.drop(sample, axis = 1, inplace = True)
        elif IDs[sample]:
            pool = [ID for ID in IDs if IDs[ID] == IDs[sample]]
            if len(pool) > 1 and IDs[sample] not in t:
                #Pool samples with the same sample ID.
                print('\n\t- The samples {} and {} were combined into one sample with ID {}.\n'.format(', '.join(pool[:-1]), pool[-1], IDs[sample]))
                t[IDs[sample]] = t.loc[:,pool].astype(float).max(axis = 1, numeric_only = True)
                t.drop(pool, axis = 1, inplace = True)
            else:
                #Rename samples if a new ID is specified.
                t.rename(columns={sample:IDs[sample]}, inplace = True)
    
    #Reorder the columns in the dataframe based on the order of the samples in IDs.
    order = list()
    for old, new in IDs.items():
        if new:
            if new not in order:
                order.append(new)
        else:
            order.append(old)
    t = t[order]
    return t

def get_colour_IDs(ref_samples, colours = args['colours']):
    #Define a list of colours.
    colour_list = colours.split(',')
    
    #Iterate over the reference samples.
    for index, sample in enumerate(ref_samples.keys()):
        if ref_samples[sample]:
            #Check if the colour name is in the colour list.
            if ref_samples[sample].lower() in mpl.colors.CSS4_COLORS:
                colour = ref_samples[sample].lower()
            
            #Get the name of colours specified in Hexadecimal code (HEX).
            elif ref_samples[sample].upper() in mpl.colors.CSS4_COLORS.values() or ref_samples[sample].upper() in colours:
                colour = ref_samples[sample].upper()
                
            #Get the name from a colour specified in Red - Green - Blue code (RGB).
            else:
                colour = tuple([int(x.strip('()')) for x in ref_samples[sample].split(',')])
                assert len(colour) == 3, 'The colour {} of sample {} in the --reference_sample list is not valid. Please provide a CSS colour name as listed on https://matplotlib.org/stable/gallery/color/named_colors.html, a hexadecimal (HEX) colour code starting with a #, or a rgb colour code consisting of 3 numbers ranging between 0 and 255 that are separated by a comma (e.g. 0, 0, 0).'.format(ref_samples[sample], sample)
        else:
            colour = colour_list[index % len(colour_list)]
            
        #Add the colour code to the ref_samples dictionary.
        ref_samples[sample] = colour
    
    #Add an unassigned and unknown category to the reference samples. If already present, the names of these categories are capitalized.
    if 'unassigned' not in ref_samples and 'Unassigned' not in ref_samples:
        ref_samples['Unassigned'] = 'silver'
    elif 'unassigned' in ref_samples:
        ref_samples['Unassigned'] = ref_samples.pop('unassigned')
    if 'unknown' not in ref_samples and 'Unknown' not in ref_samples:
        ref_samples['Unknown'] = 'grey'
    elif 'unknown' in ref_samples:
        ref_samples['Unknown'] = ref_samples.pop('unknown')
    return ref_samples

def check_samples(sample, type = args['type'], ploidy = args['ploidy'], unique_haplotypes = args['unique']):
    #Create the list 'counts' to store the count data.
    count_df = pd.DataFrame(index = indexes, columns = columns)
    count_df = pd.concat({sample : count_df}, names=['Sample'])
    
    #Iterate over all loci in the haplotypes table.
    for locus in loci:
        calls = t.loc[locus, sample].fillna(0).astype(float)
        
        #Only consider loci without missing data.
        if sum(calls) > 0:
            NumberOfHaplotypes = 0
            
            #Extract the reference calls for the considered locus. Remove columns with only missing data and replace all 0.0 values in the reference samples by missing data.
            ref_calls = t.loc[locus, [x for x in ref_IDs if x not in ['Unassigned', 'Unknown']]].dropna(axis = 1, how = 'all').astype(float).replace(0.0, np.nan)
            
            if type == 'haplotype_count' or len(ref_calls.columns) == len(ref_samples) - 2:
                #Add two additional columns to the dataframe for the unknown and the unassigned category and rearrange the order of the dataframe based on the ref_IDs.
                if type == 'reference_similarity':
                    ref_calls['Unassigned'] = np.nan
                    ref_calls['Unknown'] = np.nan
                    ref_calls = ref_calls[ref_IDs]
                
                #Iterate over the calls in the locus and call the number of haplotypes in the sample. Only check for the haplotype in the reference samples if the call is larger than 0.
                for index, call in calls.iteritems():
                    if not unique_haplotypes or ref_calls.loc[index].sum() == 1:
                        count = ref_calls.loc[index].tolist()
                        if call > 0:
                            NumberOfHaplotypes += 1
                            
                             #Check the number of times the haplotype is present in the reference samples. Store the count data in the counts list.
                            if type == 'reference_similarity':
                                count = check_unassigned_unknown(call, count)
                                count_df.loc[(sample, locus, index)] = count
                        elif type == 'reference_similarity':
                            if pd.isna(call):
                                count_df.loc[(sample, locus, index)] = [np.nan] * len(ref_IDs)
                            else:
                                 count_df.loc[(sample, locus, index)] = [x if np.isnan(x) else call for x in count]
                
                #Add the number of haplotypes to the counts list. Change the number of haplotypes to 'more than the ploidy level' if the NumberOfHaplotypes is larger than the ploidy level.
                if type == 'haplotype_count':
                    count_df.loc[(sample, locus)] = [0] * (ploidy + 2)
                    count_df.loc[(sample, locus)].iloc[min(NumberOfHaplotypes, ploidy + 1)] = 1
                
            #Add a list of missing values to the counts list if at least one reference sample had no data.
            else:
                for index, _ in calls.iteritems():
                    count_df.loc[(sample, locus, index)] = [np.nan] * len(ref_IDs)
                
        #Add a list of missing values to the counts list if the locus had no data.
        else:
            if type == 'reference_similarity':
                for index, _ in calls.iteritems():
                    count_df.loc[(sample, locus, index)] = [np.nan] * len(ref_IDs)
            else:
                count_df.loc[(sample, locus)] = [np.nan] * (ploidy + 2)
    return count_df

def check_unassigned_unknown(call, count):
    #Add the haplotype count/frequency to the unknown category if the haplotype is not present in any reference sample.
    if np.nansum(count) == 0:
        count_index = ref_IDs.index('Unknown')
    
    #Add the haplotype count/frequency to the unassigned category if the haplotype is present in more than one reference sample.
    elif len([x for x in count if x > 0]) > 1:
        count_index = ref_IDs.index('Unassigned')
    
    #In all other cases, convert the pandas series into a list and replace missing values by 0's.
    else:
        count_index = count.index(np.nansum(count))
    count[count_index] = call
    return count

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

def plot_barplot(count_df, table = args['table'], type = args['type'], bar_width = args['bar_width'], o = args['output_directory'], format = args['plot_format'],plot_type = args['plot_type'], plot_direction = args['plot_direction'], value_type = args['value_type'], lxy = args['legend_position'], dpi = args['plot_resolution'], colours = args['colours'], s  = args['suffix'], unique_haplotypes = args['unique'], print_sample_info = args['print_sample_info']):
    '''
    Make bar plots of per sample for the calculated values.
    '''
    #Set plot type.
    if plot_type == 'stacked':
        stacked = True
    else:
        stacked = False
    
    #Get bar colours and create bar plot.
    if type == 'haplotype_count':
        colours = [x.strip(' ').lower() for x in colours.split(',')]
    else:
        colours = ref_samples
    
    #Calculate the sum, the mean, or the percentage for each sample in the dataframe.
    stderr = None
    if value_type == 'mean':
        values = count_df.groupby(level = 0).mean()
        stderr = count_df.groupby(level = 0).sem()
    elif value_type == 'percentage':
        values = count_df.astype(float).groupby(level = 0).apply(lambda x: 100 * x / x.sum().sum()).groupby(level = 0).sum()
    elif value_type == 'median':
        values = count_df.groupby(level = 0).median()
    elif value_type == 'minimum':
        values = count_df.dropna(how = 'all').groupby(level = 0).min()
    elif value_type == 'minimum_without_zero':
        values = count_df.replace(0.0, np.nan).dropna(how = 'all').groupby(level = 0).min()
    elif value_type == 'maximum':
        values = count_df.groupby(level = 0).max()
    else:
        values = count_df.groupby(level = 0).sum()
    values = values.reindex(samples)
    
    #Print sample information (i.e. haplotype count or similarity to reference samples for each sample) to a tab-delimited output file.
    if print_sample_info:
        values.to_csv('{}/{}_{}_samples{}.txt'.format(o, type.capitalize(), value_type if value_type != 'count' else 'number', '_{}'.format(s) if s else ''), sep = '\t')
    
    if plot_direction == 'horizontal':
        #Reverse the order of the samples in the dataframe in horizontal bar plots.
        values = values.iloc[::-1]
        
        #Reverse the order of the haplotype counts in horizontal, unstacked bar plots.
        if plot_type == 'unstacked':
            values = values[values.columns[::-1]]
            colours = colours[::-1]
        ax = values.plot.barh(stacked = stacked, width = bar_width, xerr = stderr, color = colours, figsize = (len(samples), 10))
        
        #Specify x and y-axis labels.
        x, y = ax.get_xticks(), samples[::-1]
        ax.set_xlabel('{}'.format(value_type.capitalize() if value_type != 'count' else 'Number of loci'), **font_set['label'])
        ax.set_ylabel('')
    else:
        ax = values.plot.bar(stacked = stacked, width = bar_width, yerr = stderr, color = colours, figsize = (len(samples), 10))
        
        #Specify x and y-axis labels.
        x, y = samples ,ax.get_yticks()
        ax.set_xlabel('')
        ax.set_ylabel('{}'.format(value_type.capitalize() if value_type != 'count' else 'Number of loci'), **font_set['label'])
    
    #Set legend properties. Reverse the order of the legend in horizontal, unstacked bar plots.
    handles, labels = ax.get_legend_handles_labels()
    if plot_direction == 'horizontal' and plot_type == 'unstacked':
        handles, labels = handles[::-1], labels[::-1]
    ax.legend(handles, labels, bbox_to_anchor = lxy, prop = font_set['legend'])
    
    #Set title properties.
    ax.set_title('{} per sample'.format(value_type.capitalize() if value_type != 'count' else 'Number of loci'), loc = 'left', **font_set['title'])
    
    #Set axes properties.
    if value_type == 'count':
        if plot_direction == 'horizontal':
            x = [int(n) for n in x]
        else:
            y = [int(n) for n in y]
    else:
        if plot_direction == 'horizontal':
            x = [round(n, 2) for n in x]
        else:
            y = [round(n, 2) for n in y]
    
    ax.xaxis.set_ticks(ax.get_xticks())
    ax.set_xticklabels(x, **font_set['ticks'])
    ax.yaxis.set_ticks(ax.get_yticks())
    ax.set_yticklabels(y, **font_set['ticks'])
    
    #Set the values ticks to integers if value_type == number of to a float with two digits if value_type == mean or value_type == percentage.
    if plot_direction == 'horizontal':
        if value_type == 'count':
            ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
        else:
            ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
            if value_type == 'percentage':
                ax.set_xlim(0, 101)
    else:
        if value_type == 'count':
            ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
        else:
            ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
            if value_type == 'percentage':
                ax.set_ylim(0, 101)
    
    #Remove top and right edge from the MatPlotLib frame and make sure that all figure elements fit within the figure size.
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    plt.tight_layout()
    
    #Save and close figure
    plt.savefig('{}/BarPlot_{}_{}_{}_{}{}.{}'.format(o.rstrip('/'), type, value_type if value_type != 'count' else 'number', plot_type, plot_direction, '_{}'.format(s) if s else '', format), dpi = dpi, format = format)
    plt.close()
    return

def plot_heatmap(df, o = args['output_directory'], type = args['type'], dpi = args['plot_resolution'], f = args['font'], annotate = args['add_values'], cmap = args['colour_map'], plot_format = args['plot_format'], no_labels = args['disable_plot_labels'], s = args['suffix']):
    #Change values in dataframe to numeric and extract values from dataframe.
    df = df.apply(pd.to_numeric, errors='coerce')
    
    #Define the minimum and maximum value of the heatmap.
    min_value, max_value = df.min().min(), df.max().max()
    
    #Create a matplotlib figure.
    fig, ax = plt.subplots(figsize = (min(df.shape[0] + 10, 100), min(df.shape[1] + 10, 100)))
    
    #Change the font of the heatmap and add a title to the heatmap.
    plt.rcParams["font.family"] = f
    plt.title('Heatmap reference similarities', **font_set['title'])
    
    #Create heatmap with seaborn.
    ax = sns.heatmap(df, cmap = cmap, vmin = min_value, vmax = max_value, xticklabels=True, yticklabels=True)
    
    #Change label settings of the legend.
    cax = plt.gcf().axes[-1]
    cax.tick_params(labelsize = font_set['legend']['size'])
    cax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
    
    #Change font settings of x and y labels.
    if no_labels:
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
    ax.set_xticklabels(samples, fontdict = font_set['label'], rotation = 'vertical')
    ax.set_yticklabels(df.index.get_level_values(1), fontdict = font_set['label'], rotation = 'horizontal')
    ax.set_facecolor('w')
    
    #Remove axis titles.
    ax.set_ylabel('')
    ax.set_xlabel('')
    
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
        for y in range(df.shape[0]):
            for x in range(df.shape[1]):
                if float('%.2f' % df.iat[y,x]) < min_value + (max_value - min_value) / 2:
                    plt.text(x + 0.5, y + 0.5, '%.2f' % df.iat[y, x], horizontalalignment='center', verticalalignment='center', color=text_colours[0], fontdict = font_set['ticks'])
                else:
                    plt.text(x + 0.5, y + 0.5, '%.2f' % df.iat[y, x], horizontalalignment='center', verticalalignment='center', color=text_colours[1], fontdict = font_set['ticks'])
    
    #Set tide layout.
    plt.tight_layout()
    
    #Save and close heatmap.
    fig.savefig('{}/Heatmap_{}{}{}.{}'.format(o.rstrip('/'), type, '_annot' if args['add_values'] else '', '_{}'.format(s) if s else '', plot_format), format = plot_format, dpi = dpi)
    plt.close()
    return

#===============================================================================
# Script
#===============================================================================
if __name__ == '__main__':
    print_date()
    
    print(' * Loading input files ...')
    #Import master table with haplotype calls and extract sample names.
    t = import_table(args['table'])
    
    if args['additional_table']:
        #Import the additional haplotypes table.
        t2 = import_table(args['additional_table'])
        
        #Delete columns in the additional haplotypes table that are also present in the first table.
        deleted = t2[t2.columns.intersection(t.columns)].columns.tolist()
        t2.drop(deleted, axis = 1, inplace = True)
        if len(deleted) > 0:
            print('\n\t- {} samples were removed from the additional haplotypes table because they were present in the first haplotypes table:\n\t\t- {}'.format(len(deleted), '\n\t\t- '.join(deleted)))
        
        #Merge the two tables.
        t = t.merge(t2, left_index = True, right_index = True, how = 'outer')
    
    #Get (and rename if specified) sample IDs from master table.
    samples = t.columns.tolist()
    if args['samples']:
        #Pool and/or rename samples if specified in the --samples list. Update the samples list.
        t = rename_samples(t)
        samples = t.columns.tolist()
        
    #Create a selection of loci (if specified). Remove loci from the haplotypes table if they are not present in the locus selection.
    loci = list(dict.fromkeys(t.index.get_level_values('Locus').tolist()))
    if args['loci']:
        locusIDs = txt2dict(loci, args['loci'], 'Locus')
        deleted = [x for x in loci if x not in locusIDs]
        if len(deleted) > 0:
            print('\t- {} loc{} removed from the haplotypes table because {} not present in the --loci list:\n\t\t - {}\n'.format(len(deleted), 'i were' if len(deleted) > 1 else 'us was', 'it was' if len(deleted) == 1 else 'they were', '\n\t\t - '.join(deleted)))
            t.drop(deleted, level = 0, inplace = True)
        t = t.reindex(list(locusIDs.keys()), level = 0)
        loci = list(dict.fromkeys(t.index.get_level_values('Locus').tolist()))
        assert len(loci) > 0, 'No loci were retained after the removal of loci. Please disable the --loci option or list minimum one locus ID in the --loci list.'
    
    #Import the list of reference samples. If no list with reference samples is provided, take the first sample in the haplotypes table as a reference.
    if args['reference_samples']:
        ref_samples = txt2dict(samples, args['reference_samples'], 'Ref_sample')
    else:
        ref_samples = {samples[0] : None}
    
    #Get the colour ID of the reference samples.
    ref_samples = get_colour_IDs(ref_samples)
    ref_IDs = [x for x in ref_samples.keys()]
    
    #Count the number of haplotypes in each locus.
    if args['type'] == 'haplotype_count':
        print(' * Counting the number of loci with 1{} haplotype{} ... '.format(' to {}'.format(args['ploidy']) if args['ploidy'] > 1 else '', 's' if args['ploidy'] > 1 else ''))
        indexes = loci
        columns = ['{}{} haplotype{}'.format('More than ' if n > args['ploidy'] else '', min(n, args['ploidy']) if n > 0 else 'No', 's' if n != 1 else '') for n in range(args['ploidy'] + 2)]
    else:
        print(' * Calculating the haplotype similarity between each sample and the reference sample{} ... '.format('s' if len(ref_samples) - 2 > 1 else ''))
        indexes = [t.index.get_level_values('Locus'), t.index.get_level_values('Haplotypes')]
        columns = ref_IDs
    
    with mp.Pool(args['processes']) as p:
        count_dfs = p.map(check_samples, samples)
        
    #Convert the counts data into a dataframe with the samples as indexes and the counts as columns (type == count) or the samples and loci as indexes and the ref_samples as columns.
    count_df = pd.concat(count_dfs)
    
    #Remove rows with only missing data and columns with only missing values or zeros.
    count_df = count_df.loc[:, (count_df.sum(axis=0) > 0)]
    count_df.dropna(how = 'all', inplace = True)
    
    #Remove reference samples from the bar plot.
    if args['type'] == 'reference_similarity' and not args['plot_refs']:
        count_df.drop([x for x in samples if x in ref_IDs], level = 0, inplace = True)
        samples = [x for x in samples if x not in ref_IDs]
    
    #Define plot font settings.
    font_set = set_font_style()
    
    #Print locus information to a tab-delimited text file.
    if args['print_locus_info']:
        count_df.to_csv('{}/{}_{}_loci{}.txt'.format(args['output_directory'], args['type'].capitalize(), args['value_type'] if args['value_type'] != 'count' else 'number', '_{}'.format(args['suffix']) if args['suffix'] else ''), sep = '\t')
    
    #Create heatmap showing the reference similarity for each locus of a reference sample (rows) in each sample (columns).
    if args['plot_heatmap_ref_samples']:
        assert args['unique'], 'Heatmaps are only correctly displayed for datasets with unique haplotypes. Please adjust the options accordingly.'
        print(' * Creating a heatmap of the reference similarity per locus for every sample ... ')
        ref_sample_df = count_df[count_df.columns].apply(pd.to_numeric, errors='coerce').idxmax(1).to_frame(name = 'Ref_sample').reset_index().drop('Sample', axis = 1).drop_duplicates().set_index('Ref_sample').sort_index()
        heatmap_df = pd.DataFrame(index = [ref_sample_df.index.get_level_values('Ref_sample'), ref_sample_df['Locus'], ref_sample_df['Haplotypes']], columns = samples).sort_index()
        for ref_sample in count_df.columns:
            for i, call in count_df.groupby(level = [0, 1, 2]):
                heatmap_df[i[0]].loc[(ref_sample, i[1], i[2])] = count_df[ref_sample].loc[i]
        plot_heatmap(heatmap_df)
    
    #Plot count data.
    print(' * Creating plot ...')
    plot_barplot(count_df)
    
    print(' * Finished!\n')
    print_date()
