# Tutorial SMAPapp-Matrix

## Input data
The example dataset consists of 49 GBS samples (23 species from the genus *Coffea* species and one species from the genus *Tricalysia*) that were retrieved from Bawin *et al*. (2021).
The trimmed sequence data, which are available in the NCBI sequence read archive (BioProject PRJNA612193), and were processed with GBprocesS v3.0.4 (Schaumont, 2020). Reads were merged with PEAR (Zhang *et al*., 2014) if the length of the forward and reverse read was minimum 60 bp
and the overlap between the forward and reverse read was minimum 10 bp. Merged reads with a mean base quality lower than 25 were discarded using the AverageQualityFilter operation in GBprocesS, whereas reads with internal restriction sites were removed with Cutadapt in GBprocesS.
Next, the reads were mapped onto the reference genome sequence of *Coffea canephora* (Denoeud *et al*., 2014) using the BWA-mem algorithm in BWA v0.7.17 (Li & Durbin, 2009) with all parameters left to default. Mapped reads were filtered for a minimum mapping quality score of 20 with Samtools v1.10.0 (Li *et al*., 2009). SNP calling was performed with the Unified Genotyper of the Genome Analysis ToolKit (GATK) v3.7 (Mckenna *et al*., 2010)
and SNPs were filtered for a minimum quality score of 20 and a minimum minor allele count of 2 using VCFtools (Danecek *et al*., 2011). Genotype calls were masked if the genotype quality was lower than 30, if their read depth was lower than 30, or if the allele depth was lower than 3 using VCFtools and a custom python script.
In addition, non-polymorphic and multiallelic SNPs were removed with GATK.

## SMAP

Haplotypes were called using SMAP v4.2.2. GBS stack coordinates were inferred using SMAP delineate with the BAM files filtered for a minimum MAPQ score of 20 as input. The parameters of SMAP delineate were set as follows:

 * -r ignore (ignore mapping orientation)
 * -q 20 (minimum MAPQ of 20) 
 * -x 5 (minimum stack read depth of 5 reads)
 * -y 1500 (maximum stack read depth of 1500)
 * -w 5 (completeness percentage of 5)
 * -c 10 (minimum cluster depth of 10 reads)
 * -s 20 (maximum number of 20 SMAPs)
 * -n Coffea (label of the dataset)

SMAP delineate was run using the following command:

`smap delineate . -r ignore -p 32 -q 20 -x 5 -y 1500 -w 5 -c 10 -s 20 -n Coffea`

The output BED file from SMAP delineate was used as input for SMAP haplotype-sites together with the BAM files filtered for a minimum MAPQ score of 20 and the filtered VCF file. The parameters of SMAP haplotype-sites were set as follows:

 * -r ignore (ignore mapping orientation)
 * -partial include (include haplotypes with partial alignments)
 * --no_indels (haplotypes with indels on SNP positions are excluded)
 * -q 20 (minimum mapping quality of 20)
 * -c 10 (minimum total read depth of 10 reads)
 * -f 5 (minimum haplotype frequency of 5 percent)
 * -e dominant (haplotypes table converted into presence/absence data)
 * --frequency_interval_bounds 10 (A frequency interval bounds of 10 and 90 for transforming haplotype frequencies into discrete calls)
 * --plot all (create all plots)

SMAP haplotype-sites was run using the following command:

`smap haplotype-sites . final_stack_positions_Coffea_C5.0_SMAP20_CL0_inf.bed Coffea_SNPs.vcf -r ignore  -a include -p 32 --no_indels  -q 20 -c 10 -f 5 -e dominant -i 10 --plot all`

The resulting haplotypes table contained 59 179 haplotypes in 6 098 loci that were called in 49 samples.

## Aim
The aim of this tutorial is to determine which species are genetically most similar to the allotetraploid species *Coffea arabica* (Bawin et al., 2021) and to select a small set of loci that are highly discriminatory among the species in the dataset.

## Data exploration (first use SMAPapp-Matrix)
We first used SMAPapp-Matrix to construct line curves of the genetic distance values between a *C. arabica* sample and all other samples in function of a cumulative number of shared loci. This first run can be considered as a dry run to check how many loci are needed to obtain stable genetic distance estimates between *C. arabica* and all other samples. Most parameters were left to default, and only the following parameter settings were used:

- The (Jaccard) genetic similarity values were converted to (inversed) genetic distances (distance = 1 - Jaccard similarity).
- The locus completeness was set to 0.10, so that only loci that are present in at least 10% of the samples (minimum 5 samples) were considered.

Prior to the construction of the line curves, all samples in the haplotypes table were renamed using the --samples option. The --samples list (Sample_names.txt) is a tab-delimited text file containing the original sample names (as listed in the haplotypes table) in the first column and, if preferred, the new sample names in the second column, as shown below (Table 1). We chose to label each sample by its species identity followed by a decimal number. Technical GBS library replicates are indicated by the same unit before the decimal point but a different tenth after the decimal point (e.g. 1.1 and 1.2). Biological replicates of the same species are indicated by the same species name, but a different ID number (e.g. 1.1 and 2.1).

**Table 1.** Example of a list of sample names (left column = original sample name, right column = new sample name) that can be parsed to SMAPapp-Matrix script via the --samples option.
|     |     |
| --- | --- |
| SRR11294243 | Coffea arabica 1.1            |
| SRR11294283 | Coffea arabica 1.2            |
| SRR11294274 | Coffea arabica 2.1            |
| SRR11294273 | Coffea arabica 3.1            |
| SRR12571240 | Tricalysia 1.1                |
| SRR12571239 | Tricalysia 1.2                |
| SRR12571233 | Coffea mannii 1.1             |
| SRR12571254 | Coffea lebruniana 1.1         |
| SRR12571230 | Coffea charrieriana 1.1       |
| SRR12571229 | Coffea charrieriana 1.2       |
| SRR12571258 | Coffea charrieriana 2.1       |
| SRR12571243 | Coffea pseudozanguebariae 1.1 |
| SRR12571244 | Coffea salvatrix 1.1          |
| SRR12571248 | Coffea pocsii 1.1             |
| SRR12571238 | Coffea pocsii 1.2             |
| SRR12571245 | Coffea sessiliflora 1.1       |
| SRR12571237 | Coffea sessiliflora 1.2       |
| SRR12571241 | Coffea racemosa 1.1           |
| SRR12571234 | Coffea racemosa 1.2           |
| SRR12571251 | Coffea heterocalyx 1.1        |
| SRR11294286 | Coffea anthonyi 1.1           |
| SRR11294272 | Coffea anthonyi 1.2           |
| SRR12571252 | Coffea anthonyi 2.1           |
| SRR12571256 | Coffea kivuensis 1.1          |
| SRR11294278 | Coffea kivuensis 2.1          |
| SRR11294249 | Coffea eugenioides 1.1        |
| SRR12571228 | Coffea eugenioides 2.1        |
| SRR12571227 | Coffea eugenioides 2.2        |
| SRR12571255 | Coffea lulandoensis 1.1       |
| SRR12571253 | Coffea lulandoensis 2.1       |
| SRR12571235 | Coffea kapakata 1.1           |
| SRR12571247 | Coffea stenophylla 1.1        |
| SRR12571236 | Coffea stenophylla 2.1        |
| SRR12571246 | Coffea humilis 1.1            |
| SRR11294265 | Coffea liberica 1.1           |
| SRR11294270 | Coffea liberica 1.2           |
| SRR12571231 | Coffea liberica 2.1           |
| SRR11294285 | Coffea brevipes 1.1           |
| SRR11294271 | Coffea brevipes 1.2           |
| SRR12571259 | Coffea brevipes 2.1           |
| SRR11294254 | Coffea congensis 1.1          |
| SRR11294284 | Coffea congensis 1.2          |
| SRR12571242 | Coffea congensis 2.1          |
| SRR12571232 | Coffea canephora 1.1          |
| SRR11294238 | Coffea canephora 2.1          |
| SRR11294277 | Coffea canephora 2.2          |
| SRR12571250 | Coffea macrocarpa 1.1         |
| SRR12571249 | Coffea macrocarpa 1.2         |
| SRR12571257 | Coffea dubardii 1.1           |

The order of the sample names in the --samples list defines the order of the samples in the resulting matrix. For this tutorial, samples were ordered according to their assumed phylogenetic relationships (Hamon *et al*., 2017). The user can choose to reorder but not to rename the samples in the output files by only specifying one column with the original sample names in the --samples list. Sample names that are not listed in the --samples list are removed from the haplotypes table prior to all calculations done by the script. The user can also choose to rename some samples while keeping the original sample name of other samples. If, for example, the user would only like to work with samples from the species *C. arabica*, and the first *C. arabica* sample must retain its original sample name, the following --samples list can be used (Table 2).

**Table 2.** Example of a list of sample names (left column = original sample name, right column = new sample name) containing only samples in the haplotypes table belonging to the species *C. arabica*. The first sample in the list with ID SRR11294243 keeps its original sample ID, while the remaining three samples are renamed so that they have a more convenient sample name in the output files.
|     |     |
| --- | --- |
| SRR11294243 | <NA>               |
| SRR11294283 | Coffea arabica 1.2 |
| SRR11294274 | Coffea arabica 2.1 |
| SRR11294273 | Coffea arabica 3.1 |

Samples belonging to the same group (e.g. species) can be merged into one unit by specifying the same new sample name in the --samples list for all samples that must be merged into the same unit. For instance, all samples in our dataset of 49 samples that belong to the same species could be merged into species units (Table 3).

**Table 3.** Example of a list of sample names (left column = original sample name, right column = new sample name). All samples belonging to the same species are merged into species units.
|     |     |
| --- | --- |
| SRR11294243 | Coffea arabica            |
| SRR11294283 | Coffea arabica            |
| SRR11294274 | Coffea arabica            |
| SRR11294273 | Coffea arabica            |
| SRR12571240 | Tricalysia                |
| SRR12571239 | Tricalysia                |
| SRR12571233 | Coffea mannii             |
| SRR12571254 | Coffea lebruniana         |
| SRR12571230 | Coffea charrieriana       |
| SRR12571229 | Coffea charrieriana       |
| SRR12571258 | Coffea charrieriana       |
| SRR12571243 | Coffea pseudozanguebariae |
| SRR12571244 | Coffea salvatrix          |
| SRR12571248 | Coffea pocsii             |
| SRR12571238 | Coffea pocsii             |
| SRR12571245 | Coffea sessiliflora       |
| SRR12571237 | Coffea sessiliflora       |
| SRR12571241 | Coffea racemosa           |
| SRR12571234 | Coffea racemosa           |
| SRR12571251 | Coffea heterocalyx        |
| SRR11294286 | Coffea anthonyi           |
| SRR11294272 | Coffea anthonyi           |
| SRR12571252 | Coffea anthonyi           |
| SRR12571256 | Coffea kivuensis          |
| SRR11294278 | Coffea kivuensis          |
| SRR11294249 | Coffea eugenioides        |
| SRR12571228 | Coffea eugenioides        |
| SRR12571227 | Coffea eugenioides        |
| SRR12571255 | Coffea lulandoensis       |
| SRR12571253 | Coffea lulandoensis       |
| SRR12571235 | Coffea kapakata           |
| SRR12571247 | Coffea stenophylla        |
| SRR12571236 | Coffea stenophylla        |
| SRR12571246 | Coffea humilis            |
| SRR11294265 | Coffea liberica           |
| SRR11294270 | Coffea liberica           |
| SRR12571231 | Coffea liberica           |
| SRR11294285 | Coffea brevipes           |
| SRR11294271 | Coffea brevipes           |
| SRR12571259 | Coffea brevipes           |
| SRR11294254 | Coffea congensis          |
| SRR11294284 | Coffea congensis          |
| SRR12571242 | Coffea congensis          |
| SRR12571232 | Coffea canephora          |
| SRR11294238 | Coffea canephora          |
| SRR11294277 | Coffea canephora          |
| SRR12571250 | Coffea macrocarpa         |
| SRR12571249 | Coffea macrocarpa         |
| SRR12571257 | Coffea dubardii           |

Line curves are constructed by setting the --plot_line_curves option. By default, a line curve is constructed for every sample pair. The first sample in every sample pair is considered to be the reference sample, and curves of pairs with the same reference sample are plotted in the same graph. However, the user can choose to plot only a subset of sample pairs that is specified in a list. That list (List.txt) must be parsed to the script via the --list_line_curves option as a tab-delimited text file with each row containing a pair of (new) sample names. The sample name in the first column is the reference sample, while the sample in the second column is shown in the line plot.
An example of such a list of sample pairs is included below (Table 4).

**Table 4.** Example of a list containing sample pairs to construct line curves. The sample Coffea arabica 1.1 was used as reference sample.
|     |     |
| --- | --- |
| Coffea arabica 1.1 | Coffea arabica 1.2            |
| Coffea arabica 1.1 | Coffea arabica 2.1            |
| Coffea arabica 1.1 | Coffea arabica 3.1            |
| Coffea arabica 1.1 | Tricalysia 1.1                |
| Coffea arabica 1.1 | Tricalysia 1.2                |
| Coffea arabica 1.1 | Coffea mannii 1.1             |
| Coffea arabica 1.1 | Coffea lebruniana 1.1         |
| Coffea arabica 1.1 | Coffea charrieriana 1.1       |
| Coffea arabica 1.1 | Coffea charrieriana 1.2       |
| Coffea arabica 1.1 | Coffea charrieriana 2.1       |
| Coffea arabica 1.1 | Coffea pseudozanguebariae 1.1 |
| Coffea arabica 1.1 | Coffea salvatrix 1.1          |
| Coffea arabica 1.1 | Coffea pocsii 1.1             |
| Coffea arabica 1.1 | Coffea pocsii 2.1             |
| Coffea arabica 1.1 | Coffea sessiliflora 1.1       |
| Coffea arabica 1.1 | Coffea sessiliflora 2.1       |
| Coffea arabica 1.1 | Coffea racemosa 1.1           |
| Coffea arabica 1.1 | Coffea racemosa 2.1           |
| Coffea arabica 1.1 | Coffea heterocalyx 1.1        |
| Coffea arabica 1.1 | Coffea anthonyi 1.1           |
| Coffea arabica 1.1 | Coffea anthonyi 1.2           |
| Coffea arabica 1.1 | Coffea anthonyi 2.1           |
| Coffea arabica 1.1 | Coffea kivuensis 1.1          |
| Coffea arabica 1.1 | Coffea kivuensis 2.1          |
| Coffea arabica 1.1 | Coffea eugenioides 1.1        |
| Coffea arabica 1.1 | Coffea eugenioides 2.1        |
| Coffea arabica 1.1 | Coffea eugenioides 2.2        |
| Coffea arabica 1.1 | Coffea lulandoensis 1.1       |
| Coffea arabica 1.1 | Coffea lulandoensis 2.1       |
| Coffea arabica 1.1 | Coffea kapakata 1.1           |
| Coffea arabica 1.1 | Coffea stenophylla 1.1        |
| Coffea arabica 1.1 | Coffea stenophylla 2.1        |
| Coffea arabica 1.1 | Coffea humilis 1.1            |
| Coffea arabica 1.1 | Coffea liberica 1.1           |
| Coffea arabica 1.1 | Coffea liberica 1.2           |
| Coffea arabica 1.1 | Coffea liberica 2.1           |
| Coffea arabica 1.1 | Coffea brevipes 1.1           |
| Coffea arabica 1.1 | Coffea brevipes 1.2           |
| Coffea arabica 1.1 | Coffea brevipes 2.1           |
| Coffea arabica 1.1 | Coffea congensis 1.1          |
| Coffea arabica 1.1 | Coffea congensis 1.2          |
| Coffea arabica 1.1 | Coffea congensis 2.1          |
| Coffea arabica 1.1 | Coffea canephora 1.1          |
| Coffea arabica 1.1 | Coffea canephora 2.1          |
| Coffea arabica 1.1 | Coffea canephora 2.2          |
| Coffea arabica 1.1 | Coffea macrocarpa 1.1         |
| Coffea arabica 1.1 | Coffea macrocarpa 1.2         |
| Coffea arabica 1.1 | Coffea dubardii 1.1           |

In the first run of SMAPapp-Matrix, Table 1 and Table 5 were parsed to the script using the --samples option and the --list_line_curves option, respectively. The font size of the legend was set to 3.5 so that the entire legend would fit within the figure frame. The output format of the plots was set to png. The command line parameters are listed below.

`python3 SMAPapp-Matrix.py --table haplotypes_R10_F1_call.txt -s Jaccard --distance --distance_method Inversed -lc 0.10 --samples Sample_names.txt --plot_line_curves --list_line_curves List.txt --plot_format png --legend_fontsize 3.5`

808 loci were discarded due to a locus completeness lower than 0.10, so that all calculations were performed on 5290 loci in 49 samples. The line curves show that genetic distance estimates between *C. arabica* and the other samples seem to stabilize after 500 shared loci, suggesting that 500 shared loci are required to estimate genetic distances between *C. arabica* and the other samples (Fig. 1).

![Line curves showing genetic distance estimates between *C. arabica* and the other samples](/Tutorials/Line_curves.png)
**Figure 1.** Line curves showing genetic distance estimates between *C. arabica* and the other samples in function of the number of shared loci.

## Genetic distance estimation
SMAPapp-Matrix was rerun with a sample completeness of 500. All samples had data for more than 500 loci, so the sample completeness filter did not remove samples from the haplotypes table. The font sizes of the title, labels, and legend were set to 40, 30, and 40, respectively, to increase the readability of the plot. The command line parameters are listed below.

`python3 SMAPapp-Matrix.py --table haplotypes_R10_F1_call.txt -s Jaccard --distance --distance_method Inversed -lc 0.10 --samples Sample_names.txt --plot_format png -sc 500 --print_sample_information All --title_fontsize 40 --label_fontsize 30 --legend_fontsize 40`

The genetic distance matrix shows that all comparisons between *C. arabica* samples and *C. eugenioides* samples on the one hand and between *C. arabica* samples and *C. canephora* samples on the other have a lighter colour than the comparisons between *C. arabica* samples and samples from other species (Fig. 2). From all species in the dataset, *C. eugenioides* and *C. canephora* seem to be genetically most similar to *C. arabica*. The genetic distance between technical replicates was also near zero, showing that the applied procedure is highly reproducible. 

![Original genetic distance matrix showing the genetic distance between each pair of samples in the dataset](/Tutorials/Genetic_distance_matrix_1.png)
**Figure 2.** Genetic distance matrix showing the genetic distance between each pair of samples in the dataset.

The minimum number of shared loci in a sample pair was slightly lower than 800. These numbers were found for all comparisons between a *Tricalysia* sample and a *Coffea* sample, confirming the relatively low genetic similarity between species from both genera compared to the similarity between two *Coffea* species. In general, a higher number of shared loci was found between species belonging to the same cluster in the *Coffea* phylogeny. No clear differences in the number of shared loci were observed between species that were classified into the same cluster.

![The number of shared loci between each pair of species in the dataset.](/Tutorials/NumberOfSharedLoci.png)
**Figure 3.** The number of shared loci between each pair of species in the dataset.

The genetic distance matrix can be further customized. The user can annotate all cells in the matrix with their corresponding value using the --annotate_matrix option. Moreover, the user can mask all cells using the --mask option to exclude cells from the matrix with redundant values. Either all cells above (mask upper) or below (mask lower) the main diagonal can be masked. Fig. 4 is a copy of Fig. 2, but with genetic distance values added to their corresponding cells and with all cells above the main diagonal masked. The font size of the distance values (i.e. the ticks in the graph) was set to 20. The command line parameters are listed below.

`python3 SMAPapp-Matrix.py --table haplotypes_R10_F1_call.txt -s Jaccard --distance --distance_method Inversed -lc 0.10 --samples Sample_names.txt --plot_format png -sc 500 --print_sample_information All --title_fontsize 40 --label_fontsize 30 --legend_fontsize 40 --mask Upper --annotate_matrix --tick_fontsize 20`

The higher genetic similarity between *C. arabica* samples and *C. eugenioides* samples or *C. canephora* samples was confirmed by the annotated values, as *C. eugenioides* and *C. canephora* were the only two species with a genetic distance with *C. arabica* below 0.70 (Fig. 4).

![Genetic distance matrix with annotations of distance values between each pair of samples in the dataset.](Genetic_distance_matrix_2.png)
**Figure 4.** Genetic distance matrix with annotations of distance values between each pair of samples in the dataset. The cells above the main diagonal are masked.

The genetic similarity matrix can be simplified by merging all samples that belong to the same species (Table 3). The used command line options are provided below.

`python3 SMAPapp-Matrix.py --table haplotypes_R10_F1_call.txt -s Jaccard --distance --distance_method Inversed -lc 0.10 --samples Sample_names_merged.txt --plot_format png -sc 500 --print_sample_information All --title_fontsize 40 --label_fontsize 30 --legend_fontsize 40 --mask Upper --annotate_matrix --tick_fontsize 20`

As the genetic similarity between samples from the same species was very high (Fig. 2 and Fig. 4), merging samples into species units resulted in a genetic distance matrix that was similar to the original matrix, but with a lower number of samples (Fig. 5).

![Genetic distance matrix of species units with annotations of distance values between each pair of samples in the dataset.](Genetic_distance_matrix_3.png)
**Figure. 5.** Genetic distance matrix of species units with annotations of distance values between each pair of samples in the dataset. The cells above the main diagonal are masked.

## Inference of locus information

SMAPapp-Matrix also allows the user to gain insight in the information content of all loci in the haplotypes table. The information content of a locus may depend on the research aims. For some applications (e.g. the inference of genetic similarity among samples), a locus is more informative if it contains one or multiple haplotypes that are shared by both samples in a sample pair. For other applications (e.g. sample discrimination), a locus might be more informative if it contains one or multiple haplotypes that are unique to one of the samples in the sample pair. The locus information content can also be expressed for each individual locus as the number of sample pairs wherein a locus was informative compared to the total number of sample pairs wherein the locus was considered. SMAPapp-Matrix was run using haplotypes table with the *Coffea* samples to obtain a matrix with the number of partially shared loci (i.e. the number of loci with at least one shared haplotype) in every sample pair. The command line parameters are provided below.

`python3 SMAPapp-Matrix.py --table haplotypes_R10_F1_call.txt -lc 0.10 --samples Sample_names.txt --plot_format png -sc 500 --print_locus_information All --partial --title_fontsize 40 --label_fontsize 30 --legend_fontsize 40 --mask Upper --annotate_matrix --tick_fontsize 20`

The resulting matrix shows that *C. eugenioides* and *C. canephora* samples had the highest number of loci with at least one shared haplotype with the *C. arabica* samples, which confirms the high genetic similarity between those species in a slightly different way than the Jaccard distance values (Fig. 6). Out of 5290 loci in the dataset, 9 loci had no shared haplotypes in any sample pair. However, these loci were only considered in 6 sample pairs, indicating a high amount of missing data (Table 5). In contrast, 317 loci had at least one shared locus in every sample pair wherein they were considered. The number of sample pairs in which these 317 loci were considered, ranged between 6 and 1 176 (i.e. the maximum amount in a set of 49 samples) (Table 5).

![Locus information matrix showing the number of loci with at least one shared haplotype in every sample pair.](Locus_information_matrix_1.png)
**Figure 6.** Locus information matrix showing the number of loci with at least one shared haplotype in every sample pair.

**Table 5.** Locus information content of all loci with no shared haplotypes in any considered sample pair (first six loci) and with at least one shared haplotype in all considered sample pairs (last 317 loci). The locus ID is shown in the first column followed by the number of sample pairs (NumberOfComparisons, second column), the number of sample pairs with at least one shared haplotype (NumberOfComparisonsWithSharedHaplotypes, third column), and the proportion of the second and third column (ProportionOfComparisonsWithSharedHaplotypes, fourth column).
|         Locus_ID          | NumberOfComparisons | NumberOfComparisonsWithSharedHaplotypes | ProportionOfComparisonsWithSharedHaplotypes |
|         --------          | ------------------- | --------------------------------------- | ------------------------------------------- |
| chr0:62933070-62933161_+  | 6                   | 0                                       | 0                                           |
| chr0:106074087-106074273_+| 6                   | 0                                       | 0                                           |
| chr1:27591904-27592060_+  | 6                   | 0                                       | 0                                           |
| chr4:820489-820641_+      | 6                   | 0                                       | 0                                           |
| chr4:6991265-6991483_+    | 6                   | 0                                       | 0                                           |
| chr5:16945704-16945927_+  | 6                   | 0                                       | 0                                           |
| chr8:5499084-5499242_+    | 6                   | 0                                       | 0                                           |
| chr10:4859293-4859450_+   | 6                   | 0                                       | 0                                           |
| chr11:18380373-18380496_+ | 6                   | 0                                       | 0                                           |
| chr0:2249938-2250100_+    | 15                  | 15                                      | 1                                           |
| chr0:8863005-8863146_+    | 10                  | 10                                      | 1                                           |
| chr0:10594156-10594366_+  | 6                   | 6                                       | 1                                           |
| chr0:10595501-10595707_+  | 10                  | 10                                      | 1                                           |
| chr0:10617024-10617071_+  | 136                 | 136                                     | 1                                           |
| chr0:10760202-10760352_+  | 15                  | 15                                      | 1                                           |
| chr0:10828291-10828428_+  | 6                   | 6                                       | 1                                           |
| chr0:10877734-10877954_+  | 253                 | 253                                     | 1                                           |
| chr0:10930739-10930928_+  | 36                  | 36                                      | 1                                           |
| chr0:10965804-10965970_+  | 10                  | 10                                      | 1                                           |
| chr0:10996242-10996359_+  | 496                 | 496                                     | 1                                           |
| chr0:11086580-11086704_+  | 105                 | 105                                     | 1                                           |
| chr0:11175626-11175729_+  | 6                   | 6                                       | 1                                           |
| chr0:11223148-11223246_+  | 28                  | 28                                      | 1                                           |
| chr0:11231054-11231213_+  | 6                   | 6                                       | 1                                           |
| chr0:11267377-11267442_+  | 36                  | 36                                      | 1                                           |
| chr0:11289946-11290108_+  | 55                  | 55                                      | 1                                           |
| chr0:11293319-11293394_+  | 6                   | 6                                       | 1                                           |
| chr0:11368007-11368118_+  | 6                   | 6                                       | 1                                           |
| chr0:14754719-14754895_+  | 946                 | 946                                     | 1                                           |
| chr0:32273340-32273400_+  | 10                  | 10                                      | 1                                           |
| chr0:33664152-33664347_+  | 6                   | 6                                       | 1                                           |
| chr0:40005040-40005224_+  | 28                  | 28                                      | 1                                           |
| chr0:40637051-40637124_+  | 28                  | 28                                      | 1                                           |
| chr0:42246793-42246855_+  | 120                 | 120                                     | 1                                           |
| chr0:52315303-52315375_+  | 21                  | 21                                      | 1                                           |
| chr0:52359906-52360079_+  | 10                  | 10                                      | 1                                           |
| chr0:61061596-61061843_+  | 6                   | 6                                       | 1                                           |
| chr0:83941767-83941954_+  | 6                   | 6                                       | 1                                           |
| chr0:84547247-84547551_+  | 6                   | 6                                       | 1                                           |
| chr0:84951706-84951822_+  | 10                  | 10                                      | 1                                           |
| chr0:85882908-85883047_+  | 6                   | 6                                       | 1                                           |
| chr0:90149297-90149465_+  | 10                  | 10                                      | 1                                           |
| chr0:92032747-92032894_+  | 6                   | 6                                       | 1                                           |
| chr0:95586151-95586225_+  | 6                   | 6                                       | 1                                           |
| chr0:105912989-105913092_+| 21                  | 21                                      | 1                                           |
| chr0:116333871-116334077_+| 10                  | 10                                      | 1                                           |
| chr0:117990597-117990804_+| 6                   | 6                                       | 1                                           |
| chr0:119293826-119294057_+| 6                   | 6                                       | 1                                           |
| chr0:124409536-124409672_+| 78                  | 78                                      | 1                                           |
| chr0:128754621-128754809_+| 21                  | 21                                      | 1                                           |
| chr0:132296554-132296768_+| 6                   | 6                                       | 1                                           |
| chr0:133363707-133363971_+| 21                  | 21                                      | 1                                           |
| chr0:144852314-144852429_+| 6                   | 6                                       | 1                                           |
| chr0:148124186-148124330_+| 21                  | 21                                      | 1                                           |
| chr0:158450585-158450698_+| 66                  | 66                                      | 1                                           |
| chr0:173198954-173199216_+| 36                  | 36                                      | 1                                           |
| chr0:183208274-183208442_+| 6                   | 6                                       | 1                                           |
| chr0:183891233-183891416_+| 6                   | 6                                       | 1                                           |
| chr0:191044820-191044948_+| 6                   | 6                                       | 1                                           |
| chr0:192817044-192817128_+| 10                  | 10                                      | 1                                           |
| chr0:194916476-194916566_+| 120                 | 120                                     | 1                                           |
| chr0:195852127-195852319_+| 6                   | 6                                       | 1                                           |
| chr1:2095153-2095386_+    | 55                  | 55                                      | 1                                           |
| chr1:2098543-2098687_+    | 6                   | 6                                       | 1                                           |
| chr1:2227975-2228025_+    | 6                   | 6                                       | 1                                           |
| chr1:9507622-9507690_+    | 6                   | 6                                       | 1                                           |
| chr1:24041005-24041141_+  | 6                   | 6                                       | 1                                           |
| chr1:24382850-24382952_+  | 6                   | 6                                       | 1                                           |
| chr1:25838353-25838413_+  | 21                  | 21                                      | 1                                           |
| chr1:26052730-26052858_+  | 10                  | 10                                      | 1                                           |
| chr1:27179485-27179633_+  | 6                   | 6                                       | 1                                           |
| chr1:27452893-27453056_+  | 6                   | 6                                       | 1                                           |
| chr1:28281032-28281190_+  | 6                   | 6                                       | 1                                           |
| chr1:28760197-28760479_+  | 6                   | 6                                       | 1                                           |
| chr1:29863576-29863725_+  | 6                   | 6                                       | 1                                           |
| chr1:31244728-31244798_+  | 6                   | 6                                       | 1                                           |
| chr1:31290674-31290811_+  | 6                   | 6                                       | 1                                           |
| chr1:31730556-31730782_+  | 6                   | 6                                       | 1                                           |
| chr1:31807792-31807977_+  | 6                   | 6                                       | 1                                           |
| chr1:32884920-32885053_+  | 15                  | 15                                      | 1                                           |
| chr1:33362824-33362865_+  | 91                  | 91                                      | 1                                           |
| chr1:35675542-35675688_+  | 10                  | 10                                      | 1                                           |
| chr1:35714114-35714315_+  | 6                   | 6                                       | 1                                           |
| chr1:37227877-37227955_+  | 15                  | 15                                      | 1                                           |
| chr1:37368726-37368879_+  | 10                  | 10                                      | 1                                           |
| chr1:37891791-37891906_+  | 6                   | 6                                       | 1                                           |
| chr2:1668971-1669045_+    | 1176                | 1176                                    | 1                                           |
| chr2:2695992-2696126_+    | 6                   | 6                                       | 1                                           |
| chr2:2793148-2793348_+    | 45                  | 45                                      | 1                                           |
| chr2:3360987-3361232_+    | 21                  | 21                                      | 1                                           |
| chr2:3488841-3488934_+    | 6                   | 6                                       | 1                                           |
| chr2:3844586-3844628_+    | 78                  | 78                                      | 1                                           |
| chr2:4313522-4313741_+    | 6                   | 6                                       | 1                                           |
| chr2:4770474-4770656_+    | 6                   | 6                                       | 1                                           |
| chr2:6006953-6007169_+    | 6                   | 6                                       | 1                                           |
| chr2:7215938-7216083_+    | 171                 | 171                                     | 1                                           |
| chr2:7352296-7352489_+    | 6                   | 6                                       | 1                                           |
| chr2:7352496-7352736_+    | 6                   | 6                                       | 1                                           |
| chr2:9450052-9450161_+    | 10                  | 10                                      | 1                                           |
| chr2:9753164-9753216_+    | 6                   | 6                                       | 1                                           |
| chr2:9807310-9807387_+    | 120                 | 120                                     | 1                                           |
| chr2:11902612-11902769_+  | 78                  | 78                                      | 1                                           |
| chr2:12864038-12864147_+  | 45                  | 45                                      | 1                                           |
| chr2:12872756-12872816_+  | 6                   | 6                                       | 1                                           |
| chr2:12872823-12872874_+  | 6                   | 6                                       | 1                                           |
| chr2:14119410-14119609_+  | 6                   | 6                                       | 1                                           |
| chr2:15335469-15335621_+  | 21                  | 21                                      | 1                                           |
| chr2:15700634-15700849_+  | 21                  | 21                                      | 1                                           |
| chr2:17052827-17052921_+  | 45                  | 45                                      | 1                                           |
| chr2:17785751-17786057_+  | 6                   | 6                                       | 1                                           |
| chr2:18404954-18405090_+  | 6                   | 6                                       | 1                                           |
| chr2:18658549-18658700_+  | 6                   | 6                                       | 1                                           |
| chr2:18742720-18742857_+  | 6                   | 6                                       | 1                                           |
| chr2:19173664-19173765_+  | 15                  | 15                                      | 1                                           |
| chr2:19662548-19662713_+  | 6                   | 6                                       | 1                                           |
| chr2:22404894-22405145_+  | 6                   | 6                                       | 1                                           |
| chr2:22448076-22448214_+  | 10                  | 10                                      | 1                                           |
| chr2:24359515-24359571_+  | 6                   | 6                                       | 1                                           |
| chr2:27011118-27011193_+  | 21                  | 21                                      | 1                                           |
| chr2:31017229-31017485_+  | 15                  | 15                                      | 1                                           |
| chr2:38086109-38086315_+  | 15                  | 15                                      | 1                                           |
| chr2:39453026-39453186_+  | 10                  | 10                                      | 1                                           |
| chr2:39533686-39533722_+  | 21                  | 21                                      | 1                                           |
| chr2:41196401-41196544_+  | 6                   | 6                                       | 1                                           |
| chr2:41871940-41872130_+  | 6                   | 6                                       | 1                                           |
| chr2:43433758-43433835_+  | 153                 | 153                                     | 1                                           |
| chr2:44565835-44566026_+  | 66                  | 66                                      | 1                                           |
| chr2:45416746-45416913_+  | 21                  | 21                                      | 1                                           |
| chr2:51282658-51282894_+  | 6                   | 6                                       | 1                                           |
| chr2:51395648-51395889_+  | 6                   | 6                                       | 1                                           |
| chr2:51425407-51425482_+  | 10                  | 10                                      | 1                                           |
| chr2:52416296-52416522_+  | 120                 | 120                                     | 1                                           |
| chr2:52547347-52547565_+  | 10                  | 10                                      | 1                                           |
| chr2:53053249-53053413_+  | 6                   | 6                                       | 1                                           |
| chr2:53308470-53308686_+  | 6                   | 6                                       | 1                                           |
| chr2:53393529-53393687_+  | 21                  | 21                                      | 1                                           |
| chr2:53455413-53455505_+  | 6                   | 6                                       | 1                                           |
| chr2:54471422-54471569_+  | 6                   | 6                                       | 1                                           |
| chr3:674480-674634_+      | 6                   | 6                                       | 1                                           |
| chr3:1134857-1134912_+    | 28                  | 28                                      | 1                                           |
| chr3:3178976-3179047_+    | 6                   | 6                                       | 1                                           |
| chr3:4723550-4723687_+    | 6                   | 6                                       | 1                                           |
| chr3:4772070-4772145_+    | 6                   | 6                                       | 1                                           |
| chr3:4782080-4782149_+    | 10                  | 10                                      | 1                                           |
| chr3:5295788-5295998_+    | 6                   | 6                                       | 1                                           |
| chr3:5380731-5380937_+    | 6                   | 6                                       | 1                                           |
| chr3:8410901-8411008_+    | 10                  | 10                                      | 1                                           |
| chr3:12054581-12054726_+  | 36                  | 36                                      | 1                                           |
| chr3:21859176-21859238_+  | 6                   | 6                                       | 1                                           |
| chr3:26235695-26235754_+  | 15                  | 15                                      | 1                                           |
| chr3:26608090-26608284_+  | 6                   | 6                                       | 1                                           |
| chr3:27900756-27900893_+  | 6                   | 6                                       | 1                                           |
| chr3:29614774-29614814_+  | 28                  | 28                                      | 1                                           |
| chr3:29806239-29806301_+  | 36                  | 36                                      | 1                                           |
| chr3:30460871-30461071_+  | 10                  | 10                                      | 1                                           |
| chr3:30855317-30855586_+  | 6                   | 6                                       | 1                                           |
| chr3:31700310-31700555_+  | 21                  | 21                                      | 1                                           |
| chr4:2192806-2193068_+    | 6                   | 6                                       | 1                                           |
| chr4:3479854-3480028_+    | 21                  | 21                                      | 1                                           |
| chr4:3647420-3647576_+    | 6                   | 6                                       | 1                                           |
| chr4:4491842-4491969_+    | 6                   | 6                                       | 1                                           |
| chr4:5033469-5033508_+    | 15                  | 15                                      | 1                                           |
| chr4:7033861-7033909_+    | 6                   | 6                                       | 1                                           |
| chr4:7169912-7170167_+    | 6                   | 6                                       | 1                                           |
| chr4:7774958-7775086_+    | 66                  | 66                                      | 1                                           |
| chr4:9754049-9754177_+    | 6                   | 6                                       | 1                                           |
| chr4:9754183-9754284_+    | 6                   | 6                                       | 1                                           |
| chr4:10304491-10304705_+  | 6                   | 6                                       | 1                                           |
| chr4:13445127-13445306_+  | 6                   | 6                                       | 1                                           |
| chr4:13606846-13606953_+  | 6                   | 6                                       | 1                                           |
| chr4:13972388-13972553_+  | 10                  | 10                                      | 1                                           |
| chr4:15190834-15191073_+  | 6                   | 6                                       | 1                                           |
| chr4:15357219-15357415_+  | 15                  | 15                                      | 1                                           |
| chr4:19810967-19811135_+  | 6                   | 6                                       | 1                                           |
| chr4:23709668-23709760_+  | 10                  | 10                                      | 1                                           |
| chr4:24451313-24451510_+  | 6                   | 6                                       | 1                                           |
| chr4:25352355-25352471_+  | 6                   | 6                                       | 1                                           |
| chr4:25846809-25847075_+  | 10                  | 10                                      | 1                                           |
| chr5:1378597-1378756_+    | 6                   | 6                                       | 1                                           |
| chr5:1990067-1990126_+    | 15                  | 15                                      | 1                                           |
| chr5:3204906-3205042_+    | 10                  | 10                                      | 1                                           |
| chr5:5387187-5387434_+    | 6                   | 6                                       | 1                                           |
| chr5:8276008-8276230_+    | 10                  | 10                                      | 1                                           |
| chr5:16240885-16241108_+  | 6                   | 6                                       | 1                                           |
| chr5:17054461-17054569_+  | 6                   | 6                                       | 1                                           |
| chr5:19607404-19607455_+  | 21                  | 21                                      | 1                                           |
| chr5:20181951-20182175_+  | 28                  | 28                                      | 1                                           |
| chr5:23632183-23632245_+  | 15                  | 15                                      | 1                                           |
| chr5:24752734-24752834_+  | 10                  | 10                                      | 1                                           |
| chr5:26458155-26458409_+  | 21                  | 21                                      | 1                                           |
| chr5:26469248-26469431_+  | 6                   | 6                                       | 1                                           |
| chr5:27810399-27810439_+  | 6                   | 6                                       | 1                                           |
| chr6:222444-222657_+      | 10                  | 10                                      | 1                                           |
| chr6:251066-251332_+      | 10                  | 10                                      | 1                                           |
| chr6:2619373-2619529_+    | 6                   | 6                                       | 1                                           |
| chr6:2908507-2908597_+    | 120                 | 120                                     | 1                                           |
| chr6:3599559-3599799_+    | 6                   | 6                                       | 1                                           |
| chr6:3960527-3960696_+    | 10                  | 10                                      | 1                                           |
| chr6:4280746-4280798_+    | 10                  | 10                                      | 1                                           |
| chr6:5162767-5162992_+    | 15                  | 15                                      | 1                                           |
| chr6:6445268-6445413_+    | 28                  | 28                                      | 1                                           |
| chr6:6764186-6764253_+    | 15                  | 15                                      | 1                                           |
| chr6:6964479-6964665_+    | 6                   | 6                                       | 1                                           |
| chr6:7149663-7149776_+    | 21                  | 21                                      | 1                                           |
| chr6:7296746-7296829_+    | 15                  | 15                                      | 1                                           |
| chr6:7619793-7620028_+    | 6                   | 6                                       | 1                                           |
| chr6:7874158-7874404_+    | 6                   | 6                                       | 1                                           |
| chr6:8910652-8910804_+    | 21                  | 21                                      | 1                                           |
| chr6:9300281-9300348_+    | 105                 | 105                                     | 1                                           |
| chr6:9802617-9802738_+    | 15                  | 15                                      | 1                                           |
| chr6:12592217-12592306_+  | 78                  | 78                                      | 1                                           |
| chr6:12896381-12896497_+  | 15                  | 15                                      | 1                                           |
| chr6:13933817-13933995_+  | 15                  | 15                                      | 1                                           |
| chr6:16487425-16487489_+  | 45                  | 45                                      | 1                                           |
| chr6:16835171-16835311_+  | 15                  | 15                                      | 1                                           |
| chr6:17681034-17681167_+  | 21                  | 21                                      | 1                                           |
| chr6:20482801-20483016_+  | 15                  | 15                                      | 1                                           |
| chr6:20919123-20919271_+  | 78                  | 78                                      | 1                                           |
| chr6:30723602-30723697_+  | 6                   | 6                                       | 1                                           |
| chr6:34406979-34407033_+  | 6                   | 6                                       | 1                                           |
| chr6:35783226-35783290_+  | 15                  | 15                                      | 1                                           |
| chr6:35862271-35862540_+  | 6                   | 6                                       | 1                                           |
| chr6:35915611-35915648_+  | 6                   | 6                                       | 1                                           |
| chr7:2286019-2286300_+    | 6                   | 6                                       | 1                                           |
| chr7:3680772-3680900_+    | 36                  | 36                                      | 1                                           |
| chr7:3999253-3999612_+    | 10                  | 10                                      | 1                                           |
| chr7:4036895-4036983_+    | 6                   | 6                                       | 1                                           |
| chr7:4090051-4090226_+    | 6                   | 6                                       | 1                                           |
| chr7:4759288-4759419_+    | 6                   | 6                                       | 1                                           |
| chr7:5283769-5283843_+    | 10                  | 10                                      | 1                                           |
| chr7:5426461-5426650_+    | 21                  | 21                                      | 1                                           |
| chr7:8152463-8152567_+    | 21                  | 21                                      | 1                                           |
| chr7:8382698-8382945_+    | 6                   | 6                                       | 1                                           |
| chr7:8621721-8621774_+    | 28                  | 28                                      | 1                                           |
| chr7:8752883-8753092_+    | 15                  | 15                                      | 1                                           |
| chr7:8787200-8787414_+    | 21                  | 21                                      | 1                                           |
| chr7:8964027-8964147_+    | 15                  | 15                                      | 1                                           |
| chr7:9358565-9358631_+    | 6                   | 6                                       | 1                                           |
| chr7:9371326-9371485_+    | 6                   | 6                                       | 1                                           |
| chr7:10159908-10160066_+  | 6                   | 6                                       | 1                                           |
| chr7:10747440-10747593_+  | 45                  | 45                                      | 1                                           |
| chr7:11797266-11797427_+  | 6                   | 6                                       | 1                                           |
| chr7:12209545-12209712_+  | 78                  | 78                                      | 1                                           |
| chr7:13026711-13026953_+  | 6                   | 6                                       | 1                                           |
| chr7:15576933-15577072_+  | 6                   | 6                                       | 1                                           |
| chr7:16604321-16604572_+  | 6                   | 6                                       | 1                                           |
| chr7:18994551-18994612_+  | 6                   | 6                                       | 1                                           |
| chr7:19008532-19008631_+  | 6                   | 6                                       | 1                                           |
| chr8:164037-164240_+      | 6                   | 6                                       | 1                                           |
| chr8:1661956-1662023_+    | 78                  | 78                                      | 1                                           |
| chr8:3603195-3603273_+    | 10                  | 10                                      | 1                                           |
| chr8:3990772-3990856_+    | 21                  | 21                                      | 1                                           |
| chr8:6067283-6067398_+    | 10                  | 10                                      | 1                                           |
| chr8:8441637-8441710_+    | 28                  | 28                                      | 1                                           |
| chr8:10610807-10610966_+  | 45                  | 45                                      | 1                                           |
| chr8:11033403-11033482_+  | 6                   | 6                                       | 1                                           |
| chr8:16590528-16590698_+  | 28                  | 28                                      | 1                                           |
| chr8:17400875-17401071_+  | 6                   | 6                                       | 1                                           |
| chr8:18873450-18873626_+  | 15                  | 15                                      | 1                                           |
| chr8:20818466-20818511_+  | 6                   | 6                                       | 1                                           |
| chr8:21526977-21527126_+  | 21                  | 21                                      | 1                                           |
| chr8:22088319-22088535_+  | 21                  | 21                                      | 1                                           |
| chr8:22189345-22189602_+  | 15                  | 15                                      | 1                                           |
| chr8:22505730-22505794_+  | 28                  | 28                                      | 1                                           |
| chr8:22815344-22815514_+  | 171                 | 171                                     | 1                                           |
| chr8:24663464-24663598_+  | 6                   | 6                                       | 1                                           |
| chr8:25028058-25028204_+  | 6                   | 6                                       | 1                                           |
| chr8:25989531-25989612_+  | 6                   | 6                                       | 1                                           |
| chr8:25989619-25989685_+  | 6                   | 6                                       | 1                                           |
| chr8:27944102-27944212_+  | 10                  | 10                                      | 1                                           |
| chr8:28117888-28117949_+  | 15                  | 15                                      | 1                                           |
| chr8:29022591-29022696_+  | 21                  | 21                                      | 1                                           |
| chr8:30029825-30029944_+  | 465                 | 465                                     | 1                                           |
| chr8:31418852-31418955_+  | 21                  | 21                                      | 1                                           |
| chr9:526599-526669_+      | 15                  | 15                                      | 1                                           |
| chr9:776759-776888_+      | 6                   | 6                                       | 1                                           |
| chr9:892723-892989_+      | 78                  | 78                                      | 1                                           |
| chr9:1021185-1021365_+    | 6                   | 6                                       | 1                                           |
| chr9:1510005-1510270_+    | 15                  | 15                                      | 1                                           |
| chr9:2321008-2321167_+    | 6                   | 6                                       | 1                                           |
| chr9:4199345-4199412_+    | 15                  | 15                                      | 1                                           |
| chr9:7018869-7018994_+    | 6                   | 6                                       | 1                                           |
| chr9:7225731-7225911_+    | 6                   | 6                                       | 1                                           |
| chr9:7349607-7349745_+    | 15                  | 15                                      | 1                                           |
| chr9:7505220-7505323_+    | 6                   | 6                                       | 1                                           |
| chr9:15127101-15127209_+  | 10                  | 10                                      | 1                                           |
| chr9:15606721-15607058_+  | 10                  | 10                                      | 1                                           |
| chr9:16425516-16425586_+  | 6                   | 6                                       | 1                                           |
| chr9:20658246-20658393_+  | 6                   | 6                                       | 1                                           |
| chr9:21288360-21288432_+  | 15                  | 15                                      | 1                                           |
| chr10:341770-341974_+     | 6                   | 6                                       | 1                                           |
| chr10:700036-700218_+     | 6                   | 6                                       | 1                                           |
| chr10:3640170-3640328_+   | 6                   | 6                                       | 1                                           |
| chr10:7848538-7848700_+   | 6                   | 6                                       | 1                                           |
| chr10:8317778-8317822_+   | 6                   | 6                                       | 1                                           |
| chr10:10558783-10558867_+ | 6                   | 6                                       | 1                                           |
| chr10:13330672-13330766_+ | 15                  | 15                                      | 1                                           |
| chr10:20319234-20319354_+ | 15                  | 15                                      | 1                                           |
| chr10:20805775-20805842_+ | 55                  | 55                                      | 1                                           |
| chr10:22250135-22250305_+ | 6                   | 6                                       | 1                                           |
| chr10:22494859-22495118_+ | 153                 | 153                                     | 1                                           |
| chr10:23046997-23047157_+ | 78                  | 78                                      | 1                                           |
| chr10:23335769-23335822_+ | 15                  | 15                                      | 1                                           |
| chr10:25991804-25991895_+ | 6                   | 6                                       | 1                                           |
| chr10:26441198-26441310_+ | 6                   | 6                                       | 1                                           |
| chr10:26905981-26906183_+ | 6                   | 6                                       | 1                                           |
| chr11:6459-6655_+         | 6                   | 6                                       | 1                                           |
| chr11:653341-653449_+     | 6                   | 6                                       | 1                                           |
| chr11:4601688-4601773_+   | 15                  | 15                                      | 1                                           |
| chr11:7694004-7694126_+   | 136                 | 136                                     | 1                                           |
| chr11:11705636-11705833_+ | 6                   | 6                                       | 1                                           |
| chr11:12437123-12437368_+ | 6                   | 6                                       | 1                                           |
| chr11:18304572-18304644_+ | 136                 | 136                                     | 1                                           |
| chr11:21087860-21088085_+ | 15                  | 15                                      | 1                                           |
| chr11:23797255-23797318_+ | 15                  | 15                                      | 1                                           |
| chr11:24570515-24570652_+ | 36                  | 36                                      | 1                                           |
| chr11:25355024-25355235_+ | 10                  | 10                                      | 1                                           |
| chr11:25385565-25385773_+ | 6                   | 6                                       | 1                                           |
| chr11:25733100-25733185_+ | 6                   | 6                                       | 1                                           |
| chr11:27455428-27455668_+ | 6                   | 6                                       | 1                                           |
| chr11:30041620-30041821_+ | 6                   | 6                                       | 1                                           |
| chr11:31663047-31663245_+ | 21                  | 21                                      | 1                                           |
| chr11:31799568-31799762_+ | 6                   | 6                                       | 1                                           |
| chr11:32219676-32219872_+ | 15                  | 15                                      | 1                                           |
| chr11:32367940-32368000_+ | 561                 | 561                                     | 1                                           |
| chr11:33507671-33507791_+ | 91                  | 91                                      | 1                                           |

To select a list of loci with for every species at least one locus with only unique haplotypes, a list with the locus information content was created using the following command line options.

`python3 SMAPapp-Matrix.py --table haplotypes_R10_F1_call.txt -lc 0.10 --samples Sample_names_merged.txt --plot_format png -sc 500 --print_locus_information List -lic Unique --title_fontsize 40 --label_fontsize 30 --legend_fontsize 40 --mask Upper --annotate_matrix --tick_fontsize 20`

Six loci were selected because of their high discriminatory power among the *Coffea* and *Tricalysia* species in the dataset. All loci were considered in every comparison between species units (276 in total). Every locus had species-specific haplotypes for at least nine species; with locus chr1:30787505-30787768_+ having species-specific haplotypes for 16 out of 24 species.

**Table 6.** List of loci with only unique haplotypes in at least 9 out of 24 species in the dataset. The locus ID is shown in the first column followed by the number of sample pairs (NumberOfComparisons, second column), the number of sample pairs with only unique haplotypes (NumberOfComparisonsWithUniqueHaplotypes, third column), the proportion of the second and third column (ProportionOfComparisonsWithUniqueHaplotypes, fourth column), and the sample names with only unique haplotypes in every comparison (SamplesWithUniqueHaplotypes, fifth column).
|         Locus_ID         | NumberOfComparisons | NumberOfComparisonsWithUniqueHaplotypes | ProportionOfComparisonsWithUniqueHaplotypes |                                                                                                                                             SamplesWithUniqueHaplotypes                                                                                                                                              |
|         --------         | ------------------- | ----------------------------------      | --------------------------------------      |                                                                                                                                             ---------------------------                                                                                                                                              |
| chr1:30787505-30787768_+ | 276                 | 270                                     | 0.98                                        | Coffea anthonyi, Coffea brevipes, Coffea charrieriana, Coffea dubardii, Coffea kapakata, Coffea lebruniana, Coffea liberica, Coffea lulandoensis, Coffea macrocarpa, Coffea mannii, Coffea pocsii, Coffea pseudozanguebariae, Coffea racemosa, Coffea salvatrix, Coffea sessiliflora, Coffea stenophylla, Tricalysia |
| chr2:1006717-1006879_+   | 276                 | 218                                     | 0.79                                        | Coffea brevipes, Coffea charrieriana, Coffea congensis, Coffea eugenioides, Coffea kapakata, Coffea liberica, Coffea mannii, Coffea stenophylla, Tricalysia                                                                                                                                                          |
| chr3:2510192-2510423_+   | 276                 | 261                                     | 0.95                                        | Coffea charrieriana, Coffea congensis, Coffea humilis, Coffea kapakata, Coffea kivuensis, Coffea lebruniana, Coffea liberica, Coffea macrocarpa, Coffea mannii, Coffea racemosa, Coffea salvatrix, Coffea stenophylla, Tricalysia                                                                                    |
| chr3:3816695-3816946_+   | 276                 | 267                                     | 0.97                                        | Coffea canephora, Coffea charrieriana, Coffea dubardii, Coffea humilis, Coffea kapakata, Coffea lebruniana, Coffea liberica, Coffea macrocarpa, Coffea mannii, Coffea pseudozanguebariae, Coffea racemosa, Coffea salvatrix, Coffea stenophylla, Tricalysia                                                          |
| chr7:20976690-20976906_+ | 276                 | 268                                     | 0.97                                        | Coffea anthonyi, Coffea brevipes, Coffea charrieriana, Coffea dubardii, Coffea heterocalyx, Coffea kapakata, Coffea lebruniana, Coffea lulandoensis, Coffea macrocarpa, Coffea mannii, Coffea pseudozanguebariae, Coffea racemosa, Coffea salvatrix, Coffea stenophylla, Tricalysia                                  |
| chr10:2239337-2239448_+  | 276                 | 220                                     | 0.8                                         | Coffea anthonyi, Coffea arabica, Coffea dubardii, Coffea humilis, Coffea kapakata, Coffea lebruniana, Coffea liberica, Coffea macrocarpa, Coffea mannii, Coffea pseudozanguebariae, Coffea stenophylla, Tricalysia                                                                                                   |

A Jaccard genetic distance matrix and a locus information matrix was created based on the haplotype data of the six selected loci using the following command line options.

`python3 SMAPapp-Matrix.py -s Jaccard --table haplotypes_R10_F1_call.txt --distance --distance_method Inversed -lc 0.10 --samples Sample_names_merged.txt --loci Loci.txt --plot_format png --print_sample_information Plot --print_locus_information Plot -lic Unique --title_fontsize 40 --label_fontsize 30 --legend_fontsize 40 --mask Upper --annotate_matrix --tick_fontsize 20`

The Jaccard genetic similarity matrix (Fig. 7) shows that the majority of the species units had a very high genetic distance value. The lowest value (0.67) was found between *C. sessiliflora* and *C. pocsii*. Hence, the haplotype diversity in the set of six loci is highly suitable to discriminate the 24 species in the dataset from each other.

![Jaccard genetic distance matrix of species units based on the selection of six loci.](Genetic_distance_matrix_4.png)
**Figure 7.** Jaccard genetic distance matrix of species units based on the selection of six loci.

All pairwise comparisons between species units had at least one locus with only unique haplotypes between the two samples (Fig. 8). Many comparisons had only unique haplotypes in all six loci. Consequently, SMAPapp-Matrix provides a fast and accurate strategy to select loci in a large haplotypes table with a high discriminatory power among the samples in the Table. Such a set of loci could be targeted in a larger sampling using Sanger sequencing or targeted high-throughput sequencing methods (e.g. multiplex amplicon sequencing).

![Locus information matrix of species units based on the selection of six loci.](Locus_information_matrix_2.png)
**Figure 8.** Locus information matrix of species units based on the selection of six loci.

## Phylogenetic tree reconstruction

The pairwise genetic distance matrix with the species units was recreated in Phylip format together with 100 bootstrap replicates using the following command line options.

`python3 SMAPapp-Matrix.py --table haplotypes_R10_F1_call.txt -s Jaccard --distance --distance_method Inversed -lc 0.10 --samples Sample_names_merged.txt -sc 500 --print_sample_information Matrix -b 100`

The resulting distance matrix was used to infer a neighbor joining phylogenetic tree using the *nj* function from the R package ape (Paradis & Schliep, 2018) in R v4.1.2 (R Core Team, 2021).
Bootstrap values were calculated using the 100 bootstrap replicates and plotted on the nodes of the neighbor joining phylogenetic tree that was created with the original distance matrix.
The resulting phylogenetic tree shows well-supported phylogenetic relationships between the included *Coffea* species (Fig. 9). Due to its allotetraploid origin, the phylogenetic relationships between *C. arabica* and the other *Coffea* species could not be correctly displayed in the phylogenetic tree. The lower bootstrap value assigned to the node that separated the *C. canephora* clade and *C. eugenioides* clade is likely caused by the presence of the *C. arabica* sample in the *C. canephora* clade.

![Neighbor joining phylogenetic tree of *Coffea* with *Tricalysia* sp. as outgroup species. Node labels indicate bootstrap values.](/Tutorials/NJ_phylogenetic_tree.png)
**Figure 9.** Neighbor joining phylogenetic tree of *Coffea* with *Tricalysia* sp. as outgroup species. Node labels indicate bootstrap values.

## Phylogenetic network reconstruction
To demonstrate the hybrid origin of *Coffea arabica*, a phylogenetic network was reconstructed based on the Jaccard genetic distance matrix in Nexus format.
This distance matrix was first created in Nexus format using the following command line parameters. Whitespaces in sample names were replaced by underscores.

`python3 SMAPapp-Matrix.py --table haplotypes_R10_F1_call.txt -s Jaccard --distance --distance_method Inversed -lc 0.10 --samples Sample_names_merged_without_whitespace.txt -sc 500 --print_sample_information Matrix --matrix_format Nexus`

Afterwards, the matrix was imported in SplitsTree v4.14.8 (Huson & Bryant, 2006) and a phylogenetic network was reconstructed using the SplitDecomposition method.
The resulting phylogenetic network shows that *C. arabica* is a hybrid species that is genetically most similar to *C. eugenioides* and *C. canephora* (Fig 10).

![Phylogenetic network of *Coffea* that was reconstructed using the SplitDecomposition method in SplitsTree.](/Tutorials/Phylogenetic_network.png)
**Figure 10.** Phylogenetic network of *Coffea* that was reconstructed using the SplitDecomposition method in SplitsTree.

## References
Bawin, Y., Ruttink, T., Staelens, A., Haegeman, A., Stoffelen, P., Mwanga Mwanga, J.-C.I., … Janssens, S.B. (2021). Phylogenomic analysis clarifies the evolutionary origin of *Coffea arabica*. *Journal of Systematics and Evolution, 59*(5), 953-963. https://doi.org/10.1111/jse.12694

Danecek, P., Auton, A., Abecasis, G., Albers, C. A., Banks, E., DePristo, M. A., … Durbin, R. (2011). The variant call format and VCFtools. *Bioinformatics, 27*(15), 2156–2158. https://doi.org/10.1093/bioinformatics/btr330

Denoeud, F., Carretero-paulet, L., Dereeper, A., Droc, G., Guyot, R., Pietrella, M., … Lashermes, P. (2014). The coffee genome provides insight into the convergent evolution of caffeine biosynthesis. *Science, 345*(6201), 1181–1184. https://doi.org/10.1126/science.1255274

Hamon, P., Grover, C.E., Davis, A.P., Rakotomalala, J.J., Raharimalala, N.E., Albert, V.A., … Guyot, R. (2017). Genotyping-by-sequencing provides the first well-resolved phylogeny for coffee (*Coffea*) and insights into the evolution of caffeine content in its species GBS coffee phylogeny and the evolution of caffeine content. *Molecular Phylogenetics and Evolution 109*, 351–361. https://doi.org/10.1016/j.ympev.2017.02.009

Huson, D. H. & Bryant, D. (2006). Application of Phylogenetic Networks in Evolutionary Studies. *Molecular Biology and Evolution, 23*(2), 254-267. https://doi.org/10.1093/molbev/msj030

Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows–Wheeler transform. *Bioinformatics, 25*(14), 1754–1760. https://doi.org/10.1093/bioinformatics/btp324

Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., … 1000 Genome Project Data Processing Subgroup. (2009). The Sequence Alignment/Map format and SAMtools. *Bioinformatics, 25*(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352

Mckenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., … DePristo, M. A. (2010). The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Research, 20*(9), 1297–1303. https://doi.org/10.1101/gr.107524.110.20

Paradis, E. & Schliep, K. (2018). ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. *Bioinformatics 35*(3): 526-528. https://doi.org/10.1093/bioinformatics/bty633

R Core Team (2020). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/

Schaumont, D. (2020) GBprocesS: Genotyping-by-Sequencing Data Processing Toolkit [Online]. Available online at https://gitlab.com/dschaumont/GBprocesS

Zhang, J., Kobert, K., Flouri, T., & Stamatakis, A. (2014). PEAR: a fast and accurate Illumina Paired-End reAd mergeR. *Bioinformatics, 30*(5), 614–620. https://doi.org/10.1093/bioinformatics/btt593
