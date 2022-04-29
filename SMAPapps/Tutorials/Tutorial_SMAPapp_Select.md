# Tutorial SMAPapp-Select

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
The aim of this tutorial is to create two new haplotypes tables based on the original table. The first table only consists of samples from the species *Coffea arabica*, whereas the second table consists of a subset of loci with a high discriminatory power among samples from the 24 species in the original haplotypes table.

## Subset samples
To select the *C. arabica* samples in the original haplotypes table, a --samples list was created with only names from *C. arabica* samples. Sample names that were not listed in the --samples list are removed from the haplotypes table. The order of the sample names in the --samples list defines the order of the samples in the resulting matrix. The retained *C. arabica* samples were also renamed so that they had more convenient names in the new table. The --samples list (Sample_names_Carabica.txt) is a tab-delimited text file containing the original sample names (as listed in the haplotypes table) in the first column and the new sample names in the second column, as shown below (Table 1). We chose to label each sample by its species identity (i.e. *Coffea arabica*) followed by a decimal number. Technical GBS library replicates are indicated by the same unit before the decimal point but a different tenth after the decimal point (e.g. 1.1 and 1.2). Biological replicates of *C. arabica* are indicated by a different ID number (e.g. 1.1 and 2.1).

**Table 1.** List of sample names (left column = original sample name, right column = new sample name) in the haplotypes table belonging to the species *C. arabica*.
|     |     |
| --- | --- |
| SRR11294243 | Coffea arabica 1.1            |
| SRR11294283 | Coffea arabica 1.2            |
| SRR11294274 | Coffea arabica 2.1            |
| SRR11294273 | Coffea arabica 3.1            |

The user can also choose to rename some samples while keeping the original sample name of other samples. If, for example, the user would only like to retain the original sample name of the first *C. arabica* sample, the following --samples list can be used (Table 2).

**Table 2.** Example of a list of sample names (left column = original sample name, right column = new sample name) in the haplotypes table belonging to the species *C. arabica*. The first sample in the list with ID SRR11294243 keeps its original sample ID, while the remaining three samples are renamed.
|     |     |
| --- | --- |
| SRR11294243 | <NA>               |
| SRR11294283 | Coffea arabica 1.2 |
| SRR11294274 | Coffea arabica 2.1 |
| SRR11294273 | Coffea arabica 3.1 |

The command line parameters are listed below.

`python3 SMAPapp-Select.py --table haplotypes_R10_F1_call.txt --samples Sample_names_Coffea_arabica.txt`

The output haplotypes table consisted of 808 haplotypes in 292 loci that were polymorphic in the four *C. arabica* samples.

## Subset loci

The original haplotypes table was filtered for the haplotype calls in a subset of six loci (Table 3). These loci were selected with the SMAPapp-Matrix script because they had unique haplotypes for at least 9 out of 24 species (see the *inference of locus information* section in the [Tutorial of SMAPapp-Matrix](Tutorial_SMAPapp_Matrix.md)).

**Table 3.** List of loci with only unique haplotypes in at least 9 out of 24 species in the dataset. The locus ID is shown in the first column.
|     |
| --- |
| chr1:30787505-30787768_+ |
| chr2:1006717-1006879_+   |
| chr3:2510192-2510423_+   |
| chr3:3816695-3816946_+   |
| chr7:20976690-20976906_+ |
| chr10:2239337-2239448_+  |

All samples were also renamed so that they have a more convenient sample name (Table 4). The samples were ordered according to their assumed phylogenetic relationships (Hamon *et al*., 2017).

**Table 4.** List with sample names (left column = original sample name, right column = new sample name) that was parsed to SMAPapp-Select script via the --samples option.
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

The following command line options were used.

`python3 SMAPapp-Select.py --table haplotypes_R10_F1_call.txt --samples Sample_names.txt --loci Loci.txt`

The resulting haplotypes table comprised 144 haplotypes in six loci that were called in 49 samples.

## References
Bawin, Y., Ruttink, T., Staelens, A., Haegeman, A., Stoffelen, P., Mwanga Mwanga, J.-C.I., … Janssens, S.B. (2021). Phylogenomic analysis clarifies the evolutionary origin of *Coffea arabica*. *Journal of Systematics and Evolution, 59*(5), 953-963. https://doi.org/10.1111/jse.12694

Danecek, P., Auton, A., Abecasis, G., Albers, C. A., Banks, E., DePristo, M. A., … Durbin, R. (2011). The variant call format and VCFtools. *Bioinformatics, 27*(15), 2156–2158. https://doi.org/10.1093/bioinformatics/btr330

Denoeud, F., Carretero-paulet, L., Dereeper, A., Droc, G., Guyot, R., Pietrella, M., … Lashermes, P. (2014). The coffee genome provides insight into the convergent evolution of caffeine biosynthesis. *Science, 345*(6201), 1181–1184. https://doi.org/10.1126/science.1255274

Hamon, P., Grover, C.E., Davis, A.P., Rakotomalala, J.J., Raharimalala, N.E., Albert, V.A., … Guyot, R. (2017). Genotyping-by-sequencing provides the first well-resolved phylogeny for coffee (*Coffea*) and insights into the evolution of caffeine content in its species GBS coffee phylogeny and the evolution of caffeine content. *Molecular Phylogenetics and Evolution 109*, 351–361. https://doi.org/10.1016/j.ympev.2017.02.009

Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows–Wheeler transform. *Bioinformatics, 25*(14), 1754–1760. https://doi.org/10.1093/bioinformatics/btp324

Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., … 1000 Genome Project Data Processing Subgroup. (2009). The Sequence Alignment/Map format and SAMtools. *Bioinformatics, 25*(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352

Mckenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., … DePristo, M. A. (2010). The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Research, 20*(9), 1297–1303. https://doi.org/10.1101/gr.107524.110.20

Schaumont, D. (2020) GBprocesS: Genotyping-by-Sequencing Data Processing Toolkit [Online]. Available online at https://gitlab.com/dschaumont/GBprocesS"

Zhang, J., Kobert, K., Flouri, T., & Stamatakis, A. (2014). PEAR: a fast and accurate Illumina Paired-End reAd mergeR. *Bioinformatics, 30*(5), 614–620. https://doi.org/10.1093/bioinformatics/btt593
