# DArTseq_pipeline
## UNDER DEVELOPMENT - WORK IN PROGRESS

Hub containing information and scripts for the DArTseq analysis pipeline.

### DArTseq_pipeline
Contains the scripts and configuration files to process raw DArTseq FASTQ files into VCF files.
The DArTseq_pipeline.py is a python script that runs multiple commands for every FASTQ.gz file in the directory.
Gbprocess_SE.ini contains the specifics for the GBprocesS analysis. Check the available online manual to change these settings to your analysis. **Important:** files need to have the same name length, and that length also needs to be specified in the .ini file.
Broader description available in the **DArTseq_pipeline_Manual.pdf** file.

*Overview of the DArTseq pipeline: from raw FASTQ.gz files to VCF files.*
<img src="https://github.com/sanderdebacker/DArTseq_pipeline/blob/main/Images/DArTseq_pipeline.png?raw=true" />

### Information
Offers additional information about the Variant Call Format and installation guides for some tools (fastStructure).

### Phylogenetics 
Contains script for filtering reads based on missing data and merging VCF files into one multi-sample VCF file.
For MrBayes, RAxML-NG, and IQTree there are simple bash loop-scripts to automate analysis for multiple alignments.

### Rscripts
Different script for generating graphs, plots, newick files.

### SMAPapps
Third-party Python scripts to analyse SMAP output.
All credits to the creator(s).

### VcfHunter 
Third-party toolbox to process NGS and variant call data.
Some scripts were slightly adjusted to fit the available computing system.
All credits to the creator(s).

### DArTseq_pipeline_Manual.pdf
More information about the scripts.
Explains in detail multiple analyses performed on DArTseq data.

### Renamer.xlsm
Excel worksheet using Macros to rename files en masse.
