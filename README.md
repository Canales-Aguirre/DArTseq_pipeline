# DArTseq_pipeline
## UNDER DEVELOPMENT - WORK IN PROGRESS

Scripts for DArTseq_pipeline

Identify_barcodes_SdB.py is an adapted script for writing a fasta file containing the barcode of all FASTQ(.gz) files in the directory (TableOfBarcodes.fasta).

Gbprocess_SE.ini contains the specifics for the GBprocesS analysis. Check the available online manual to change these settings to your analysis. **Important:** files need to have the same name length, and that length also needs to be specified in the .ini file. 

The DArTseq_pipeline.py is a python script that runs multiple commands for every FASTQ.gz file in the directory.

Output are VCF files.

Broad description available in the DArTseq_pipeline_Manual pdf file.

![DArTseq_pipeline](https://github.com/sanderdebacker/DArTseq_pipeline/blob/main/Images/DArTseq_pipeline.png?raw=true)
