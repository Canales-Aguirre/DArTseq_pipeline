# DArTseq_pipeline
## UNDER DEVELOPMENT - WORK IN PROGRESS

Scripts for DArTseq_pipeline

Identify_barcodes.py for writing a fasta file containing the barcode for a sample. Identify_barcodes_loop is a bash script that loops over all FASTQ.gz files in the directory and writes a new line for every sample in the fasta output (TableOfBarcodes.fasta).

Gbprocess_SE.ini contains the specifics for the GBprocesS analysis. Files need to have the same name length, and that length needs to be specified in the .ini file. 

The DArTseq_pipeline.txt is a bash script that runs multiple commands for every FASTQ.gz file in the directory.

Output are VCF files.

![DArTseq_pipeline](https://user-images.githubusercontent.com/91179859/153420119-d8b21336-0e79-4842-8f80-85ed7b449bf8.png)
![DArTseq_pipeline_2](https://github.com/sanderdebacker/DArTseq_pipeline/blob/main/DArTseq_pipeline_2.png?raw=true)
