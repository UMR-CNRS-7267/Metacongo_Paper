# Metacongo_Paper


This repository contains descriptions of the tools and methods used to analyze metagenomic data from surface water sampled from the city of Pointe-Noire, Congo. The paper is published in XXX and  available here: http:blabla.com



# Prerequisites
To reproduce these analyses we need the following tools installed in your local machine (Linux)
(For our side we used : CentOS Linux release 7.5.1804 (Core) Maipo) 

**FastQC** available at : https://github.com/s-andrews/FastQC (Version used 0.11.8)

**fastp**  available at : https://github.com/OpenGene/fastp    (Version used 0.21.0)



# Getting the raw data 

* The technology used is Illumina, in paired end mode.
* The data are made available in Sequence Read Archive from NCBI and can be requested using these accessions numbers: XXXXX, YYYY
*  Feel free to use your favorite tools to get these files in local.
* If you want to reproduce these analysis, please renames these files according to this manual:

# Renaming data files (Optional)


```bash
# make a directory to store the raw data 
$ mkdir RAW_DATA
$ mv XXXX_R1.fastq.gz EPNC_R1.fq.gz
$ mv XXXX_R2.fastq.gz EPNC_R2.fq.gz

# Get Reads count (Optional) from the fastq files
$ for i in *gz; do  printf $i"\t"; gzip -cd $i | grep -c "@GWNJ-"; done
EPNC_R1.fq.gz  84886827
EPNC_R2.fq.gz  84886827


``` 


