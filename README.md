# Metacongo_Paper
![alt text](https://github.com/UMR-CNRS-7267/Metacongo_Paper/blob/a32eebb3d65b6e58c5f7adf7758bbef7ad969aa8/metacongo_figues/00_logo.png)

This repository contains descriptions of the tools and methods used to analyze metagenomic data from surface water sampled from the city of Pointe-Noire, Congo. The paper is published in XXX and  available here: http:blabla.com



# Prerequisites
To reproduce these analyses we need the following tools installed in your local machine (Linux)
(For our side we used : CentOS Linux release 7.5.1804 (Core) Maipo) 

## Tools and scripts
**bash** available by deafult in all linux distritions (version 4.2.46(2)-release )

**FastQC** available at https://github.com/s-andrews/FastQC (Version used 0.11.8)

**fastp**  available at https://github.com/OpenGene/fastp    (Version used 0.21.0)

**Bowtie 2** available at https://github.com/BenLangmead/bowtie2 (Version used  XXXX)

**Kaiju** available at https://github.com/bioinformatics-centre/kaiju (Version used  XXXX)

**Kraken 2** available at https://github.com/DerrickWood/kraken2 (Version used  XXXX)

**KronaTools** available at https://github.com/marbl/Krona (Version used  XXXX)

**Megahit** available at  https://github.com/voutcn/megahit (Version used  XXXX)

**metaSPAdes** available at https://github.com/ablab/spades (Version used  XXXX)

>some modules from MetaWRAP  Pipleine (but you have to install the whole pipeline if you want to use this ....)

**MetaWRAP** availbale  at https://github.com/bxlab/metaWRAP  (Version used  XXXX)



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


