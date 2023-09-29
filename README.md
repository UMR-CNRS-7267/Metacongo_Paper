# Metacongo_Paper

![](https://github.com/UMR-CNRS-7267/Metacongo_Paper/blob/a32eebb3d65b6e58c5f7adf7758bbef7ad969aa8/metacongo_figues/00_logo.png)

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


> If you choose to keep the original filenames, as downloaded from ncbi, feel free to track that change over the whole pipeline

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
# Raw data QC (Using fastqc)
 ```bash
$ pwd
/home/user_name/RAW_DATA

# Call fastqc program on your fastq files
# In standard Linux (Pc or simple server)
$ fastqc *gz -t X # when X is the number of threads you want to use 

# If you have a cluster with slurm (see Scripts folder for a script named fastqc_slurm.sh)

```

# Trimming and filtering

As mentionned in materiel & methods section of the paper, we used fastp to trim and filter the raw reads.

Here the script used for that purpose:

```bash
# Create a folder to hold the filtered data
$ mkdir  FASTP_FILTERING && cd FASTP_FILTERING
# Symlink files here ...
$ for  i  in ../RAW_DATA/*fg.gz; do ln -$i; done
# Check if OK
$ ll
lrwxrwxrwx 1 foo users         33 Jun 20 14:03 EPNC_R1.fq.gz -> ../RAW_DATA/EPNC_R1.fq.gz
lrwxrwxrwx 1 foo users         33 Jun 20 14:03 EPNC_R2.fq.gz -> ../RAW_DATA/EPNC_R2.fq.gz
``` 
For removing adaptors and filter by quality, please see metacongo_fastp_slurm.sh in Script folder

The slurm script used, with this --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

```bash
# Run the script
$ sbatch metacongo_fastp_slurm.sh
```
How to use fastp if you do not have slurm or HPC

```bash
# Fastp cmd line 
$ fastp -i EPNC_R1.fq.gz -o EPNC_trim_R1.fq.gz -I EPNC_R2.fq.gz -O EPNC_trim_R2.fq.gz --unpaired1 EPNC_orphan_1.fq.gz --unpaired2 EPNC_orphan_2.fq.gz  -z 4 --trim_poly_g --trim_poly_x --detect_adapter_for_pe -l 50 -c -p -h  metacongo_fastp_report.html--json metacongo_fastp_report.json --overrepresentation_analysis  -w 2   --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA  --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

```
Doing sanity check to check the output of fastp (get and idea about the filtering process by compting before/after )

```bash
############ Before fastp #################
EPNC_R1.fq.gz  84886827
EPNC_R2.fq.gz  84886827

TOTAL REAS COUNT BEFORE FASTP¨
>>> 84886827+84886827=== 169773654

############ After fastp #################

# Count was done as previous (using grep -c "@GWNJ-")
____________________________________
EPNC_orphan_1.fq.gz     4678717
EPNC_orphan_2.fq.gz     49618
EPNC_trim_R1.fq.gz      80014112
EPNC_trim_R2.fq.gz      80014112
____________________________________

TOTAL READS COUNT AFTER FASTP
>>> 4678717+49618+80014112+80014112===161779130
REMOVED READS >>> 169773654-164756559=5017095

```

Optional but it is always good to see before/after (we will count the adaptors) in files before/after

```bash
$ for i in *gz; do printf $i"\t"; zcat $i |grep -c "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"; done
EPNC_orphan_1.fq.gz     0
EPNC_orphan_2.fq.gz     0
EPNC_R1.fq.gz   157789 (original file)
EPNC_R2.fq.gz   0
EPNC_trim_R1.fq.gz      1
EPNC_trim_R2.fq.gz      0

$ for i in *gz; do printf $i"\t"; zcat $i |grep -c "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"; done
EPNC_orphan_1.fq.gz     0
EPNC_orphan_2.fq.gz     0
EPNC_R1.fq.gz   0
EPNC_R2.fq.gz   106984 (Original file)
EPNC_trim_R1.fq.gz      0
EPNC_trim_R2.fq.gz      3

```

QC the data again (Optional since fastp has report with before after)

```bash
$ mkdir QC_FILTERED_DATA && cd QC_FILTERED_DATA
#Symlink filtered data as usual (using ln -s)
$ ll
lrwxrwxrwx 1 foo users     40 Jan 22  2021 EPNC_orphan_1.fq.gz -> ../4.FASTP_FILTERING/EPNC_orphan_1.fq.gz
lrwxrwxrwx 1 foo users     40 Jan 22  2021 EPNC_orphan_2.fq.gz -> ../4.FASTP_FILTERING/EPNC_orphan_2.fq.gz
lrwxrwxrwx 1 foo users     39 Jan 22  2021 EPNC_trim_R1.fq.gz -> ../4.FASTP_FILTERING/EPNC_trim_R1.fq.gz
lrwxrwxrwx 1 foo users     39 Jan 22  2021 EPNC_trim_R2.fq.gz -> ../4.FASTP_FILTERING/EPNC_trim_R2.fq.gz

#fastqc on these files
$ fastqc *gz -t X # when X is the number of threads you want to use 

# If you have a cluster with slurm (see Scripts folder for a script named fastqc_filtered_data_slurm.sh)

```










