# Metacongo_Paper

![](https://github.com/UMR-CNRS-7267/Metacongo_Paper/blob/a32eebb3d65b6e58c5f7adf7758bbef7ad969aa8/metacongo_figues/00_logo.png)

This repository contains descriptions of the tools and methods used to analyze metagenomic data from surface water sampled from the city of Pointe-Noire, Congo. The paper is published in XXX and  available here: http:blabla.com



# Prerequisites
To reproduce these analyses we need the following tools installed in your local machine (Linux)
(For our side we used : CentOS Linux release 7.5.1804 (Core) Maipo) 

# Tools and scripts


**bash** available by deafult in all linux distritions (version 4.2.46(2)-release )

**FastQC** available at https://github.com/s-andrews/FastQC (Version used 0.11.8)

**fastp**  available at https://github.com/OpenGene/fastp    (Version used 0.21.0)

**Bowtie 2** available at https://github.com/BenLangmead/bowtie2 (Version used  XXXX)

**Kaiju** available at https://github.com/bioinformatics-centre/kaiju (Version used  XXXX)

**Kraken 2** available at https://github.com/DerrickWood/kraken2 (Version used  XXXX)

**KronaTools** available at https://github.com/marbl/Krona (Version used  XXXX)

**seqtk** available at https://github.com/lh3/seqtk (Vresion XXX)

**seqkit** available at https://github.com/shenwei356/seqkit/ (Version XXX)

**Samtools** available at https://github.com/samtools/samtools (Version used v1.9)

**bamtools** available at https://github.com/pezmaster31/bamtools (Version used 2.5.1)

**Megahit** available at  https://github.com/voutcn/megahit (Version used  XXXX)

**metaSPAdes** available at https://github.com/ablab/spades (Version used  XXXX)

>some modules from MetaWRAP  Pipleine (but you have to install the whole pipeline if you want to use this ....)

**MetaWRAP** availbale  at https://github.com/bxlab/metaWRAP  (Version used  XXXX)


# Databases 

**kaiju databases** 

(Feel free to use recent databases if you want), Newer versions are available from kaiju webserver https://kaiju.binf.ku.dk/server

>> Versions Used In the Paper are :

kaiju_db_nr_euk_2022-03-10 (71GB) available at https://kaiju-idx.s3.eu-central-1.amazonaws.com/2022/kaiju_db_nr_euk_2022-03-10.tgz

kaiju_db_plasmids_2022-04-10 (966MB) available at https://kaiju-idx.s3.eu-central-1.amazonaws.com/2022/kaiju_db_plasmids_2022-04-10.tgz

kaiju_db_rvdb_2022-04-07 (983MB)  available at  https://kaiju-idx.s3.eu-central-1.amazonaws.com/2022/kaiju_db_rvdb_2022-04-07.tgz


>> NB: for Kaiju databses we need 3 files for each db: database.fmi, nodes.dmp and names.dmp


**kraken2 databases**





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

# If you have a cluster with slurm (see Scripts folder for a script named fastqc_filtered_data_slurm.sh) and sbatch fastqc_filtered_data_slurm.sh to run it

```
# Human sequences/reads removal (in silico decontamination)

```bash
# Folers/subfolders
$ mkdir HUMAN_CONTA_REMOVAL && cd HUMAN_CONTA_REMOVAL
$ mkdir REF && cd REF

#Bowtie  index was downnloaded from here
https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip

$ wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip
#Extact and delete unwatned folders, files
$ unzip GRCh38_noalt_as.zip
$ mv GRCh38_noalt_as/GRCh38_noalt_as.* .
$ rm -rf GRCh38_noalt_as GRCh38_noalt_as.zip
$ cd ../
#Symlink fastq files
$ for i in ../QC_FILTERED_DATA/*gz; do ln -s $i ; done
$ ll
drwxr-xr-x 2 foo users 216 Nov 21 10:48 6.1.REF
lrwxrwxrwx 1 foo users  41 Nov 21 10:58 EPNC_orphan_1.fq.gz -> ../5.QC_FILTERED_DATA/EPNC_orphan_1.fq.gz
lrwxrwxrwx 1 foo users  41 Nov 21 10:58 EPNC_orphan_2.fq.gz -> ../5.QC_FILTERED_DATA/EPNC_orphan_2.fq.gz
lrwxrwxrwx 1 foo users  40 Nov 21 10:58 EPNC_trim_R1.fq.gz -> ../5.QC_FILTERED_DATA/EPNC_trim_R1.fq.gz
lrwxrwxrwx 1 foo users  40 Nov 21 10:58 EPNC_trim_R2.fq.gz -> ../5.QC_FILTERED_DATA/EPNC_trim_R2.fq.gz

# Get data for input (Note orphan reads from R1 and R2 are merged with simple cat command)
$ pwd 
/home/foo/Metacongo_Paper/HUMAN_CONTA_REMOVAL

# Symlink files as usual  (using ln -s)
$ for i in ../QC_FILTERED_DATA/*gz; do ln -s $i ; done
$ ll
drwxr-xr-x 2 foo users        216 Nov 21  2022 6.1.REF
lrwxrwxrwx 1 foo users         41 Nov 21  2022 EPNC_orphan_1.fq.gz -> ../5.QC_FILTERED_DATA/EPNC_orphan_1.fq.gz
lrwxrwxrwx 1 foo users         41 Nov 21  2022 EPNC_orphan_2.fq.gz -> ../5.QC_FILTERED_DATA/EPNC_orphan_2.fq.gz
-rw-r--r-- 1 foo users  468447598 Nov 21  2022 EPNC_orphan.fq.gz #This file is from cat EPNC_orphan_1.fq.gz EPNC_orphan_2.fq.gz> EPNC_orphan.fq.gz 
lrwxrwxrwx 1 foo users         40 Nov 21  2022 EPNC_trim_R1.fq.gz -> ../QC_FILTERED_DATA/EPNC_trim_R1.fq.gz
lrwxrwxrwx 1 foo users         40 Nov 21  2022 EPNC_trim_R2.fq.gz -> ../QC_FILTERED_DATA/EPNC_trim_R2.fq.gz

# TO FACILITATE THIS STEP WE HAVE TO EXPORT SOME VAR DIRECTLY FROM CMD (IN BASH)
$ F_READS='EPNC_trim_R1.fq.gz'
$ R_READS='EPNC_trim_R2.fq.gz'
$ ORPHAN_READS='EPNC_orphan.fq.gz'
$ HUMAN_BW2_INDEX='REF/GRCh38_noalt_as'

$ OUT_SAM='Read_vs_human.sam'
$ OUT_BAM='Read_vs_human.bam'
$ OUT_UNMAPPED_BAM='Unmapped.bam'

$ OUT_BAM_SORTED='Unmapped_sorted.bam'
$ UNMAPPED_LIST='Unmapped.list'



# Run bowtie2 to map reads to human genome index

$ bowtie2 --very-sensitive-local -x $HUMAN_BW2_INDEX -1 $F_READS  -2 $R_READS -U $ORPHAN_READS -S $OUT_SAM -p $SLURM_CPUS_PER_TASK

# Converting sam to bam 

$ samtools view -S -bh $OUT_SAM > $OUT_BAM

# Sort bam file

$ samtools sort $OUT_BAM -o $OUT_BAM_SORTED

# Getting some stats from the whole BAM file  (Optional)

$ bamtools stats -in $OUT_BAM_SORTED >mapping_stat_from_whole_bam.stats

# Getting unmapped reads from the bam file ....

$ samtools view -b -f 4 $OUT_BAM_SORTED > $OUT_UNMAPPED_BAM

# This line produce the unmapped list BUT with duplication
# Have to remove duplicated reads at the end using seqkit

$ samtools view $OUT_UNMAPPED_BAM  |awk '{print $1}' >$UNMAPPED_LIST

# Using seqtk to get unmapped reads in fastq format ...

$ seqtk subseq $F_READS  $UNMAPPED_LIST  | gzip > EPNC_trim_no_human_R1.fq.gz
$ seqtk subseq $R_READS  $UNMAPPED_LIST  | gzip > EPNC_trim_no_human_R2.fq.gz
$ seqtk subseq $ORPHAN_READS   $UNMAPPED_LIST | gzip   > EPNC_orphan_no_human.fq.gz

# Rename files (Optional)
$ mv EPNC_trim_no_human_R1.fq.gz  EPNC_trim_ready_R1.fq.gz
$ mv EPNC_trim_no_human_R2.fq.gz  EPNC_trim_ready_R2.fq.gz
$ mv EPNC_orphan_no_human.fq.gz   EPNC_orphan_ready.fq.gz

# Puts these files in another folder
$ mkdir READY_FASTQ_FILES
mv EPNC_trim_ready_R1.fq.gz EPNC_trim_ready_R2.fq.gz EPNC_orphan_ready.fq.gz 6.2.READY_FASTQ_FILES

# If you have a cluster with slurm (see Scripts folder for a script named fastqc_filtered_data_slurm.sh) and sbatch bowtie2_vs_human_slurm.sh to run it
```


>> At this step data are ready to analyse, we will profile the metagenomic  read using two profiler (Kaiju and Kraken2).

>> kaiju use amino acid database, so there is six-frames translation of our reads then compare to proteins DB

>> kraken2 use  nucleotide database. 

>> All reads those are not assigned/classified by kaiju were passed to kraken2

>> ***This step  needs databases so be sur to have enough place in your local server ...***








