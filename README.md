# Metacongo_Paper

![](https://github.com/UMR-CNRS-7267/Metacongo_Paper/blob/a32eebb3d65b6e58c5f7adf7758bbef7ad969aa8/metacongo_figues/00_logo.png)



 **_The repository contains descriptions of the tools and methods used to analyze metagenomic data from surface water sampled from the city of Pointe-Noire, Congo. These analyses were published in this paper :_**



# Prerequisites
To reproduce these analyses you need the following tools installed in your local machine (Linux)
(For our side we used : CentOS Linux release 7.5.1804 (Core) Maipo) ....

# Tools and scripts


**bash** available by deafult in all linux distritions (version used  4.2.46(2)-release )

**FastQC** available at https://github.com/s-andrews/FastQC (Version used 0.11.8)

**fastp**  available at https://github.com/OpenGene/fastp (Version used 0.23.0)

**Bowtie 2** available at https://github.com/BenLangmead/bowtie2 (Version used  2.3.4.3)

**Kaiju** available at https://github.com/bioinformatics-centre/kaiju (Version used  v1.9.0)

**Kraken2** available at https://github.com/DerrickWood/kraken2 (Version used  2.1.2)

**KronaTools** available at https://github.com/marbl/Krona (Version used  v 2.8.1)

**seqtk** available at https://github.com/lh3/seqtk (Vresion v 1.3-r106)

**seqkit** available at https://github.com/shenwei356/seqkit/ (Version used v2.5.1 )

**Samtools** available at https://github.com/samtools/samtools (Version used v1.9)

**bamtools** available at https://github.com/pezmaster31/bamtools (Version used 2.5.1)

**Megahit** available at  https://github.com/voutcn/megahit (Version used  1.2.9)

**MetaWRAP** availbale  at https://github.com/bxlab/metaWRAP (Version used  v1.1.2)

**gtdb** available at https://github.com/Ecogenomics/GTDBTk (version used 1.4.1)

**MMseqs2** available at https://github.com/soedinglab/MMseqs2 (Version used release 14-7e284) 

**metabat2** available at 2.12.1 

**Maxbin2**  available at 2.2.7

>some modules from MetaWRAP  Pipleine (but you have to install the whole pipeline if you want to use this ....)

**metaSPAdes** available at https://github.com/ablab/spades (Version used  3.15.2)


# Databases 

**kaiju databases** 

Please, feel free to use recent databases if you want , Newer versions are available from kaiju webserver https://kaiju.binf.ku.dk/server

>> Versions Used In the Paper are :

kaiju_db_nr_euk_2022-03-10 (71GB) available at https://kaiju-idx.s3.eu-central-1.amazonaws.com/2022/kaiju_db_nr_euk_2022-03-10.tgz

kaiju_db_plasmids_2022-04-10 (966MB) available at https://kaiju-idx.s3.eu-central-1.amazonaws.com/2022/kaiju_db_plasmids_2022-04-10.tgz

kaiju_db_rvdb_2022-04-07 (983MB)  available at  https://kaiju-idx.s3.eu-central-1.amazonaws.com/2022/kaiju_db_rvdb_2022-04-07.tgz


>> NB: for Kaiju databses we need 3 files for each db: database.fmi, nodes.dmp and names.dmp


**kraken2 databases**

The Version used in this study is (PlusPF) which means Standard plus Refeq protozoa & fungi available at https://benlangmead.github.io/aws-indexes/k2,
we used March 14, 2023 version but feel free to use recent version if availbale ...

**Human genome reference index for bowtie2  decontamination**

The index used is available here : https://benlangmead.github.io/aws-indexes/bowtie





# Getting the raw data 

* The technology used is Illumina, in paired end mode (2X150).
* The data are made available in Sequence Read Archive from NCBI from bioproject : PRJNA1021800
* Feel free to use your favorite tools to get these files in local.
* If you want to reproduce these analysis, please renames these files according to this manual:


```bash

# Feel free to organize your project as you want ...
$ mkdir -p workspaces/Metacongo_Paper  && cd workspaces/Metacongo_Paper

```
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
lrwxrwxrwx 1 foo users  41 Nov 21 10:58 EPNC_orphan_1.fq.gz -> ../QC_FILTERED_DATA/EPNC_orphan_1.fq.gz
lrwxrwxrwx 1 foo users  41 Nov 21 10:58 EPNC_orphan_2.fq.gz -> ../QC_FILTERED_DATA/EPNC_orphan_2.fq.gz
lrwxrwxrwx 1 foo users  40 Nov 21 10:58 EPNC_trim_R1.fq.gz -> ../QC_FILTERED_DATA/EPNC_trim_R1.fq.gz
lrwxrwxrwx 1 foo users  40 Nov 21 10:58 EPNC_trim_R2.fq.gz -> ../QC_FILTERED_DATA/EPNC_trim_R2.fq.gz

# Get data for input (Note orphan reads from R1 and R2 are merged with simple cat command)
$ pwd 
/home/foo/Metacongo_Paper/HUMAN_CONTA_REMOVAL

# Symlink files as usual  (using ln -s)
$ for i in ../QC_FILTERED_DATA/*gz; do ln -s $i ; done
$ ll
drwxr-xr-x 2 foo users        216 Nov 21  2022 6.1.REF
lrwxrwxrwx 1 foo users         41 Nov 21  2022 EPNC_orphan_1.fq.gz -> ../QC_FILTERED_DATA/EPNC_orphan_1.fq.gz
lrwxrwxrwx 1 foo users         41 Nov 21  2022 EPNC_orphan_2.fq.gz -> ../QC_FILTERED_DATA/EPNC_orphan_2.fq.gz
-rw-r--r-- 1 foo users  468447598 Nov 21  2022 EPNC_orphan.fq.gz #This file is from cat EPNC_orphan_1.fq.gz EPNC_orphan_2.fq.gz> EPNC_orphan.fq.gz 
lrwxrwxrwx 1 foo users         40 Nov 21  2022 EPNC_trim_R1.fq.gz -> ../QC_FILTERED_DATA/EPNC_trim_R1.fq.gz
lrwxrwxrwx 1 foo users         40 Nov 21  2022 EPNC_trim_R2.fq.gz -> ../QC_FILTERED_DATA/EPNC_trim_R2.fq.gz

# TO FACILITATE THIS STEP WE HAVE TO EXPORT SOME VAR DIRECTLY FROM CMD (IN BASH)
$ F_READS='EPNC_trim_R1.fq.gz'
$ R_READS='EPNC_trim_R2.fq.gz'
$ ORPHAN_READS='EPNC_orphan.fq.gz'
$ HUMAN_BW2_INDEX='REF/GRCh38_noalt_as'

OUT_SAM='Read_vs_human.sam'
OUT_BAM='Read_vs_human.bam'
OUT_UNMAPPED_BAM='Unmapped.bam'

OUT_BAM_SORTED='Unmapped_sorted.bam'
UNMAPPED_LIST='Unmapped.list'




# Run bowtie2 to map reads to human genome index

$ bowtie2 --very-sensitive-local -x REF/GRCh38_noalt_as -1 EPNC_trim_R1.fq.gz  -2 EPNC_trim_R2.fq.gz -U EPNC_orphan.fq.gz -S Read_vs_human.sam -p 4

# Converting sam to bam 

$ samtools view -S -bh Read_vs_human.sam > Read_vs_human.bam

# Sort bam file

$ samtools sort Read_vs_human.bam -o Read_vs_human_sorted.bam

# Getting some stats from the whole BAM file  (Optional)

$ bamtools stats -in Read_vs_human_sorted.bam >mapping_stat_from_whole_bam.stats

# Getting unmapped reads from the bam file ....

$ samtools view -b -f 4 Read_vs_human_sorted.bam > Unmapped.bam

# This line produce the unmapped list BUT with duplication
# Have to remove duplicated reads at the end using seqkit

$ samtools view  Unmapped.bam  |awk '{print $1}' > Unmapped.list

# Using seqtk to get unmapped reads in fastq format ...

$ seqtk subseq EPNC_trim_R1.fq.gz   Unmapped.list  | gzip > EPNC_trim_no_human_R1.fq.gz
$ seqtk subseq EPNC_trim_R2.fq.gz   Unmapped.list  | gzip > EPNC_trim_no_human_R2.fq.gz
$ seqtk subseq EPNC_orphan.fq.gz    Unmapped.list | gzip   > EPNC_orphan_no_human.fq.gz

# Rename files (Optional)
$ mv EPNC_trim_no_human_R1.fq.gz  EPNC_trim_ready_R1.fq.gz
$ mv EPNC_trim_no_human_R2.fq.gz  EPNC_trim_ready_R2.fq.gz
$ mv EPNC_orphan_no_human.fq.gz   EPNC_orphan_ready.fq.gz

# Puts these files in another folder
$ mkdir READY_FASTQ_FILES
mv EPNC_trim_ready_R1.fq.gz EPNC_trim_ready_R2.fq.gz EPNC_orphan_ready.fq.gz READY_FASTQ_FILES

# If you have a cluster with slurm (see Scripts folder for a script named fastqc_filtered_data_slurm.sh) and sbatch bowtie2_vs_human_slurm.sh to run it
```


>> At this step data are ready to analyse, we will profile the metagenomic  read using two profiler (Kaiju and Kraken2).

>> kaiju use amino acid database, so there is six-frames translation of our reads then compare to proteins DB

>> kraken2 use  nucleotide database. 

>> All reads those are not assigned/classified by kaiju were passed to kraken2

>> ***This step  needs databases so be sur to have enough place in your local server ...***

# PROFILING 

## Kaiju Profiling 

```bash
$ mkdir PROFILING && cd PROFILING
$ mkdir USING_KAIJU USING_KRAKEN2 && cd USING_KAIJU

# Symlink files as usual

$ ln -s ../../HUMAN_CONTA_REMOVAL/READY_FASTQ_FILES_CLEAN/EPNC_orphan_ready_clean.fq.gz
$ ln -s ../../HUMAN_CONTA_REMOVAL/READY_FASTQ_FILES_CLEAN/EPNC_trim_ready_clean_R1.fq.gz 
$ ln -s ../../HUMAN_CONTA_REMOVAL/READY_FASTQ_FILES_CLEAN/EPNC_trim_ready_clean_R2.fq.gz

# For informations (Read count)

EPNC_orphan_ready_clean.fq.gz   4699049
EPNC_trim_ready_clean_R1.fq.gz  79712668
EPNC_trim_ready_clean_R2.fq.gz  79712668

# Run kaiju on the data using databases: kaiju_db_nr_euk_2022-03-10 then the unclassified using this database will be used with 
# kaiju_db_plasmids_2022-04-10 then the unclassfied will be used with kaiju_db_rvdb_2022-04-07 database.
# If you have a cluster with slurm (see Scripts folder for a script named kaiju_profiling_all_slurm__last.sh) and sbatch kaiju_profiling_all_slurm__last.sh to run it 

```

> Using  ############# Using NR_EUK DATABASES #####################

```bash
# Running   kaiju-multi on paired end reads ...
# KAIJU_NR_EUK_DB='/kaiju_db/kaiju_db_nr_euk_2022-03-10/kaiju_db_nr_euk.fmi'
# KAIJU_NR_EUK_NODES='/kaiju_db/kaiju_db_nr_euk_2022-03-10/nodes.dmp'
# KAIJU_NR_EUK_NAMES='/kaiju_db/kaiju_db_nr_euk_2022-03-10/names.dmp'

# F_READS='EPNC_trim_ready_clean_R1.fq.gz'
# R_READS='EPNC_trim_ready_clean_R2.fq.gz'
# ORPHAN_READS='EPNC_orphan_ready_clean.fq.gz'

#Output for the profiling using NR_EUK
# OUTPUT_PE_NR_EUK='Metacongo_kaiju__PE_NR_EUK.out'
# OUTPUT_SE_NR_EUK='Metacongo_kaiju__SE_NR_EUK.out'
# OUTPUT_ALL_NR_EUK='Metacongo_kaiju__ALL_NR_EUK.out'

# After merging the output of PE and SE have to create krona file
# OUTPUT_ALL_NR_EUK_KR='Metacongo_kaiju__ALL_NR_EUK.krona'


$ kaiju-multi -z 4 -E 0.01 -t /kaiju_db/kaiju_db_nr_euk_2022-03-10/nodes.dmp \
                           -f /kaiju_db/kaiju_db_nr_euk_2022-03-10/kaiju_db_nr_euk.fmi \
                           -i EPNC_trim_ready_clean_R1.fq.gz\
                           -j EPNC_trim_ready_clean_R2.fq.gz  > Metacongo_kaiju__PE_NR_EUK.out


# Running   kaiju (not multi) on orphan merged reads ...
$ kaiju -z 4 -E 0.01 -t /kaiju_db/kaiju_db_nr_euk_2022-03-10/nodes.dmp \
                     -f /kaiju_db/kaiju_db_nr_euk_2022-03-10/kaiju_db_nr_euk.fmi \
                     -i EPNC_orphan_ready_clean.fq.gz  > $Metacongo_kaiju__SE_NR_EUK.out


# Combining the output for PE and ORPHAN (SE) ...

$ cat  Metacongo_kaiju__PE_NR_EUK.out Metacongo_kaiju__SE_NR_EUK.out > Metacongo_kaiju__ALL_NR_EUK.out


# Adding full taxa names ... to output ...........

$ kaiju-addTaxonNames -p  -t /kaiju_db/kaiju_db_nr_euk_2022-03-10/nodes.dmp \
                          -n /kaiju_db/kaiju_db_nr_euk_2022-03-10/names.dmp \
                          -i Metacongo_kaiju__ALL_NR_EUK.out  -o Metacongo_kaiju__ALL_NR_EUK.out_with_name.tsv"


# Get read count for classified reads (Optional)

$ echo "TOTAL REDS  COUNT: " $(wc -l  Metacongo_kaiju__ALL_NR_EUK.out |awk '{print $1}')

$ echo "READ CLASSIFIED COUNT: " $(grep -w -c "C" Metacongo_kaiju__ALL_NR_EUK.out ) 


# Converting to kaiju output to krona file

$ kaiju2krona  -t /kaiju_db/kaiju_db_nr_euk_2022-03-10/nodes.dmp \
             -n  /kaiju_db/kaiju_db_nr_euk_2022-03-10/names.dmp \
              -i  Metacongo_kaiju__ALL_NR_EUK.out -o Metacongo_kaiju__ALL_NR_EUK.krona


# Creating html from krona file

$ ktImportText -o Metacongo_kaiju__ALL_NR_EUK.krona.html  Metacongo_kaiju__ALL_NR_EUK.krona

# Creating classification summary for phylum, class, order family, genus and species ...
# using loop in bash

$ for i in phylum class order family genus species; do kaiju2table \
                                            -t /kaiju_db/kaiju_db_nr_euk_2022-03-10/nodes.dmp \
                                            -n /kaiju_db/kaiju_db_nr_euk_2022-03-10/names.dmp \
                                            -r $i -o Metacongo_kaiju__ALL_NR_EUK.out"_"$i"__summary.tsv" Metacongo_kaiju__ALL_NR_EUK.out; done




# At this stage profiling the data using kaiju_db_nr_euk_2022-03-10 is done we will get un_profiled data from this step and proceed with another kaijudb

``` 

>  ################# Using RVDB DATABASE ##################################

> we are going to extract the unclassified reads from the output and re_run on Virus db


#Names after First run for  RVBD classification
F_Unc_from_nr_euk='F_Unc_for_rvbd.fq.gz'
R_Unc_from_nr_euk='R_Unc_for_rvbd.fq.gz'
O_Unc_from_nr_euk='O_Unc_for_rvbd.fq.gz'

#Output for the profiling using RVBD
OUTPUT_PE_RVDB='Metacongo_kaiju__PE_RVDB.out'
OUTPUT_SE_RVDB='Metacongo_kaiju__SE_RVDB.out'
OUTPUT_ALL_RVDB='Metacongo_kaiju__ALL_RVDB.out'

#After merging the output of PE and SE have to create krona file
OUTPUT_ALL_RVDB_KR='Metacongo_kaiju__ALL_RVDB.krona'

KAIJU_RVDB_DB='/home/databases/kaiju_db/kaiju_db_rvdb_2022-04-07/kaiju_db_rvdb.fmi'
KAIJU_RVDB_NODES='/home/databases/kaiju_db/kaiju_db_rvdb_2022-04-07/nodes.dmp'
KAIJU_RVDB_NAMES='/home/databases/kaiju_db/kaiju_db_rvdb_2022-04-07/names.dmp'

```bash

# Getting list (Read ID)

$ grep -w 'U' Metacongo_kaiju__ALL_NR_EUK.out |awk '{print $2}' > Unclassified_from_nr_euk.list

# Extracting Unclassified reads from original fastq

$ seqtk subseq  EPNC_trim_ready_clean_R1.fq.gz  Unclassified_from_nr_euk.list  | gzip > F_Unc_for_rvbd.fq.gz
$ seqtk subseq  REPNC_trim_ready_clean_R2.fq.gz  Unclassified_from_nr_euk.list  | gzip > $R_Unc_for_rvbd.fq.gz
$ seqtk subseq  EPNC_orphan_ready_clean.fq.gz  Unclassified_from_nr_euk.list | gzip > O_Unc_for_rvbd.fq.gz

# Files are ready for second round of profiling 

# Running   kaiju-multi on paired end reads .....

$ kaiju-multi -z 2 -E 0.01 -t /kaiju_db/kaiju_db_rvdb_2022-04-07/nodes.dmp \
                         -f /kaiju_db/kaiju_db_rvdb_2022-04-07/kaiju_db_rvdb.fmi \
                        -i F_Unc_for_rvbd.fq.gz -j R_Unc_for_rvbd.fq.gz  > Metacongo_kaiju__PE_RVDB.out




$ Running  kaiju (not multi) on orphan merged reads ......

$ kaiju -z 4  -E 0.01 -t /kaiju_db/kaiju_db_rvdb_2022-04-07/nodes.dmp \
                    -f /kaiju_db/kaiju_db_rvdb_2022-04-07/kaiju_db_rvdb.fmi \
                    -i O_Unc_for_rvbd.fq.gz  > Metacongo_kaiju__SE_RVDB.out


# Combining the output for PE and ORPHAN (SE) here... 

$ cat Metacongo_kaiju__PE_RVDB.out Metacongo_kaiju__SE_RVDB.out > Metacongo_kaiju__ALL_RVDB.out


#  Adding full taxa names ... to output

$ kaiju-addTaxonNames -p  -t /kaiju_db/kaiju_db_rvdb_2022-04-07/nodes.dmp \
                        -n /kaiju_db/kaiju_db_rvdb_2022-04-07/names.dmp \
                        -i Metacongo_kaiju__ALL_RVDB.out  -o Metacongo_kaiju__ALL_RVDB.out"_with_name.tsv"



# How many reads are classified....

$ echo "TOTAL REDS  COUNT: " $(wc -l  Metacongo_kaiju__ALL_RVDB.out |awk '{print $1}') 

$ echo "READ CLASSIFIED COUNT: " $(grep -w -c "C" Metacongo_kaiju__ALL_RVDB.out ) |tee -a analysis.log


# converting to kaiju output to krona file"

$ kaiju2krona  -t /kaiju_db/kaiju_db_rvdb_2022-04-07/nodes.dmp \
             -n  /kaiju_db/kaiju_db_rvdb_2022-04-07/names.dmp \
             -i Metacongo_kaiju__ALL_RVDB.out -o Metacongo_kaiju__ALL_RVDB.krona



# Creating html from krona file ...

$ ktImportText -o Metacongo_kaiju__ALL_RVDB.krona.html  Metacongo_kaiju__ALL_RVDB.krona


# Creating classification summary for phylum, class, order family, genus and species ..

$ for i in phylum class order family genus species; do kaiju2table  \
                        -t /kaiju_db/kaiju_db_rvdb_2022-04-07/nodes.dmp \
                        -n /kaiju_db/kaiju_db_rvdb_2022-04-07/names.dmp \
                        -r $i -o Metacongo_kaiju__ALL_RVDB.out"_"$i"__summary.tsv" Metacongo_kaiju__ALL_RVDB.out; done




```

>  ################################## Using RVDB Plasmid databse ##################################

Using PL (for plasmids) DATABASE

Going to extract the unclassified reads from  previous analysis output and re_run on Plasmid db


KAIJU_PL_DB='/home/databases/kaiju_db/kaiju_db_plasmids_2022-04-10/kaiju_db_plasmids.fmi'
KAIJU_PL_NODES='/home/databases/kaiju_db/kaiju_db_plasmids_2022-04-10/nodes.dmp'
KAIJU_PL_NAMES='/home/databases/kaiju_db/kaiju_db_plasmids_2022-04-10/names.dmp'


#Names after First run for PLASMID classification
F_Unc_from_rvbd='F_Unc_for_pl.fq.gz'
R_Unc_from_rvbd='R_Unc.for_pl.fq.gz'
O_Unc_from_rvbd='O_Unc.for_pl.fq.gz'

#Output for the profiling using PL
OUTPUT_PE_PL='Metacongo_kaiju__PE_PL.out'
OUTPUT_SE_PL='Metacongo_kaiju__SE_PL.out'
OUTPUT_ALL_PL='Metacongo_kaiju__ALL_PL.out'

#After merging the output of PE and SE have to create krona file
OUTPUT_ALL_PL_KR='Metacongo_kaiju__ALL_PL.krona'




```bash

# getting reads list"

$ grep -w 'U' Metacongo_kaiju__ALL_RVDB.out |awk '{print $2}' > Unclassified_from_rvbd.list

$ # Extracting Unclassified reads from original fastq 

# Extract reads

$ seqtk subseq EPNC_trim_ready_clean_R1.fq.gz  Unclassified_from_rvbd.list  | gzip > F_Unc_for_pl.fq.gz
$ seqtk subseq EPNC_trim_ready_clean_R2.fq.gz   Unclassified_from_rvbd.list  | gzip > R_Unc.for_pl.fq.gz
$ seqtk subseq EPNC_orphan_ready_clean.fq.gz  Unclassified_from_rvbd.list | gzip   > O_Unc.for_pl.fq.gz

# Files ready for second round using rvdb..........

echo "1:Runing  kaiju-multi on paired end reads ..." 

$ kaiju-multi -z 4 -E 0.01 -t /kaiju_db/kaiju_db_plasmids_2022-04-10/nodes.dmp \
                         -f /kaiju_db/kaiju_db_plasmids_2022-04-10/kaiju_db_plasmids.fmi \
                        -i F_Unc_for_pl.fq.gz -j R_Unc.for_pl.fq.gz  > Metacongo_kaiju__PE_PL.out

# Runing  kaiju (not multi) on orphan merged reads ...



$ kaiju -z 4 -E 0.01 -t /kaiju_db/kaiju_db_plasmids_2022-04-10/nodes.dmp \
                     -f /kaiju_db/kaiju_db_plasmids_2022-04-10/kaiju_db_plasmids.fmi \
                    -i O_Unc.for_pl.fq.gz  > Metacongo_kaiju__SE_PL.out



# Combining the output for PE and ORPHAN (SE) here ....

$ cat Metacongo_kaiju__PE_PL.out  Metacongo_kaiju__SE_PL.out > Metacongo_kaiju__ALL_PL.out

# Adding full taxa names ... to output ...


$ kaiju-addTaxonNames -p  -t /kaiju_db/kaiju_db_plasmids_2022-04-10/nodes.dmp \
                          -n /kaiju_db/kaiju_db_plasmids_2022-04-10/names.dmp
                          -i Metacongo_kaiju__ALL_PL.out  -o Metacongo_kaiju__ALL_PL.out"_with_name.tsv"


# How many reads are classified .....

$ echo "TOTAL REDS  COUNT: " $(wc -l  $Metacongo_kaiju__ALL_PL.out |awk '{print $1}') 

$ echo "READ CLASSIFIED COUNT: " $(grep -w -c "C" $Metacongo_kaiju__ALL_PL.out ) 


# Converting to kaiju output to krona file

$ kaiju2krona  -t /kaiju_db/kaiju_db_plasmids_2022-04-10/nodes.dmp \
             -n  /home/databases/kaiju_db/kaiju_db_plasmids_2022-04-10/names.dmp \
             -i Metacongo_kaiju__ALL_PL.out -o Metacongo_kaiju__ALL_PL.krona



# Creating html from krona file...

$ ktImportText -o Metacongo_kaiju__ALL_PL.krona.html  Metacongo_kaiju__ALL_PL.krona


# Creating classification summary for phylum, class, order family, genus and species ...


$ for i in phylum class order family genus species; do kaiju2table -t /kaiju_db/kaiju_db_plasmids_2022-04-10/nodes.dmp \
                                                                  -n /kaiju_db/kaiju_db_plasmids_2022-04-10/names.dmp \
                                                                  -r $i -o Metacongo_kaiju__ALL_PL.out"_"$i"__summary.tsv" Metacongo_kaiju__ALL_PL.out; done





# All runs finished   ..... 

```

> From here all the profiling step was done using kaiju, will get un profiled reads fomr this step then convert them to fastq
> and profile them by another profiler with is Kraken2

```bash
# Extracting reads form PL outout

# Going to extract reads not assigned from PL profiling ...


$ grep -w 'U' Metacongo_kaiju__ALL_PL.out  |awk '{print $2}' > Unclassified_from_PL.list

echo "Extracting Unclassified reads from original fastq" |tee -a analysis.log



seqtk subseq EPNC_trim_ready_clean_R1.fq.gz  Unclassified_from_PL.list  | gzip > R1_to_kraken2.fq.gz
seqtk subseq EPNC_trim_ready_clean_R1.fq.gz  Unclassified_from_PL.list  | gzip > R2_to_kraken2.fq.gz
seqtk subseq EPNC_orphan_ready_clean.fq.gz  Unclassified_from_PL.list | gzip   > orphan_to_kraken2.fq.gz








```
#################################################################

## kraken2 Profiling

We will use these input for Kraken2

* kraken2_db_path='/home/databases/kraken2_bracken_ref_seq/standard_plus_PF/'
* R1='R1_to_kraken2.fq.gz'
* R2='R2_to_kraken2.fq.gz'
* orphan='orphan_to_kraken2.fq.gz'

#kraken2_report='metacongo_remaining_krake2.report'
#kraken2_output='metacongo_remaining_krake2.out'

#classified and unclassified reads
#uncla_reads='unclassified#.fq'
#cla_reads='classified#.fq'
#$SLURM_CPUS_PER_TASK


```bash
# Profiling the remaining reads from kaiju using kraken2 Using default parameters

# PAIRED
$ kraken2 --use-names --threads 10 --db /home/databases/kraken2_bracken_ref_seq/standard_plus_PF/ \
        --report metacongo_remaining_krake2.report  \
        --gzip-compressed --use-names \
        --paired  R1_to_kraken2.fq.gz R2_to_kraken2.fq.gz \
        --output metacongo_remaining_krake2.out \
        --unclassified-out unclassified#.fq --classified-out classified#.fq

# ORPHANS


$ kraken2  --db home/databases/kraken2_bracken_ref_seq/standard_plus_PF/ \
             orphan_to_kraken2.fq.gz  --use-names --report orphan_kraken.report \
            --output orphan_kraken2.output \
            --gzip-compressed --threads 10 \
            --classified_from_orphan.fq  unclassified_from_orphan.fq

# zipp fq files

gzip *fq

```

> Now, using kreport2krona.py from krakentools (to be added in REF)

```bash

# Paired
$ kreport2krona.py -r metacongo_remaining_kraken2.report --intermediate-ranks -o metacongo_remaining_kraken2.KRONA

# Single/orphans

$ kreport2krona.py -r orphan_kraken2.report --intermediate-ranks  -o orphan_kraken2.KRONA

#Generate html for visualization (you can use the script available with kaiju package ktImportText from KronaTools ...)

$ ktImportText -o metacongo_remaining_kraken2.KRONA.html  metacongo_remaining_kraken2.KRONA

$ ktImportText -o orphan_kraken2.KRONA.html orphan_kraken2.KRONA

```

# Assembly (another approach to analyse the data)

* Here we will assemble the data using two assembler, and bin assemblies using two binning tools. 
* Megahit and Metaspades will be used as assemblers
* Maxbin2 and Metabat2 will be used as binning tools


## Megahit assembly
```bash

$ mkdir -p ASSEMBLY/MEGAHIT_ASSEMBLY && cd ASSEMBLY/MEGAHIT_ASSEMBLY/

#Get files symlinked here ....

$ ln -s  ../../HUMAN_CONTA_REMOVAL/READY_FASTQ_FILES_CLEAN/EPNC_orphan_ready_clean.fq.gz
$ ln -s ../../HUMAN_CONTA_REMOVAL/READY_FASTQ_FILES_CLEAN/EPNC_trim_ready_clean_R1.fq.gz
$ ln -s ../../HUMAN_CONTA_REMOVAL/READY_FASTQ_FILES_CLEAN/EPNC_trim_ready_clean_R2.fq.gz

#Run megahit

$ megahit -1 EPNC_trim_ready_clean_R1.fq.gz  -2 EPNC_trim_ready_clean_R2.fq.gz -r EPNC_orphan_ready_clean.fq.gz \
           -o metacongo_final_assembly \
           -t 10  \
           --min-count 2 --k-min 21 --k-max 127 --k-step 2 --min-contig-len 200 --tmp-dir .


```
> Bonus (assembly statistics)


## Metaspades assembly 

```bash

$ mkdir METASPADES_ASSEMBLY  && cd METASPADES_ASSEMBLY
$ pwd 
/workspaces/Metacongo_Paper/ASSEMBLY/METASPADES_ASSEMBLY

# Symlink files as usual 

$ ln -s  ../../HUMAN_CONTA_REMOVAL/READY_FASTQ_FILES_CLEAN/EPNC_orphan_ready_clean.fq.gz
$ ln -s ../../HUMAN_CONTA_REMOVAL/READY_FASTQ_FILES_CLEAN/EPNC_trim_ready_clean_R1.fq.gz
$ ln -s ../../HUMAN_CONTA_REMOVAL/READY_FASTQ_FILES_CLEAN/EPNC_trim_ready_clean_R2.fq.gz


#Run metaspades
$ k_mer='21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,75,77,79,81,83,85,87,91,93,95,97,99,101,103,105,107,109,111,113,115,117,119,121,123,125,127'

metaspades.py --pe1-1 EPNC_trim_ready_clean_R1.fq.gz  --pe1-2 EPNC_trim_ready_clean_R2.fq.gz \
              -s  EPNC_orphan_ready_clean.fq.gz  -k $k_mer  -o metaspades_final_assembly_results


```
> Bonus (assembly statistics)

# Binning assemblies 

## Megahit assembly binning

```bash
$ mkdir MEGAHIT_ASSEM_BIN

# Symlink  fastq files and assembly

$ ln -s   ../../ASSEMBLY/MEGAHIT_ASSEMBLY/EPNC_orphan_ready_clean.fq.gz
$ ln -s   ../../ASSEMBLY/MEGAHIT_ASSEMBLY/EPNC_trim_ready_clean_R1.fq.gz
$ ln -s   ../../ASSEMBLY/MEGAHIT_ASSEMBLY/EPNC_trim_ready_clean_R2.fq.gz
$ ln -s   ../../ASSEMBLY/MEGAHIT_ASSEMBLY/metacongo_final_assembly/final.contigs.fa

# Symlink assembly
$ mv final.contigs.fa metacong_megahit_assembly.fsa



# Reference ... assembly file  Indexing ...............

$ bowtie2-build -f metacong_megahit_assembly.fsa  metacong_megahit_assembly.fsa__indexed --threads 10

# MAPP READS AGAINST THE INDEX

$ bowtie2 -p 10 --very-sensitive-local  -x metacong_megahit_assembly.fsa__indexed -1 EPNC_trim_ready_clean_R1.fq.gz  -2 EPNC_trim_ready_clean_R2.fq.gz \
                                        -U EPNC_orphan_ready_clean.fq.gz  |samtools view -bS - > metacongo.bam

# SORTING bam file
$ samtools sort metacongo.bam -o metacongo_sorted.bam

#Indexing ...
$ samtools index metacongo_sorted.bam

# removing unwanted file  
$ rm metacongo.bam

# Generating stats from the  mapping ..............

bamtools stats -in metacongo_sorted.bam >mapping.stats



# Start binning contigs using metabat2


$ runMetaBat.sh -m 1500 -t 10   metacong_megahit_assembly.fsa  metacongo_sorted.bam 


# Start binning using maxbin2 ...........


$ run_MaxBin.pl -contig metacong_megahit_assembly.fsa  -out maxbin2_results \
              -reads EPNC_trim_ready_clean_R1.fq.gz \
              -reads2 EPNC_trim_ready_clean_R2.fq.gz  \
              -reads3 EPNC_orphan_ready_clean.fq.gz \
              -min_contig_length 200 -thread 10 \
              -prob_threshold 0.7 -plotmarker -markerset 40 -verbose 

$ mkdir Maxbin2_binning_results  && mv maxbin* Maxbin2_binning_results

```

## Metaspades assembly binning

```bash

# Index the ref

$ bowtie2-build --large-index -f metacongo_metaspades_assembly.fsa  metacongo_metaspades_assembly.fsa__indexed --threads 10

# Map reads against assembly

$ bowtie2 -p 10 --very-sensitive-local   -x metacongo_metaspades_assembly.fsa__indexed -1 EPNC_trim_ready_clean_R1.fq.gz \
                                                                                     -2 EPNC_trim_ready_clean_R2.fq.gz \
                                                                                     -U EPNC_orphan_ready_clean.fq.gz |samtools view -bS - > metacongo.bam
#  Sort mapping file
$ samtools sort metacongo.bam -o metacongo_sorted.bam

# Index bam file
$ samtools index metacongo_sorted.bam

# Delete unwanted bam
rm metacongo.bam

# Get stats from bam file ....
$ bamtools stats -in metacongo_sorted.bam >mapping.stats

# start binning contigs using metabat2

$ runMetaBat.sh -m 1500 -t 10  $REF  metacongo_sorted.bam |tee -a  binning_analysis.log


$ mkdir Metabat2_binning_results && mv metacongo_metaspades_assembly.fsa.depth.txt \
        metacongo_metaspades_assembly.fsa.paired.txt *bt2l \
        metacongo_sorted.bam metacongo_sorted.bam.bai \
        mapping.stats metacongo_metaspades_assembly.fsa.metabat-bins20  Metabat2_binning_results




# Start binning using maxbin2 ...........

run_MaxBin.pl -contig metacongo_metaspades_assembly.fsa -out maxbin2_results -reads EPNC_trim_ready_clean_R1.fq.gz \
                                                                             -reads2 EPNC_trim_ready_clean_R2.fq.gz \
                                                                              -reads3 EPNC_orphan_ready_clean.fq.gz \
                                                                              -min_contig_length 200 -thread 10 \
                                                                              -prob_threshold 0.7 -plotmarker -markerset 40  

```




# Bins refinement 
For bins refinement we used metawrap (some modules).
The refinement was done on each assembly-bin

## Megahit Assembly Refined 

```bash
# Set the env for metawrap 
$ source ~/miniconda2/bin/activate
$ conda activate metawrap-env

$ metawrap  bin_refinement -o MEGAHIT_ASSEMB_REFINED -A 8.1.MEGAHIT_ASSEM_BIN/Maxbin2_binning_results -B 8.1.MEGAHIT_ASSEM_BIN/Metabat2_binning_results/ metacong_megahit_assembly.fsa.metabat-bins32


```


## Metaspades Assembly Refined

```bash
$ source ~/miniconda2/bin/activate
$ conda activate metawrap-env

$ metawrap  bin_refinement -o METASPADES_ASSEMB_REFINED/ -A 8.2.METASPADES_ASSEMB_BIN/Maxbin2_binning_results -B 8.2.METASPADES_ASSEMB_BIN/Metabat2_binning_results/metacongo_metaspades_assembly.fsa.metabat-bins20/

```

>> Refining refined == No results ... give an explanations



## Megahit Bins Refined renamed

form bin.XX.fa to mghit_bin.XX.fa


## Metaspades Bins Refined Renamed

form bin.XX.fa to mtspades_bin.XX.fa

## Renamed Refined Merged

cp the two folders (files of mghit_bin.XX.fa   && mtspades_bin.XX.fa) in one folder


# Bins Quantification

Metawrap quant_bins use only F and R or 1 and 2 files for reads, so we extrcated each pair from orphan and appended
to each file (1 or 2)

```
# Extract reads from orphan files ... and save them according to pairing info ..

cat EPNC_orphan_ready_clean.fastq|paste - - - -  |awk '$2~ /1:N/ {print $1,$2"\n"$3"\n"$4"\n"$5}' >1.fq
cat EPNC_orphan_ready_clean.fastq|paste - - - -  |awk '$2~ /2:N/ {print $1,$2"\n"$3"\n"$4"\n"$5}' >2.fq

# Concatenate files ...

cat EPNC_trim_ready_clean_1.fastq 1.fq> test_1.fastq 
cat EPNC_trim_ready_clean_2.fastq 2.fq> test_2.fastq

# Finally remove unwanted files ....... and rename 

rm 1.fq 2.fq EPNC_orphan_ready_clean.fastq EPNC_trim_ready_clean_1.fastq EPNC_trim_ready_clean_2.fastq

# Rename ........
mv test_1.fastq All_EPNC_trim_ready_clean_1.fastq
mv test_2.fastq All_EPNC_trim_ready_clean_2.fastq

```
>> Then we quantify bins by the module quant_bins from  metawrap pipeline ...


```bash

metawrap quant_bins  -b MERGED_REFINED_RENAMED -o QUANT_MERGED_REFINED_RENAMED All_EPNC_trim_ready_clean_1.fastq All_EPNC_trim_ready_clean_2.fastq


```




# BLOBOLOGY

metawrap blobology -a QUANT_MERGED_REFINED_RENAMED/assembly.fa -o BLOBOLOGY --bins MERGED_REFINED_RENAMED/ All_EPNC_trim_ready_clean_1.fastq All_EPNC_trim_ready_clean_2.fastq

# Bin Classification

#metawrap classify_bins -b MERGED_REFINED_RENAMED -o CLASSIFY_BIN -t 40

metaWRAP annotate_bins -o FUNCT_ANNOT -t 40 -b MERGED_REFINED_RENAMED

# Bins Taxonomy using GTDB Database

GENOME_DIR='../8.3.MTAWRAP_FROM_THEBEAST/8.3.5.MERGED_REINED_RENAMED/'
OUTDIR='GTDB_TAX'

EXTENSION='fa'


#Get conda ENV

source  /home/bioinf/apps/anaconda3/bin/activate

#Activate gtdb

conda activate gtdbtk

gtdbtk classify_wf --genome_dir $GENOME_DIR  --out_dir $OUTDIR --extension $EXTENSION   --cpus $SLURM_CPUS_PER_TASK --pplacer_cpus 1


## Clustering assembly 


> For this step we want to explore the presence of AMR genes but in all assembly

```bash

$ mkdir WHOLE_ASSEMBLY && cd WHOLE_ASSEMBLY

# Symlink assemblies (from megahit and metaspades ...)

$ ln -s ../../ASSEMBLY/MEGAHIT_ASSEMBLY/metacongo_final_assembly/final.contigs.fa
$ ln -s ../../ASSEMBLY/METASPADES_ASSEMBLY/metaspades_final_assembly_results/contigs.fasta

#Renamed assemblies file for good tracking

$ mv final.contigs.fa megahit_assembly.fa && mv contigs.fasta metaspades_assembly.fa

# Concatenate assembly files
cat megahit_assembly.fa metaspades_assembly.fa > all_assembly.fasta


$ mkdir Clustered_using_mmseq && cd Clustered_using_mmseq

$ ln -s ../all_assembly.fasta

# Run mmseq2 to clusterize the assembly (remove redundant contigs)
$ mmseqs easy-cluster all_assembly.fasta clusterRes tmp --min-seq-id 0.8 -c 0.8 --cov-mode 1
$ ll 
clusterRes_all_seqs.fasta
clusterRes_cluster.tsv
clusterRes_rep_seq.fasta (This file will be used to scan for AMR)

```

## AMRs detection

```bash
$ mkdir AMR_DETECT && cd AMR_DETECT
$ ln -s ../WHOLE_ASSEMBLY/clusterRes_rep_seq.fasta

# Run abricate
abricate clusterRes_rep_seq.fasta  --threads 10  > Abricate_ouput.tab

# Generate summary
abricate --summary  Abricate_ouput.tab> Abricate_ouput_summary.tab



```
 

# SUPPLEMENTARY COMMENTS



































