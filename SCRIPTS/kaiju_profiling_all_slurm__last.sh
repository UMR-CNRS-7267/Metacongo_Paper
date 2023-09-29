#!/bin/bash
#SBATCH -J kaiju_bench
#SBATCH -o kaiju_bench_out_%j.txt
#SBATCH -e kaiju_bench_err_%j.txt
#SBATCH --partition=large
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --workdir=.




#SBATCH --mail-type=ALL
#SBATCH --mail-user=bouziane.moumen@univ-poitiers.fr

set -e
set -x
#PREPARE ENV
source /home/bioinf/apps/anaconda3/bin/activate
conda activate kaiju

module load bioinf
module load seqtk

#Exit if error 


#Get databases for uses

#KAIJU_NR_DB='/home/databases/kaiju_db/kaiju_db_nr_2022-03-10/kaiju_db_nr.fmi'
#KAIJU_NR_NODES='/home/databases/kaiju_db/kaiju_db_nr_2022-03-10/nodes.dmp'
#KAIJU_NR_NAMES='/home/databases/kaiju_db/kaiju_db_nr_2022-03-10/names.dmp'

KAIJU_NR_EUK_DB='/home/databases/kaiju_db/kaiju_db_nr_euk_2022-03-10/kaiju_db_nr_euk.fmi'
KAIJU_NR_EUK_NODES='/home/databases/kaiju_db/kaiju_db_nr_euk_2022-03-10/nodes.dmp'
KAIJU_NR_EUK_NAMES='/home/databases/kaiju_db/kaiju_db_nr_euk_2022-03-10/names.dmp'

KAIJU_PL_DB='/home/databases/kaiju_db/kaiju_db_plasmids_2022-04-10/kaiju_db_plasmids.fmi'
KAIJU_PL_NODES='/home/databases/kaiju_db/kaiju_db_plasmids_2022-04-10/nodes.dmp'
KAIJU_PL_NAMES='/home/databases/kaiju_db/kaiju_db_plasmids_2022-04-10/names.dmp'

#KAIJU_REFSEQ_DB='/home/databases/kaiju_db/kaiju_db_refseq_2022-03-23/kaiju_db_refseq.fmi'
#KAIJU_REFSEQ_NODES='/home/databases/kaiju_db/kaiju_db_refseq_2022-03-23/nodes.dmp'
#KAIJU_REFSEQ_NAMES='/home/databases/kaiju_db/kaiju_db_refseq_2022-03-23/names.dmp'

KAIJU_RVDB_DB='/home/databases/kaiju_db/kaiju_db_rvdb_2022-04-07/kaiju_db_rvdb.fmi'
KAIJU_RVDB_NODES='/home/databases/kaiju_db/kaiju_db_rvdb_2022-04-07/nodes.dmp'
KAIJU_RVDB_NAMES='/home/databases/kaiju_db/kaiju_db_rvdb_2022-04-07/names.dmp'

# Get data for input (Note orphan reads from R1 and R2 are merged with simple cat command)

#This is the original names 
F_READS='EPNC_trim_ready_clean_R1.fq.gz'
R_READS='EPNC_trim_ready_clean_R2.fq.gz'
ORPHAN_READS='EPNC_orphan_ready_clean.fq.gz'

#Output for the profiling using NR_EUK
OUTPUT_PE_NR_EUK='Metacongo_kaiju__PE_NR_EUK.out'
OUTPUT_SE_NR_EUK='Metacongo_kaiju__SE_NR_EUK.out'
OUTPUT_ALL_NR_EUK='Metacongo_kaiju__ALL_NR_EUK.out'

#After merging the output of PE and SE have to create krona file
OUTPUT_ALL_NR_EUK_KR='Metacongo_kaiju__ALL_NR_EUK.krona'


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











#OUTPUT_NR='Metacongo_kaiju__NR.out'
#OUTPUT_NR_KR='Metacongo_kaiju__NR.krona'





#OUTPUT_PL='Metacongo_kaiju__PL.out'
#OUTPUT_PL_KR='Metacongo_kaiju__PL.krona'

#OUTPUT_REFSEQ='Metacongo_kaiju__REFSEQ.out'
#OUTPUT_REFSEQ_KR='Metacongo_kaiju__REFSEQ.krona'

#OUTPUT_RVDB='Metacongo_kaiju__RVDB.out'
#OUTPUT_RVDB_KR='Metacongo_kaiju__RVDB.krona'



echo "###################STRTING PROFILING OF YOUR DATA#################"  |tee -a analysis.log


echo  "Analysis starts at:" $(date) |tee -a analysis.log


echo "Using NR_EUK DATABASES" |tee -a analysis.log



echo "1:Runing  kaiju-multi on paired end reads ..." |tee -a analysis.log

kaiju-multi -z $SLURM_CPUS_PER_TASK -E 0.01 -t $KAIJU_NR_EUK_NODES  -f $KAIJU_NR_EUK_DB -i $F_READS -j $R_READS  > $OUTPUT_PE_NR_EUK

echo "2:Runing  kaiju (not multi) on orphan merged reads ..." |tee -a analysis.log

kaiju -z $SLURM_CPUS_PER_TASK -E 0.01 -t $KAIJU_NR_EUK_NODES  -f $KAIJU_NR_EUK_DB -i $ORPHAN_READS  > $OUTPUT_SE_NR_EUK


echo "Combining the output for PE and ORPHAN (SE) here" |tee -a analysis.log

cat $OUTPUT_PE_NR_EUK $OUTPUT_SE_NR_EUK > $OUTPUT_ALL_NR_EUK  


echo "Adding full taxa names ... to output" |tee -a analysis.log

kaiju-addTaxonNames -p  -t $KAIJU_NR_EUK_NODES -n $KAIJU_NR_EUK_NAMES  -i $OUTPUT_ALL_NR_EUK  -o $OUTPUT_ALL_NR_EUK"_with_name.tsv"


echo "How many reads are classified" |tee -a analysis.log

echo "TOTAL REDS  COUNT: " $(wc -l  $OUTPUT_ALL_NR_EUK |awk '{print $1}') |tee -a analysis.log

echo "READ CLASSIFIED COUNT: " $(grep -w -c "C" $OUTPUT_ALL_NR_EUK ) |tee -a analysis.log


echo "converting to kaiju output to krona file" |tee -a analysis.log

kaiju2krona  -t $KAIJU_NR_EUK_NODES -n  $KAIJU_NR_EUK_NAMES -i $OUTPUT_ALL_NR_EUK -o $OUTPUT_ALL_NR_EUK_KR



echo "creating html from krona file" |tee -a analysis.log

ktImportText -o $OUTPUT_ALL_NR_EUK_KR".html"  $OUTPUT_ALL_NR_EUK_KR


echo "creating classification summary for phylum, class, order family, genus and species" |tee -a analysis.log

#kaiju2table -t $KAIJU_NR_NODES -n $KAIJU_NR_NAMES -r genus -o $OUTPUT_NR"__summary" $OUTPUT_NR
for i in phylum class order family genus species; do kaiju2table -t $KAIJU_NR_EUK_NODES -n $KAIJU_NR_EUK_NAMES -r $i -o $OUTPUT_ALL_NR_EUK"_"$i"__summary.tsv" $OUTPUT_ALL_NR_EUK; done




echo "Analysis sing NR_EUK DATABASES finished " |tee -a analysis.log

echo "#######################################################" |tee -a analysis.log

echo "Using RVDB DATABASE" |tee -a analysis.log


echo "Going to extract the unclassified reads from the output and re_run on Virus db" |tee -a  analysis.log


echo " getting reads list" |tee -a analysis.log

grep -w 'U' $OUTPUT_ALL_NR_EUK |awk '{print $2}' > Unclassified_from_nr_euk.list

echo "Extracting Unclassified reads from original fastq" |tee -a analysis.log



seqtk subseq  $F_READS  Unclassified_from_nr_euk.list  | gzip > $F_Unc_from_nr_euk
seqtk subseq  $R_READS  Unclassified_from_nr_euk.list  | gzip > $R_Unc_from_nr_euk
seqtk subseq  $ORPHAN_READS  Unclassified_from_nr_euk.list | gzip > $O_Unc_from_nr_euk

echo "..............Files ready for second round using rvdb..........."



echo "1:Runing  kaiju-multi on paired end reads ..." |tee -a analysis.log

kaiju-multi -z $SLURM_CPUS_PER_TASK -E 0.01 -t $KAIJU_RVDB_NODES  -f $KAIJU_RVDB_DB -i $F_Unc_from_nr_euk -j $R_Unc_from_nr_euk  > $OUTPUT_PE_RVDB

echo "2:Runing  kaiju (not multi) on orphan merged reads ..." |tee -a analysis.log

kaiju -z $SLURM_CPUS_PER_TASK -E 0.01 -t $KAIJU_RVDB_NODES  -f $KAIJU_RVDB_DB -i $O_Unc_from_nr_euk  > $OUTPUT_SE_RVDB


echo "Combining the output for PE and ORPHAN (SE) here" |tee -a analysis.log

cat $OUTPUT_PE_RVDB $OUTPUT_SE_RVDB > $OUTPUT_ALL_RVDB 


echo "Adding full taxa names ... to output" |tee -a analysis.log

kaiju-addTaxonNames -p  -t $KAIJU_RVDB_NODES -n $KAIJU_RVDB_NAMES  -i $OUTPUT_ALL_RVDB  -o $OUTPUT_ALL_RVDB"_with_name.tsv"


echo "How many reads are classified" |tee -a analysis.log

echo "TOTAL REDS  COUNT: " $(wc -l  $OUTPUT_ALL_RVDB |awk '{print $1}') |tee -a analysis.log

echo "READ CLASSIFIED COUNT: " $(grep -w -c "C" $OUTPUT_ALL_RVDB ) |tee -a analysis.log


echo "converting to kaiju output to krona file" |tee -a analysis.log

kaiju2krona  -t $KAIJU_RVDB_NODES -n  $KAIJU_RVDB_NAMES -i $OUTPUT_ALL_RVDB -o $OUTPUT_ALL_RVDB_KR



echo "creating html from krona file" |tee -a analysis.log

ktImportText -o $OUTPUT_ALL_RVDB_KR".html"  $OUTPUT_ALL_RVDB_KR


echo "creating classification summary for phylum, class, order family, genus and species" |tee -a analysis.log

#kaiju2table -t $KAIJU_NR_NODES -n $KAIJU_NR_NAMES -r genus -o $OUTPUT_NR"__summary" $OUTPUT_NR
for i in phylum class order family genus species; do kaiju2table -t $KAIJU_RVDB_NODES -n $KAIJU_RVDB_NAMES -r $i -o $OUTPUT_ALL_RVDB"_"$i"__summary.tsv" $OUTPUT_ALL_RVDB; done



echo "Analysis finished for RVDB profiling  at : " echo$(date)




echo "#######################################################" |tee -a analysis.log



echo "Using PL (for plasmids) DATABASE" |tee -a analysis.log


echo "Going to extract the unclassified reads from  previous analysis output and re_run on Plasmid db" |tee -a  analysis.log


echo " getting reads list" |tee -a analysis.log



grep -w 'U' $OUTPUT_ALL_RVDB |awk '{print $2}' > Unclassified_from_rvbd.list

echo "Extracting Unclassified reads from original fastq" |tee -a analysis.log



seqtk subseq $F_READS  Unclassified_from_rvbd.list  | gzip > $F_Unc_from_rvbd
seqtk subseq $R_READS   Unclassified_from_rvbd.list  | gzip > $R_Unc_from_rvbd
seqtk subseq $ORPHAN_READS  Unclassified_from_rvbd.list | gzip   > $O_Unc_from_rvbd

echo "..............Files ready for second round using rvdb..........."



echo "1:Runing  kaiju-multi on paired end reads ..." |tee -a analysis.log

kaiju-multi -z $SLURM_CPUS_PER_TASK -E 0.01 -t $KAIJU_PL_NODES  -f $KAIJU_PL_DB -i $F_Unc_from_rvbd -j $R_Unc_from_rvbd  > $OUTPUT_PE_PL

echo "2:Runing  kaiju (not multi) on orphan merged reads ..." |tee -a analysis.log

kaiju -z $SLURM_CPUS_PER_TASK -E 0.01 -t $KAIJU_PL_NODES  -f $KAIJU_PL_DB -i $O_Unc_from_rvbd  > $OUTPUT_SE_PL


echo "Combining the output for PE and ORPHAN (SE) here" |tee -a analysis.log

cat $OUTPUT_PE_PL $OUTPUT_SE_PL > $OUTPUT_ALL_PL 


echo "Adding full taxa names ... to output" |tee -a analysis.log

kaiju-addTaxonNames -p  -t $KAIJU_PL_NODES -n $KAIJU_PL_NAMES  -i $OUTPUT_ALL_PL  -o $OUTPUT_ALL_PL"_with_name.tsv"


echo "How many reads are classified" |tee -a analysis.log

echo "TOTAL REDS  COUNT: " $(wc -l  $OUTPUT_ALL_PL |awk '{print $1}') |tee -a analysis.log

echo "READ CLASSIFIED COUNT: " $(grep -w -c "C" $OUTPUT_ALL_PL ) |tee -a analysis.log


echo "converting to kaiju output to krona file" |tee -a analysis.log

kaiju2krona  -t $KAIJU_PL_NODES -n  $KAIJU_PL_NAMES -i $OUTPUT_ALL_PL -o $OUTPUT_ALL_PL_KR



echo "creating html from krona file" |tee -a analysis.log

ktImportText -o $OUTPUT_ALL_PL_KR".html"  $OUTPUT_ALL_PL_KR


echo "creating classification summary for phylum, class, order family, genus and species" |tee -a analysis.log

#kaiju2table -t $KAIJU_NR_NODES -n $KAIJU_NR_NAMES -r genus -o $OUTPUT_NR"__summary" $OUTPUT_NR
for i in phylum class order family genus species; do kaiju2table -t $KAIJU_PL_NODES -n $KAIJU_PL_NAMES -r $i -o $OUTPUT_ALL_PL"_"$i"__summary.tsv" $OUTPUT_ALL_PL; done



echo "Analysis finished for PL profiling  at : " echo$(date)  |tee -a analysis.log



echo "############################################################################"

echo "all runs finsihed ate :" echo $(date) |tee -a analysis.log



#kaiju2krona and kaiju2table need nodes and 


#Extracting reads form PL outout

echo "Going to extract reads not assigned from PL profiling ..."  |tee -a analysis.log


grep -w 'U' $OUTPUT_ALL_PL  |awk '{print $2}' > Unclassified_from_PL.list

echo "Extracting Unclassified reads from original fastq" |tee -a analysis.log



seqtk subseq $F_READS  Unclassified_from_PL.list  | gzip > R1_to_kraken2.fq.gz
seqtk subseq $R_READS  Unclassified_from_PL.list  | gzip > R2_to_kraken2.fq.gz
seqtk subseq $ORPHAN_READS  Unclassified_from_PL.list | gzip   > orphan_to_kraken2.fq.gz



echo "ALL PROFLING DONE ................" |tee -a analysis.log


echo "#################################################################"




