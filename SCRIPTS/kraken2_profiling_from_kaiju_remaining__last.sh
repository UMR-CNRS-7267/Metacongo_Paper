#!/bin/bash
#SBATCH -J kraken2_profiling_from_kaiju_remaining__last
#SBATCH -o kraken2_profiling_from_kaiju_remaining__last_out_%j.txt
#SBATCH -e kraken2_profiling_from_kaiju_remaining__last_err_%j.txt
#SBATCH --partition=large
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --workdir=.




#SBATCH --mail-type=ALL
#SBATCH --mail-user=bouziane.moumen@univ-poitiers.fr



#PREPARE ENV
source /home/bioinf/apps/anaconda3/bin/activate
conda activate kraken2
module load bioinf
module load seqtk



##The reads count  which are not assigned from kaiju 

##([bmoumen@volcano 9.2.USING_KRAKEN2]$ for i in orphan_to_kraken2.fq.gz R1_to_kraken2.fq.gz R2_to_kraken2.fq.gz ; do printf $i"\t"; gzip -cd $i |grep -c "@GWNJ"; done |awk '{sum+=$2} END {print sum}'
## 93192655
## [bmoumen@volcano 9.2.USING_KRAKEN2]$ for i in *gz; do printf $i"\t"; gzip -cd $i |grep -c "@GWNJ"; done
## orphan_to_kraken2.fq.gz 3017613
## R1_to_kraken2.fq.gz     45087521
## R2_to_kraken2.fq.gz     45087521




#VAR
kraken2_db_path='/home/databases/kraken2_bracken_ref_seq/standard_plus_PF/'
R1='R1_to_kraken2.fq.gz'
R2='R2_to_kraken2.fq.gz'
orphan='orphan_to_kraken2.fq.gz'

kraken2_report='metacongo_remaining_krake2.report'
kraken2_output='metacongo_remaining_krake2.out'

#classified and unclassified reads
uncla_reads='unclassified#.fq'
cla_reads='classified#.fq'


$SLURM_CPUS_PER_TASK

echo " Profiling the remaining reads from kaiju using kraken2" |tee -a analysis.log


echo "Using default parameters ..............." |tee -a analysis.log


kraken2 --use-names --threads $SLURM_CPUS_PER_TASK  --db $kraken2_db_path --report $kraken2_report  --gzip-compressed  --paired  $R1 $R2 --output $kraken2_output --unclassified-out $uncla_reads --classified-out $cla_reads



echo "kraken 2 profiling done using ...." |tee -a analysis.log





echo "will profile the remaining orphan files alone ...." |tee -a analysis.log



kraken2  --db $kraken2_db_path  orphan_to_kraken2.fq.gz  --use-names --report orphan_kraken.report  --output orphan_kraken.output --gzip-compressed --threads $SLURM_CPUS_PER_TASK --classified_from_orphan.fq  unclassified_from_orphan.fq


#zipp fq files

echo "we will zip files in fq ....uncl"

gzip *fq



echo "Analysis terminated ...." |tee -a analysis.log




#using confidence 0.9
#kraken2 --use-names --threads 60 --confidence 0.9 --db ../CONGO_ASCEl_TO_BE_COPIED_AFTER/Assembled_contigs_on_volcano/KRAKEN2_C
#USTOM/  --report stemo_illumina_c0.9_report --gzip-compressed --paired stenomatis_R1.fq.gz stenomatis_R2.fq.gz --output stemo_i
#llumina_kraken_c0.9.output --unclassified-out stenomatis_uncla#.fq --classified-out stenomatis_cla#.fq

 
