#!/bin/bash
#SBATCH -J bowtie2_vs_human
#SBATCH -o bowtie2_vs_human_out_%j.txt
#SBATCH -e bowtie2_vs_human_err_%j.txt
#SBATCH --partition=large
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --workdir=.




#SBATCH --mail-type=ALL
#SBATCH --mail-user=bouziane.moumen@univ-poitiers.fr

# EXIT ON ERRORS
set -e
set -x
#PREPARE ENV

module load bioinf
module load seqtk
module load bowtie2/2.3.4.3
module load samtools/1.9
module load bamtools/2.5.1

source  /home/bioinf/apps/anaconda3/bin/activate
conda activate seqkit


# Get data for input (Note orphan reads from R1 and R2 are merged with simple cat command)

#This is the original names 
F_READS='EPNC_trim_R1.fq.gz'
R_READS='EPNC_trim_R2.fq.gz'
ORPHAN_READS='EPNC_orphan.fq.gz'

HUMAN_BW2_INDEX='6.1.REF/GRCh38_noalt_as'

OUT_SAM='Read_vs_human.sam'
OUT_BAM='Read_vs_human.bam'
OUT_UNMAPPED_BAM='Unmapped.bam'

OUT_BAM_SORTED='Unmapped_sorted.bam'
UNMAPPED_LIST='Unmapped.list'


echo  "Analysis starts at:" $(date) |tee -a bowtie2_mapping_analysis.log


bowtie2 --very-sensitive-local -x $HUMAN_BW2_INDEX -1 $F_READS  -2 $R_READS -U $ORPHAN_READS -S $OUT_SAM -p $SLURM_CPUS_PER_TASK 



echo "...........Mapping done ............."|tee -a bowtie2_mapping_analysis.log

echo "Converting sam to bam " |tee -a bowtie2_mapping_analysis.log

samtools view -S -bh $OUT_SAM > $OUT_BAM

echo "Sort bam file .................." |tee -a bowtie2_mapping_analysis.log

samtools sort $OUT_BAM -o $OUT_BAM_SORTED


echo "Getting some stats from the whole BAM file ...." |tee -a bowtie2_mapping_analysis.log

bamtools stats -in $OUT_BAM_SORTED >mapping_stat_from_whole_bam.stats


echo "Getting unmapped reads from the bam file ...."|tee -a bowtie2_mapping_analysis.log

samtools view -b -f 4 $OUT_BAM_SORTED > $OUT_UNMAPPED_BAM

echo "Get a list from the unmapped file " |tee -a bowtie2_mapping_analysis.log


# This line produce the unmapped list BUT with duplication
# Have to remove duplicated reads at the end using seqkit 
samtools view $OUT_UNMAPPED_BAM  |awk '{print $1}' >$UNMAPPED_LIST


echo "Using seqtk to get unmapped reads in fastq format ..."|tee -a bowtie2_mapping_analysis.log


seqtk subseq $F_READS  $UNMAPPED_LIST  | gzip > EPNC_trim_no_human_R1.fq.gz
seqtk subseq $R_READS  $UNMAPPED_LIST  | gzip > EPNC_trim_no_human_R2.fq.gz
seqtk subseq $ORPHAN_READS   $UNMAPPED_LIST | gzip   > EPNC_orphan_no_human.fq.gz





echo "CLEANING ...................................." |tee -a bowtie2_mapping_analysis.log


rm $OUT_BAM $OUT_SAM $OUT_BAM_SORTED $OUT_UNMAPPED_BAM 




echo  "	ACTIVATING CONDA ENV FOR SEQKIT TO DO THE LAST JOB = REMOVING REMAINING HUMAN READS " |tee -a bowtie2_mapping_analysis.log




echo "Removing last remaining human reads from previous kraken run ..." |tee -a bowtie2_mapping_analysis.log

seqkit grep -f human_from_kraken2.list -v EPNC_trim_no_human_R1.fq.gz -o EPNC_trim_ready_R1.fq.gz

seqkit grep -f human_from_kraken2.list -v EPNC_trim_no_human_R2.fq.gz -o EPNC_trim_ready_R2.fq.gz

seqkit grep -f human_from_kraken2.list -v EPNC_orphan_no_human.fq.gz -o EPNC_orphan_ready.fq.gz


mkdir 6.2.READY_FASTQ_FILES

mv EPNC_trim_ready_R1.fq.gz EPNC_trim_ready_R2.fq.gz EPNC_orphan_ready.fq.gz 6.2.READY_FASTQ_FILES




echo  "Analysis terminated without error  at:" $(date) |tee -a bowtie2_mapping_analysis.log


#NB: DO NOT  FORGET TO REMOVE DUPLICATED READS FROM THESES FILES/
	#EPNC_trim_ready_R1.fq.gz
	#EPNC_trim_ready_R2.fq.gz
	#EPNC_orphan_ready.fq.gz
