#!/bin/bash
#SBATCH -J metacongo_megahit_binning
#SBATCH -o metacongo_megahit_binning_out_%j.txt
#SBATCH -e metacongo_megahit_binning_err_%j.txt
#SBATCH --partition=large
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --workdir=.


#SBATCH --mail-type=ALL
#SBATCH --mail-user=bouziane.moumen@univ-poitiers.fr


set -e


module load bioinf
module load bowtie2/2.3.4.3
module load samtools/1.9
module load bamtools/2.5.1
module load METABAT_2/metabat_2



#DEFFINE VAR


REF='metacong_megahit_assembly.fsa'
INDEX=$REF"__indexed"
R1='EPNC_trim_ready_clean_R1.fq.gz'
R2='EPNC_trim_ready_clean_R2.fq.gz'
ORPHAN='EPNC_orphan_ready_clean.fq.gz'
MAPPING_MODE='--very-sensitive-local'


#INDEX THE REF

echo "Reference ___assembly__file__ Indexing ..............." |tee -a megahit_binning_analysis.log

bowtie2-build -f $REF  $INDEX --threads $SLURM_CPUS_PER_TASK

echo "Indexing done ..............."


#MAPP READS AGAINST THE INDEX

echo "Mapping reads to indexed assembly bw2 index,.... be patient" |tee -a megahit_binning_analysis.log

bowtie2 -p $SLURM_CPUS_PER_TASK $MAPPING_MODE  -x $INDEX -1 $R1 -2 $R2 -U $ORPHAN  |samtools view -bS - > metacongo.bam

echo "Mapping step done .............." | tee -a megahit_binning_analysis.log

#SORTING bam file

echo "Sorting and indexing  bam file ............" | tee -a megahit_binning_analysis.log
samtools sort metacongo.bam -o metacongo_sorted.bam

samtools index metacongo_sorted.bam
echo "Sorting bam file done ....." | tee -a megahit_binning_analysis.log

rm metacongo.bam


echo "Generating stats from the  mapping .............." | tee -a megahit_binning_analysis.log
bamtools stats -in metacongo_sorted.bam >mapping.stats

echo "MAPPING COMPLETLY DONE ................." | tee -a megahit_binning_analysis.log



echo "############################################################################################"


echo "STRAT BINNING CONTIGS USING METABAT2" | tee -a megahit_binning_analysis.log

runMetaBat.sh -m 1500 -t $SLURM_CPUS_PER_TASK  $REF  metacongo_sorted.bam |tee -a  megahit_binning_analysis.log


echo "Metabat2 binning done .... moving results to new folder"

mkdir Metabat2_binning_results && mv metacong_megahit_assembly.fsa.depth.txt metacong_megahit_assembly.fsa.paired.txt metacong_megahit_assembly.fsa.metabat-bins32 *bt2  metacongo_sorted.bam metacongo_sorted.bam.bai mapping.stats  Metabat2_binning_results

echo "Moving  results done .........." | tee -a megahit_binning_analysis.log

echo "############################################################################################"

echo " Binning using maxbin2 ............" | tee -a megahit_binning_analysis.log


# Activate conda and conda env" 

echo "Activating env for maxbin2"  | tee -a megahit_binning_analysis.log

source /home/bioinf/apps/anaconda3/bin/activate 


conda activate maxbin2


run_MaxBin.pl -contig $REF -out maxbin2_results -reads $R1  -reads2 $R2   -reads3 $ORPHAN  -min_contig_length 200 -thread $SLURM_CPUS_PER_TASK -prob_threshold 0.7 -plotmarker -markerset 40 -verbose |tee -a megahit_binning_analysis.log


mkdir Maxbin2_binning_results  && mv maxbin* Maxbin2_binning_results



conda deactivate 
conda deactivate


echo "maxbin binning done ..............." | tee -a megahit_binning_analysis.log







