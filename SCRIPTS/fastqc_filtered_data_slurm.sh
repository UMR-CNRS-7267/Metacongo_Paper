#!/bin/bash
#SBATCH -J fastqc_filtered_data
#SBATCH -o fastqc_fastqc_filtered_data_out_%j.txt
#SBATCH -e fastqc_fastqc_filtered_data_err_%j.txt
#SBATCH --partition=large
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --workdir=.


#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_mail_here

##Get some Variables

module load bioinf

module load FastQC/0.11.8


for i in *gz; do fastqc  -t $SLURM_CPUS_PER_TASK -q $i; done 








