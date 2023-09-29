#!/bin/bash
#SBATCH -J megahit_metacongo_fil_trim_assembly
#SBATCH -o megahit_metacongo_fil_trim_assembly_out_%j.txt
#SBATCH -e megahit_metacongo_fil_trim_assembly_err_%j.txt
#SBATCH --partition=large
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --workdir=.


#SBATCH --mail-type=ALL
#SBATCH --mail-user=bouziane.moumen@univ-poitiers.fr




module load bioinf

module load megahit/1.2.9

megahit -1 EPNC_trim_ready_clean_R1.fq.gz  -2 EPNC_trim_ready_clean_R2.fq.gz -r EPNC_orphan_ready_clean.fq.gz  -o metacongo_final_assembly -t $SLURM_CPUS_PER_TASK  --min-count 2 --k-min 21 --k-max 127 --k-step 2 --min-contig-len 200 --tmp-dir .




