#!/bin/bash
#SBATCH -J metaspades_metacongo_fil_trim_assembly
#SBATCH -o metaspades_metacongo_fil_trim_assembly_out_%j.txt
#SBATCH -e metaspades_metacongo_fil_trim_assembly_err_%j.txt
#SBATCH --partition=large
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --workdir=.


#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_mail_here




module load bioinf

module load SPAdes/3.15.2

ASSEMBLY_OUTPUT='metaspades_final_assembly_results'

F_READS='EPNC_trim_ready_clean_R1.fq.gz'
R_READS='EPNC_trim_ready_clean_R2.fq.gz'
ORPHAN_READS='EPNC_orphan_ready_clean.fq.gz'


k_mer='21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,75,77,79,81,83,85,87,91,93,95,97,99,101,103,105,107,109,111,113,115,117,119,121,123,125,127'





metaspades.py --pe1-1 $F_READS  --pe1-2 $R_READS  -s  $ORPHAN_READS  -k $k_mer  -o $ASSEMBLY_OUTPUT


