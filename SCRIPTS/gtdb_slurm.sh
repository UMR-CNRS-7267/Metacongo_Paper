#!/bin/bash
#SBATCH -J gtdb__tax
#SBATCH -o gtdb_out_%j.txt
#SBATCH -e gtdb_err_%j.txt
#SBATCH --partition=large
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --workdir=.
#SBATCH --nodelist=node03
#SBATCH --mem-per-cpu=380000MB

#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_mail_here

#LOAD WANTED MODULES

##Get some Variables




#VARIABLE=''

GENOME_DIR='/MERGED_REINED_RENAMED/'
OUTDIR='GTDB_TAX'

EXTENSION='fa'


#Get conda ENV

source  /home/bioinf/apps/anaconda3/bin/activate

#Activate gtdb

conda activate gtdbtk


echo "STARTING at ......................."

echo $(date)


gtdbtk classify_wf --genome_dir $GENOME_DIR  --out_dir $OUTDIR --extension $EXTENSION   --cpus $SLURM_CPUS_PER_TASK --pplacer_cpus 1



echo "gtdb finished at ...................."

echo $(date)







