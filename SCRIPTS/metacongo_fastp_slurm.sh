#!/bin/bash
#SBATCH -J metacongo_fastp
#SBATCH -o metacongo_fastp_out_%j.txt
#SBATCH -e metacongo_fastp_err_%j.txt
#SBATCH --partition=large
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --workdir=.


#SBATCH --mail-type=ALL
#SBATCH --mail-user=bouziane.moumen@univ-poitiers.fr

#LOAD WANTED MODULES
module load bioinf
module load Fastp/0.21.0

##Get some Variables

R1='EPNC_R1.fq.gz'
R1_TRIMMED='EPNC_trim_R1.fq.gz'
R1_ORPHAN='EPNC_orphan_1.fq.gz'


R2='EPNC_R2.fq.gz'
R2_TRIMMED='EPNC_trim_R2.fq.gz'
R2_ORPHAN='EPNC_orphan_2.fq.gz'


MINIMUM_LENGTH='50'

REPORT_FILE='metacongo_fastp_report.html'
JSON_REPORT_FILES='metacongo_fastp_report.json'




fastp -i $R1 -o $R1_TRIMMED -I $R2 -O $R2_TRIMMED --unpaired1 $R1_ORPHAN --unpaired2 $R2_ORPHAN  -z 4 --trim_poly_g --trim_poly_x  --detect_adapter_for_pe -l $MINIMUM_LENGTH -c -p -h  $REPORT_FILE --json $JSON_REPORT_FILES --overrepresentation_analysis  -w $SLURM_CPUS_PER_TASK   --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA  --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT




