#!/bin/bash
#SBATCH -J abricate
#SBATCH -o abricate_out_%j.txt
#SBATCH -e abricate_err_%j.txt
#SBATCH --partition=large
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --workdir=.




#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_mail_here

set -e
set -x
#PREPARE ENV
source /home/bioinf/apps/anaconda3/bin/activate
conda activate abricate


#Exit if error 


INPUT='clusterRes_rep_seq.fasta'
RESULTS='Abricate_ouput.tab'
SUMMARY='Abricate_ouput_summary.tab'

echo "Abricating ............" |tee -a abricate_analysis.log

echo "#######################################################" |tee -a abricate_analysis.log

abricate $INPUT  --threads $SLURM_CPUS_PER_TASK  > $RESULTS

echo "abricate finished .... we will generate summary ...." |tee -a abricate_analysis.log

abricate --summary $RESULTS > $SUMMARY




echo "generating summary finished ..........."

echo "#######################################################" |tee -a abricate_analysis.log





