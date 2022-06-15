#!/bin/bash
#SBATCH --job-name=Telescope		# Job name
#SBATCH --mail-type=END,FAIL		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=aidan.burn@tufts.edu
#SBATCH --partition=batch
#SBATCH --nodes=1			# Run a single task
#SBATCH --cpus-per-task=1		# Number of CPU cores per task
#SBATCH --mem-per-cpu=32000
#SBATCH --time=1:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --constraint=rhel7
#SBATCH --array=0-30%5

module load anaconda/3
module load multiqc
source activate /cluster/tufts/micr_coffin01/aburn01/condaenv/telescope_env
source /cluster/tufts/micr_coffin01/aburn01/GTEx_HML2/Variables.txt

cd ${1}
multiqc .
mkdir ${1}_report
mv *report.tsv ${1}_report
mv multiqc_report.html ${1}_report
