#!/bin/bash
#SBATCH --job-name=Telescope		# Job name
#SBATCH --mail-type=END,FAIL		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=aidan.burn@tufts.edu
#SBATCH --partition=preempt	
#SBATCH --nodes=1			# Run a single task
#SBATCH --cpus-per-task=4		# Number of CPU cores per task
#SBATCH --mem-per-cpu=32000 
#SBATCH --time=7:00:00
#SBATCH --output=/cluster/tufts/micr_coffin01/aburn01/Tools/Reports/%j.out
#SBATCH --error=/cluster/tufts/micr_coffin01/aburn01/Tools/Reports/%j.err
#SBATCH --constraint=rhel7
#SBATCH --array=0-7

module load anaconda/3
module load multiqc
source activate /cluster/tufts/micr_coffin01/aburn01/condaenv/telescope_env
source /cluster/tufts/micr_coffin01/aburn01/GTEx_HML2/Variables.txt

cd ${1}
FILE=($(ls  *_HISAT2_hg38_align.sam | rev | cut -c 23- | rev | uniq))

FILENAME=${FILE[$SLURM_ARRAY_TASK_ID]}
echo $FILENAME

telescope assign --attribute gene_id --outdir ${1} --exp_tag ${FILENAME}  ${FILENAME}_HISAT2_hg38_align.sam ${Ref_Dir}/hg38_HML2_genes_051918.gtf

source deactivate
