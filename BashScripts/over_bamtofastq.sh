#!/bin/bash
#SBATCH --job-name=bamtofastq		# Job name
#SBATCH --mail-type=ALL		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=aidan.burn@tufts.edu
#SBATCH --partition=preempt
#SBATCH --nodes=1			# Run a single task
#SBATCH --cpus-per-task=2		# Number of CPU cores per task
#SBATCH --mem-per-cpu=20g
#SBATCH --time=0-4:00:00
#SBATCH --output=/cluster/tufts/micr_coffin01/aburn01/Tools/Reports/%j.out
#SBATCH --error=/cluster/tufts/micr_coffin01/aburn01/Tools/Reports/%j.err
#SBATCH --constraint=rhel7

module load python/3.6.0
module load anaconda/3
module load samtools/1.9
module load picard/2.8.0
source /cluster/tufts/coffinlab/micr_coffin01/aburn01/GTEx_HML2/Start_Variables.txt

file=$1
sample_prefix=$2
out_dir=$3

echo "got the file $file"
echo "sample prefix is $sample_prefix"
echo "output directory is $out_dir"

python /cluster/tufts/micr_coffin01/aburn01/Tools/run_SamToFastq.py -p ${sample_prefix} -o ${out_dir}  ${file}

