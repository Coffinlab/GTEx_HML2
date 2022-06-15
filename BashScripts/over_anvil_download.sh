#!/bin/bash
#SBATCH --job-name=Download              # Job name
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=aidan.burn@tufts.edu
#SBATCH --partition=batch
#SBATCH --nodes=1                       # Run a single task
#SBATCH --cpus-per-task=1               # Number of CPU cores per task
#SBATCH --mem-per-cpu=16000
#SBATCH --time=-2:00:00
#SBATCH --output=/cluster/tufts/micr_coffin01/aburn01/Tools/Reports/%j.out
#SBATCH --error=/cluster/tufts/micr_coffin01/aburn01/Tools/Reports/%j.err

source /cluster/tufts/micr_coffin01/aburn01/GTEx_HML2/Start_Variables.txt
module load python/3.6.0
cd /cluster/tufts/micr_coffin01/aburn01/Tools

sample=$(echo ${1}.tsv | rev | cut -c 5- | rev | uniq)
echo ${sample}


cut -f 2 $sample.tsv > ${sample}_cut.tsv
tail -n +2 ${sample}_cut.tsv > ${sample}_tail.tsv
head ${sample}_tail.tsv
rm ${sample}.tsv ${sample}_cut.tsv

mkdir /cluster/tufts/coffinlab/gtex/aburn01/V8/${sample}
cat ${sample}_tail.tsv | gsutil -u gtex-terra-project cp -I  /cluster/tufts/coffinlab/gtex/aburn01/V8/${sample}
