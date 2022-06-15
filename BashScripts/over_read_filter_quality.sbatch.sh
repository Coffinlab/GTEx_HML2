#!/bin/bash
#SBATCH --job-name=Trim_QC		# Job name
#SBATCH --mail-type=END,FAIL		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=aidan.burn@tufts.edu
#SBATCH --partition=preempt
#SBATCH --nodes=1			# Run a single task
#SBATCH --cpus-per-task=2		# Number of CPU cores per task
#SBATCH --mem-per-cpu=12000
#SBATCH --time=12:00:00
#SBATCH --output=/cluster/tufts/micr_coffin01/aburn01/Tools/Reports/%j.out
#SBATCH --error=/cluster/tufts/micr_coffin01/aburn01/Tools/Reports/%j.err
#SBATCH --constraint=rhel7

module load fastqc
source /cluster/tufts/micr_coffin01/aburn01/GTEx_HML2/Variables.txt
echo $1

for f in $(find ${1} -name "*_1.fastq.gz"); do
  	FILENAME=$(echo $f | rev | cut -c 12- | rev | uniq)
	gunzip ${FILENAME}_1.fastq.gz
	gunzip ${FILENAME}_2.fastq.gz
	java -jar $Trim_Path PE ${FILENAME}_1.fastq  ${FILENAME}_2.fastq  ${FILENAME}_f_p.fastq ${FILENAME}_f_u.fastq ${FILENAME}_r_p.fastq  ${FILENAME}_r_u.fastq  \
	ILLUMINACLIP:$Clip_Path:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30  MINLEN:75 2> ${FILENAME}_trim.txt 
	
	fastqc ${FILENAME}_1.fastq ${FILENAME}_2.fastq
done	
	
 





