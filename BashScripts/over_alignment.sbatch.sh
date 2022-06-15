#!/bin/bash
#SBATCH --job-name=Align		# Job name
#SBATCH --mail-type=END,FAIL		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=aidan.burn@tufts.edu
#SBATCH --partition=batch
#SBATCH --nodes=1			# Run a single task
#SBATCH --cpus-per-task=3		# Number of CPU cores per task
#SBATCH --mem-per-cpu=32000
#SBATCH --time=10:00:00
#SBATCH --output=/cluster/tufts/micr_coffin01/aburn01/Tools/Reports/%j.out
#SBATCH --error=/cluster/tufts/micr_coffin01/aburn01/Tools/Reports/%j.err
#SBATCH --constraint=rhel7

module load hisat
source /cluster/tufts/micr_coffin01/aburn01/GTEx_HML2/Variables.txt

for f in $(find ${1} -name "*_f_p.fastq"); do
  	FILENAME=$(echo $f |  rev | cut -c 11- | rev | uniq)
  	echo ${FILENAME}
	hisat2 --phred33 --dta-cufflinks -p 2 \
	  -x $Index_Dir/genome_tran -1 ${FILENAME}_f_p.fastq -2 ${FILENAME}_r_p.fastq \
	  -U ${FILENAME}_f_u.fastq,${FILENAME}_r_u.fastq \
	  -S ${FILENAME}_HISAT2_hg38_align.sam --un ${FILENAME}_HISAT2_hg38_unalign --un-conc ${FILENAME}_HISAT2_hg38_unalign_conc_mate 2> ${FILENAME}_align.txt
done

for dir in ${1}/*/; do
	echo $dir
	cd $dir
	mv * ../
done
