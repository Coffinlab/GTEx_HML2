#!/bin/bash
#SBATCH --job-name=GTEX_Pipeline                # Job name
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=aidan.burn@tufts.edu
#SBATCH --partition=largemem
#SBATCH --nodes=1                       # Run a single task
#SBATCH --cpus-per-task=1               # Number of CPU cores per task
#SBATCH --mem-per-cpu=5000
#SBATCH --time=1-24:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --constraint=rhel7

source /cluster/tufts/micr_coffin01/aburn01/GTEx_HML2/Variables.txt
echo ${File} ${File2} ${File3} ${File4}


#echo "Running trim of ${File}"
#sbatch over_read_filter_quality.sbatch.sh $Sample_Dir

#echo "Running trim of ${File2}"
#sbatch over_read_filter_quality.sbatch.sh $Sample_Dir_2

#echo "Running trim of ${File3}"
#sbatch over_read_filter_quality.sbatch.sh $Sample_Dir_3

#echo "Running trim of ${File4}"
#sbatch -W  over_read_filter_quality.sbatch.sh $Sample_Dir_4

wait

echo "Running alignment of ${File}"
sbatch  over_alignment.sbatch.sh $Sample_Dir 

#echo "Running alignment of ${File2}"
#sbatch  over_alignment.sbatch.sh $Sample_Dir_2

#echo "Running alignment of ${File3}"
#sbatch  over_alignment.sbatch.sh $Sample_Dir_3 

#echo "Running alignment of ${File4}"
#sbatch -W  over_alignment.sbatch.sh $Sample_Dir_4 

wait

echo "Running telescope counts of ${File}"
sbatch over_Telescope.sbatch.sh $Sample_Dir

#echo "Running telescope counts of ${File2}"
#sbatch over_Telescope.sbatch.sh $Sample_Dir_2

#echo "Running telescope counts of ${File3}"
#sbatch over_Telescope.sbatch.sh $Sample_Dir_3

#echo "Running telescope counts of ${File4}"
#sbatch -W  over_Telescope.sbatch.sh $Sample_Dir_4

wait

echo "Collecting reports of ${File}"
sh over_report.sh $Sample_Dir

#echo "Collecting reports of ${File2}"
#sh over_report.sh $Sample_Dir_2

#echo "Collecting reports of ${File3}"
#sh over_report.sh $Sample_Dir_3

#echo "Collecting reports of ${File4}"
#sh over_report.sh $Sample_Dir_4
