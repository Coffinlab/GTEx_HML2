#!/bin/bash
#SBATCH --job-name=GTEX_Pipeline		# Job name
#SBATCH --mail-type=END,FAIL		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=aidan.burn@tufts.edu
#SBATCH --partition=preempt
#SBATCH --nodes=1			# Run a single task
#SBATCH --cpus-per-task=1		# Number of CPU cores per task
#SBATCH --mem-per-cpu=5000
#SBATCH --time=1-24:00:00
#SBATCH --output=/cluster/tufts/micr_coffin01/aburn01/Tools/Reports/%j.out
#SBATCH --error=/cluster/tufts/micr_coffin01/aburn01/Tools/Reports/%j.err

source /cluster/tufts/micr_coffin01/aburn01/GTEx_HML2/Start_Variables.txt
echo ${File}

#echo "Running Download of ${File}"
#sh over_anvil_download.sh $File

echo "Running bamtofastq of ${File}"
sbatch -W over_wrapper_bamtofastq.sh $Sample_Dir &

#wait

#echo "Running Download of ${File2}"
#sh over_anvil_download.sh $File2

#echo "Running bamtofastq of ${File2}"
#sbatch -W over_wrapper_bamtofastq.sh $Sample_Dir_2 &

#wait

#echo "Running Download of ${File3}"
#sh over_anvil_download.sh $File3

#echo "Running bamtofastq of ${File3}"
#sbatch -W over_wrapper_bamtofastq.sh $Sample_Dir_3 &

#wait

#echo "Running Download of ${File4}"
#sh over_anvil_download.sh $File4

#echo "Running bamtofastq of ${File4}"
#sbatch -W over_wrapper_bamtofastq.sh $Sample_Dir_4
