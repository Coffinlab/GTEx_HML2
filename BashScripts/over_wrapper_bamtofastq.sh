#!/bin/bash
source /cluster/tufts/micr_coffin01/aburn01/GTEx_HML2/Start_Variables.txt

for i in ${1}/*.Aligned.sortedByCoord.out.patched.md.bam
do
 
  base_name=$(basename $i)
  dir_name=$(dirname $i)
 
   
  sample_prefix=${base_name%.Aligned.sortedByCoord.out.patched.md.bam}
  out_dir=${dir_name}/${sample_prefix}
  mkdir -p $out_dir
  
  echo "submitting file $i, sample prefix $sample_prefix, out dir $out_dir"
  
  sbatch over_bamtofastq.sh $i $sample_prefix $out_dir

done
