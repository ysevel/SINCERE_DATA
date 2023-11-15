#!/bin/bash
#SBATCH --cpus-per-task 8
. /local/env/envconda.sh 

nbThreads="16"

### convert arguments into lists
input=$1
IFS=" " read -ra list_ids <<< $2
suffix=$3
path_out=$4
tmp_file=$5

id=${list_ids[${SLURM_ARRAY_TASK_ID}]}

conda activate checkm

## get final assembly 
path_final="${input}${id}/final_assembly/"
assembly=$(find $path_final -name "$id*_ecrem_contigs.fasta") 
echo "assembly : $assembly"

## copy it to a general directory
assem_dir="${input}analysis_on_final/final_assemblies/"
mkdir $assem_dir
assembly_name=${assembly##*/} 
assembly_name=$(echo $assembly_name | sed -e "/_ecrem_contigs/_final_assembly/g");
echo "assembly_name : $assembly_name"
cp $assembly "${assem_dir}${assembly_name}"

## checkM on it
checkm_dir="${input}analysis_on_final/checkM/"
mkdir $checkm_dir
echo "checkm lineage_wf -t $nbThreads -x fasta $assem_dir $checkm_dir"
checkm lineage_wf -t $nbThreads -x fasta $assem_dir $checkm_dir

echo "sample ${SLURM_ARRAY_TASK_ID} : $id" >> $tmp_file

