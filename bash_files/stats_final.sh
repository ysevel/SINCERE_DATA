#!/bin/bash
#SBATCH -p genouest,bigmem,ecobio
#SBATCH --cpus-per-task 8
. /local/env/envconda.sh 

nbThreads="16"

### convert arguments into lists
input=$1
IFS=" " read -ra list_ids <<< $2
suffix=$3
path_out=$4
tmp_file=$6
conda_env=$5

id=${list_ids[${SLURM_ARRAY_TASK_ID}]}

conda activate ${conda_env}

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
echo "sample ${SLURM_ARRAY_TASK_ID} : $id" >> $tmp_file
