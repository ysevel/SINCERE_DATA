#!/bin/bash

#SBATCH --cpus-per-task 8
#SBATCH --mem=300G
. /local/env/envconda.sh 

nbThreads="16"

### convert arguments into lists
IFS=" " read -ra list_path_input <<< $1
IFS=" " read -ra list_ids <<< $2
suffix=$3
path_out=$4
tmp_file=$7
conda_env=$5

path_in=${list_path_input[${SLURM_ARRAY_TASK_ID}]}
id=${list_ids[${SLURM_ARRAY_TASK_ID}]}

conda activate ${conda_env}

kraken_bdd=$6
assembly_path="${path_out}/${suffix}/${id}/decontamination/"
output_dir="${path_out}/${suffix}/analysis_assembly_on_trimmed/"
result_dir=$output_dir"kraken_results/"
report_dir=$output_dir"kraken_reports/"
mkdir $output_dir
mkdir $result_dir
mkdir $report_dir

assembly=$(find $assembly_path -name "$id*ecrem_contigs.fasta") 
echo "assembly : $assembly"
result_assem="${result_dir}${id}_${suffix}_assem_trim_result.txt"
report_assem="${report_dir}${id}_${suffix}_assem_trim_report.txt"
echo "kraken2 --db $kraken_bdd $assembly --threads $nbThreads --output $result_assem --report $report_assem --use-names"
kraken2 --db $kraken_bdd $assembly --threads $nbThreads --output $result_assem --report $report_assem --use-names

echo "sample ${SLURM_ARRAY_TASK_ID} : $id" >> $tmp_file

