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
tmp_file=$5

path_in=${list_path_input[${SLURM_ARRAY_TASK_ID}]}
id=${list_ids[${SLURM_ARRAY_TASK_ID}]}

conda activate my_pip

reads_R1=$(find $path_in -name "$id*1.fastq*") 
reads_R2=$(find $path_in -name "$id*2.fastq*") 

kraken_bdd="/groups/ecogeno/DEV/2023_kraken_db/k2_pluspfp"
output_dir="${path_out}/${suffix}/${id}/decontamination/"
result_dir=$output_dir"kraken_results/"
report_dir=$output_dir"kraken_reports/"
mkdir $result_dir
mkdir $report_dir

result="${result_dir}${id}_${suffix}_trimmed_result.txt"
report="${report_dir}${id}_${suffix}_trimmed_report.txt"
echo "kraken2 --db $kraken_bdd --paired $reads_R1 $reads_R2 --threads $nbThreads --output $result --report $report --use-names"
kraken2 --db $kraken_bdd --paired $reads_R1 $reads_R2 --threads $nbThreads --output $result --report $report --use-names

echo "sample ${SLURM_ARRAY_TASK_ID} : $id" >> $tmp_file

