#!/bin/bash
#SBATCH -p genouest,ecobio
#SBATCH --cpus-per-task 1
. /local/env/envconda.sh 

nbThreads="10"



### convert arguments into lists
IFS=" " read -ra list_path_input <<< $1
IFS=" " read -ra list_ids <<< $2
suffix=$3
path_out=$4
contaminants=$5
tmp_file=$7
conda_env=$6

conda activate ${conda_env}
path_in=${list_path_input[${SLURM_ARRAY_TASK_ID}]}
id=${list_ids[${SLURM_ARRAY_TASK_ID}]}

### input files
input_dir="${path_out}/${suffix}/${id}/decontamination/"
kraken_result="${input_dir}kraken_results/${id}_${suffix}_trimmed_result.txt"
kraken_report="${input_dir}kraken_reports/${id}_${suffix}_trimmed_report.txt"
reads_R1=$(find $path_in -name "$id*1.fastq*") 
reads_R2=$(find $path_in -name "$id*2.fastq*") 

### output files
output_dir="${path_out}/${suffix}/${id}/final_assembly/"
mkdir $output_dir
output1=$output_dir$id"_"$suffix"_decontaxo_R1.fastq"
output2=$output_dir$id"_"$suffix"_decontaxo_R2.fastq"

touch $tmp_file


echo "python /scratch/ysevellec/script/KrakenTools-master/extract_kraken_reads.py -1 $reads_R1 -2 $reads_R2 -k $kraken_result -o $output1 -o2 $output2 -t $contaminants --report $kraken_report --include-children --exclude --fastq-output"
python /scratch/ysevellec/script/KrakenTools-master/extract_kraken_reads.py -1 $reads_R1 -2 $reads_R2 -k $kraken_result -o $output1 -o2 $output2 -t $contaminants --report $kraken_report --include-children --exclude --fastq-output

echo "sample ${SLURM_ARRAY_TASK_ID} : $id" >> $tmp_file
