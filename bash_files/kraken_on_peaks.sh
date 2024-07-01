#!/bin/bash

#SBATCH --cpus-per-task 8
#SBATCH --mem=300G
. /local/env/envconda.sh 

nbThreads="16"

IFS=" " read -ra list_path_input <<< $1
IFS=" " read -ra list_ids <<< $2
suffix=$3
path_out=$4
tmp_file=$7
conda_env=$5

id=${list_ids[${SLURM_ARRAY_TASK_ID}]}

conda activate ${conda_env}

output_dir="${path_out}${suffix}/${id}/decontamination/"
kraken_bdd=$6

reads_dir="${output_dir}fastqs"
report_dir="${output_dir}kraken_reports"
result_dir="${output_dir}kraken_results"
mkdir $report_dir
mkdir $result_dir

rm ${report_dir}/*junk*.fastq

for reads1 in ${reads_dir}/*junk*R1.fastq
do 
    reads=$(echo $reads1 | sed 's/_R1//')
    reads2=$(echo $reads1 | sed 's/R1/R2/')
    report=$(echo $reads | sed 's/\.fastq$/_report.txt/' | sed 's/fastqs/kraken_reports/')
    result=$(echo $reads | sed 's/\.fastq$/_result.txt/' | sed 's/fastqs/kraken_results/')
    echo -e "$reads1 \n$reads2 \n$report \n$result"
    echo "kraken2 --db $kraken_bdd --paired $reads1 $reads2 --threads $nbThreads --output $result --report $report --use-names"
    kraken2 --db $kraken_bdd --paired $reads1 $reads2 --threads $nbThreads --output $result --report $report --use-names
done

echo "sample ${SLURM_ARRAY_TASK_ID} : $id" >> $tmp_file

