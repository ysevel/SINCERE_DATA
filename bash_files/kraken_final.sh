#!/bin/bash
#SBATCH -p genouest,bigmem,ecobio
#SBATCH --cpus-per-task 8
#SBATCH --mem=300G
. /local/env/envconda.sh 

nbThreads="16"

# "{fastq_dir} {list_ids_sh} {suffix} "

### convert arguments into lists
input=$1
IFS=" " read -ra list_ids <<< $2
suffix=$3
path_out=$4
tmp_file=$5

id=${list_ids[${SLURM_ARRAY_TASK_ID}]}

conda activate my_pip
path_final="${input}${id}/final_assembly/"
reads_R1=$(find $path_final -name "$id*1.fastq*") 
echo "reads_R1 : $reads_R1"
reads_R2=$(find $path_final -name "$id*2.fastq*") 
echo "reads_R2 : $reads_R2"

kraken_bdd="/groups/ecogeno/DEV/2023_kraken_db/k2_pluspfp"
output_dir="${path_out}${suffix}/analysis_on_final/"
result_dir=$output_dir"kraken_results_decontaxo/"
report_dir=$output_dir"kraken_reports_decontaxo/"
mkdir $result_dir
mkdir $report_dir

result_reads="${result_dir}${id}_${suffix}_reads_result.txt"
report_reads="${report_dir}${id}_${suffix}_reads_report.txt"
echo "kraken2 --db $kraken_bdd --paired $reads_R1 $reads_R2 --threads $nbThreads --output $result_reads --report $report_reads --use-names"
kraken2 --db $kraken_bdd --paired $reads_R1 $reads_R2 --threads $nbThreads --output $result_reads --report $report_reads --use-names

assembly=$(find $path_final -name "$id*_ecrem_contigs.fasta") 
echo "assembly : $assembly"
result_assem="${result_dir}${id}_${suffix}_assem_result.txt"
report_assem="${report_dir}${id}_${suffix}_assem_report.txt"
echo "kraken2 --db $kraken_bdd $assembly --threads $nbThreads --output $result_assem --report $report_assem --use-names"
kraken2 --db $kraken_bdd $assembly --threads $nbThreads --output $result_assem --report $report_assem --use-names

echo "sample ${SLURM_ARRAY_TASK_ID} : $id" >> $tmp_file

