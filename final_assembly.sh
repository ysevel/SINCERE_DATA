#!/bin/bash
#SBATCH -p genouest,ecobio
#SBATCH --cpus-per-task 8
. /local/env/envconda.sh 

nbThreads="16"

### convert arguments into lists
IFS=" " read -ra list_path_input <<< $1
IFS=" " read -ra list_ids <<< $2
suffix=$3
path_out=$4
tmp_file=$5

id=${list_ids[${SLURM_ARRAY_TASK_ID}]}
input_dir="${path_out}${suffix}/${id}/final_assembly/"

touch $tmp_file

conda activate my_pip

### assemble
reads_R1=$input_dir$id"_"$suffix"_decontaxo_R1.fastq"
reads_R2=$input_dir$id"_"$suffix"_decontaxo_R2.fastq"
out_spades="${input_dir}${id}_${suffix}_spades"
echo "spades.py --sc --careful -t ${nbThreads} -1 $reads_R1 -2 $reads_R2 -o ${input_dir}${id}_${suffix}_spades"
spades.py --sc --careful -t $nbThreads -1 $reads_R1 -2 $reads_R2 -o $out_spades ;

## FastQC on reads
run_dir="${path_out}${suffix}/"
fastqc_dir="${run_dir}analysis_on_final/fastQC/"
mkdir $fastqc_dir
echo "fastqc -t $nbThreads -o $fastqc_dir $reads_R1"
fastqc -t $nbThreads -o $fastqc_dir $reads_R1
fastqc -t $nbThreads -o $fastqc_dir $reads_R2

### ecremate assemblies 
assembly="${input_dir}${id}_${suffix}_contigs.fasta" 
cp $input_dir$id"_"$suffix"_spades/contigs.fasta" $assembly

assembly_ecrem="$input_dir${id}_${suffix}_ecrem_contigs.fasta"
echo "reformat.sh minlength=500 in=$assembly out=$assembly_ecrem overwrite=true"
reformat.sh minlength=500 in=$assembly out=$assembly_ecrem overwrite=true

## quast on it
quast_dir="${run_dir}analysis_on_final/quast/"
mkdir $quast_dir
quast $assembly_ecrem -o "$quast_dir${id}_${suffix}"

echo "sample ${SLURM_ARRAY_TASK_ID} : $id" >> $tmp_file
