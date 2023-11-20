#!/bin/bash
#SBATCH -p genouest,ecobio
#SBATCH --cpus-per-task 4
. /local/env/envconda.sh 



### convert arguments into lists
IFS=" " read -ra list_path_input <<< $1
IFS=" " read -ra list_ids <<< $2
suffix=$3
path_out=$4
tmp_file=$6
conda_env=$5

conda activate ${conda_env}
### assign matching file for sarray
path_in=${list_path_input[${SLURM_ARRAY_TASK_ID}]}
id=${list_ids[${SLURM_ARRAY_TASK_ID}]}

touch $tmp_file

### creating directories
output_dir="${path_out}/${suffix}/${id}/decontamination/"

mkdir $output_dir"bams/"
mkdir $output_dir"fastqs/"

### input files
bam_file=$output_dir$id"_"$suffix"_sorted.bam"
bed_file=$output_dir$id"_"$suffix".bed"
conda deactivate
conda activate samtools
### reading bedfile line by line, generates a bam, sorts it and 2 paired fastq files for each outlier
count_bam=0
while IFS= read -r line
do
    ### generating bam file
    bam_junk="${output_dir}bams/${id}_${suffix}_junk_${count_bam}.bam"
    echo -e "samtools view $bam_file -L <(echo -e $line) -b -h -o -P $bam_junk \n"
    samtools view "$bam_file" -L <(echo -e "$line") -b -h -P -o $bam_junk   # -U $bam_remaining

    ### sorting bam file
    bam_junk_sorted="${output_dir}bams/${id}_${suffix}_junk_${count_bam}_sorted.bam"
    echo "samtools sort -n -o $bam_junk_sorted $bam_junk"
    samtools sort -n -o $bam_junk_sorted $bam_junk

    ### converting to fastq
    fastq_junk_1="${output_dir}fastqs/${id}_${suffix}_junk_${count_bam}_R1.fastq"
    fastq_junk_2="${output_dir}fastqs/${id}_${suffix}_junk_${count_bam}_R2.fastq"
    echo -e "bedtools bamtofastq -i $bam_junk_sorted -fq $fastq_junk_1 -fq2 $fastq_junk_2 \n"
    bedtools bamtofastq -i $bam_junk_sorted -fq $fastq_junk_1 -fq2 $fastq_junk_2
    
    ((count_bam++))   # it's important to start at 0 as all sarrays start at 0 too, $count_bam will serve as an index for next operations
done < "$bed_file"

echo "sample ${SLURM_ARRAY_TASK_ID} : $id" >> $tmp_file 



