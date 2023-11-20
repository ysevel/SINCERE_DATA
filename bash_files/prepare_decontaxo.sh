#!/bin/bash
#SBATCH -p genouest,ecobio
#SBATCH --cpus-per-task 8
. /local/env/envconda.sh 

nbThreads="32"

### convert arguments into lists
IFS=" " read -ra list_path_input <<< $1
IFS=" " read -ra list_ids <<< $2
suffix=$3
path_out=$4
tmp_file=$6
conda_env=$5

### assign variables related to the sample concerned
path_in=${list_path_input[${SLURM_ARRAY_TASK_ID}]}
id=${list_ids[${SLURM_ARRAY_TASK_ID}]}

output_dir="${path_out}${suffix}/${id}/decontamination/"
mkdir $output_dir

touch $tmp_file

conda activate ${conda_env}

### assembly
reads_R1=$(find $path_in -name "${id}*1.fastq*") 
echo "find $path_in -name ${id}*1.fastq*"
reads_R2=$(find $path_in -name "${id}*2.fastq*") 
echo "find $path_in -name ${id}*2.fastq*"

echo "spades.py --sc --careful -t $nbThreads -1 $reads_R1 -2 $reads_R2 -o ${output_dir}${id}_${suffix}_spades" ;
spades.py --sc --careful -t $nbThreads -1 $reads_R1 -2 $reads_R2 -o "${output_dir}${id}_${suffix}_spades" ;



### ecremate assemblies
cp $output_dir$id"_"$suffix"_spades/contigs.fasta" $output_dir$id"_"$suffix"_contigs.fasta" 
assembly="$output_dir${id}_${suffix}_ecrem_contigs.fasta"
echo $assembly
echo "reformat.sh minlength=500 in="$output_dir${id}_${suffix}_contigs.fasta"  out=$assembly overwrite=true"
reformat.sh minlength=500 in="$output_dir${id}_${suffix}_contigs.fasta"  out=$assembly overwrite=true

### build index 
index=$output_dir$id"_"$suffix"_index"
echo "bowtie2-build -f $assembly $index"
bowtie2-build -f $assembly $index


### get reads and bowtie alignment
samfile=$output_dir$id"_"$suffix".sam"
reads_R1=$(find $path_in -name "*$id*1.fastq*") 
reads_R2=$(find $path_in -name "*$id*2.fastq*") 
echo "bowtie2 --threads $nbThreads -1 $reads_R1 -2 $reads_R2 -x $index --sensitive -S $samfile"
bowtie2 --threads $nbThreads -1 $reads_R1 -2 $reads_R2 -x $index --sensitive -S $samfile 

conda deactivate
conda activate samtools
# samtools view
bamfile=$output_dir$id"_"$suffix".bam"
echo "samtools view -bS -o $bamfile $samfile"
samtools view -bS -o $bamfile $samfile

# samtool sort
sorted_bamfile=$output_dir$id"_"$suffix"_sorted.bam"
echo "samtools sort $bamfile -T ${id}_temp_file -o $sorted_bamfile"
samtools sort $bamfile -T $id"_temp_file" -o $sorted_bamfile

# samtools index
$baifile=$output_dir$id"_"$suffix".bam.bai"
echo "samtools index -b $sorted_bamfile $baifile"
samtools index -b $sorted_bamfile $baifile

# bedtools genomecov (bedgraphs files)
bedgraphfile=$output_dir$id"_"$suffix"_bedgraph.tsv"
echo "bedtools genomecov -ibam $sorted_bamfile -bg > $bedgraphfile"
bedtools genomecov -ibam $sorted_bamfile -bg > $bedgraphfile

echo "sample ${SLURM_ARRAY_TASK_ID} : $id" >> $tmp_file
