#!/bin/bash

source Configuration_file.txt

myArray=($(<$samples))

#Number of samples to analyze
nb_samples=${#myArray[@]}

nb_array=0-${nb_samples}%${max_sim_jobs}

out_path=${log_output}/%A_%a.%N.${run_name}

step1=$(sbatch --parsable --array=${nb_array} --cpus-per-task=${nb_threads} --mem=${big_mem} -o ${out_path}_step1.out -e ${out_path}_step1.err Launch_sincere.sh reads_filter)

if [ "$splashome" = true ]
    then
        step_sp=$(sbatch --parsable --dependency=afterany:${step1} --cpus-per-task=${nb_threads} --mem=${small_mem} -o ${out_path}_step_sp.out -e ${out_path}_step_sp.err Launch_sincere.sh splash)
        step2=$(sbatch --parsable --array=${nb_array} --dependency=afterany:${step_sp} --cpus-per-task=${nb_threads} --mem=${big_mem} -o ${out_path}_step2.out -e ${out_path}_step2.err Launch_sincere.sh outlier)
    else
        step2=$(sbatch --parsable --array=${nb_array} --dependency=afterany:${step1} --cpus-per-task=${nb_threads} --mem=${big_mem} -o ${out_path}_step2.out -e ${out_path}_step2.err Launch_sincere.sh outlier)
fi

step3=$(sbatch --parsable --dependency=afterany:${step2} --cpus-per-task=${nb_threads} --mem=${small_mem} -o ${out_path}_step3.out -e ${out_path}_step3.err Launch_sincere.sh meta)
step4=$(sbatch --parsable --array=${nb_array} --dependency=afterany:${step3} --cpus-per-task=${nb_threads} --mem=${big_mem} -o ${out_path}_step4.out -e ${out_path}_step4.err Launch_sincere.sh decontam)

if [ "$assembly_filtering" = true ]
    then
        sbatch --dependency=afterany:${step4} --cpus-per-task=${nb_threads} --mem=${small_mem} -o ${out_path}_step_af.out -e ${out_path}_step_af.err Launch_sincere.sh assembly_filter
fi
