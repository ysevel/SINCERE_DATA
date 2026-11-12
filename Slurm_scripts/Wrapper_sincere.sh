#!/bin/bash

source Configuration_file.txt
cd ${sincere_path}/Slurm_scripts

myArray=($(<$samples))

#Number of samples to analyze
nb_samples=${#myArray[@]}

nb_array=0-${nb_samples}%${max_sim_jobs}

out_log=.${log_output}/%A_%a.%N.%x.out
error_log=.${log_output}/%A_%a.%N.%x.err

step1=$(sbatch --parsable --array=${nb_array} --cpus-per-task=${nb_threads} --mem=${big_mem} Launch_sincere.sh reads_filter)

if [ "$splashome" = true ]
    then
        step_sp=$(sbatch --parsable --dependency=afterany:${step1} --cpus-per-task=${nb_threads} --mem=${small_mem} Launch_sincere.sh splash)
        step2=$(sbatch --parsable --array=${nb_array} --dependency=afterany:${step_sp} --cpus-per-task=${nb_threads} --mem=${big_mem} Launch_sincere.sh outlier)
    else
        step2=$(sbatch --parsable --array=${nb_array} --dependency=afterany:${step1} --cpus-per-task=${nb_threads} --mem=${big_mem} Launch_sincere.sh outlier)
fi

step3=$(sbatch --parsable --dependency=afterany:${step2} --cpus-per-task=${nb_threads} --mem=${small_mem} Launch_sincere.sh meta)
step4=$(sbatch --parsable --array=${nb_array} --dependency=afterany:${step3} --cpus-per-task=${nb_threads} --mem=${big_mem} Launch_sincere.sh decontam)

if [ "$assembly_filtering" = true ]
    then
        sbatch --dependency=afterany:${step4} --cpus-per-task=${nb_threads} --mem=${small_mem} Launch_sincere.sh assembly_filter
fi