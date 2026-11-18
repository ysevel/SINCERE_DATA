#!/bin/bash

source Configuration_file.txt
cd ${sincere_path}/scripts

##Array for jobs parallelization 
myArray=($(<$samples))

##Parameters creation
if [ "$krk_mem" = true ]
  then
    mem="-m"
  else
    mem=""
fi

if [ "$fasta_outliers" = true ]
  then
    fo="-f"
  else
    fo=""
fi

if [ "$annotation" = true ]
  then
    bakta_cmd="-a -ba $bakta_db"
  else
    bakta_cmd=""
fi

if [ "$junk_analysis"  = true ]
  then
    platon_cmd="-j -p $platon_db"
  else
    platon_cmd=""
fi

#Command lines
##Step 1 (Reads filtering)
cmd_step1="python Sincere-data_reads_filtering.py -i $reads_dir -d ${myArray[${SLURM_ARRAY_TASK_ID}]} -o $run_dir -n $run_name -k $krk_db -b $bl_file -t $nb_threads -cs $confidence_score $mem"

##Optional step (Splashome analysis)
cmd_sp="python splashome_filter.py -i $samples -o $run_dir -n $run_name -b $bl_file -ta $threshold_anomaly"

bl_file_update="${run_dir}/meta_analysis/Blacklist_filter_update.txt"
if [ "$splashome" = true ]
  then
    bl_file=$bl_file_update #Update blacklist file after splashome analysis
    splash="--splash" #Splashome parameter that will be added to step 2
  else
    splash=""
fi

##Step 2 (Outlier prediction)
cmd_step2="python Sincere-data_outlier_predictor.py -i $run_dir -d ${myArray[${SLURM_ARRAY_TASK_ID}]} -o $run_dir -n $run_name -k $krk_db -b $bl_file -t $nb_threads -cs $confidence_score $fo $splash $mem"

##Step 3 (Meta analysis)
cmd_step3="python Sincere-data_meta_analysis.py -i $samples -o $run_dir -n $run_name -b $bl_file -ts $threshold_sample -to $threshold_outlier -e $email"

##Step 4 (Decontamination)
cmd_step4="python Sincere-data_decontamination.py -i $run_dir -d ${myArray[${SLURM_ARRAY_TASK_ID}]} -o $run_dir -n $run_name -k $krk_db -c $contaminants -t $nb_threads -cs $confidence_score $bakta_cmd $mem $platon_cmd"

##Optional step (Assembly filtering)
cmd_af="python Sincere-data_assembly_filtering.py -i $samples -o $run_dir -n $run_name -t $nb_threads -co $completeness_threshold -ct $contamination_threshold"

#Script execution
if [ "$conda_env" != "None" ]
  then
  . /local/env/envconda.sh #This command might be specific of the GenOuest cluster maybe it will need modification if working on another cluster
  conda activate $conda_env
fi

##Scripts division in the different steps
case "${1:-}" in
  reads_filter)
    $cmd_step1
    ;;
  splash)
    $cmd_sp
    ;;
  outlier)
    $cmd_step2
    ;;
  meta)
    $cmd_step3
    ;;
  decontam)
    $cmd_step4
    ;;
  assembly_filter)
    $cmd_af
    ;;
  *)
    echo "Usage: $0 {reads_filter|splash|outlier|meta|decontam|assembly_filter}"
    exit 1
    ;;
esac