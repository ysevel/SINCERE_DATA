##!/usr/bin/python
import os,sys
import glob
import time
import math
from datetime import datetime
from pathlib import Path
from textwrap import dedent as twdd
import shutil
import argparse
import subprocess
import logging
#from concurrent.futures import ProcessPoolExecutor
#from functools import partial
#from multiprocessing import cpu_count, Manager, Lock


#functions

#argument handling
def handle_program_options():
    descr = """\
    Complete pipeline of single-cell decontamination using \
    taxonomic assignation. 
    Program arguments
    -----------------"""

    parser = argparse.ArgumentParser(description=twdd(descr),
                           formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input', required=True,
                        help="path to index file of samples to process. \
                          Needs to give for each sample its directory path \
                          and names, comma separated")
    parser.add_argument('-d', '--id', required=True,
                        help="id for the sample, will be used to point to the correct files")
    parser.add_argument('-o', '--output', required=True,
                        help="Output path where directories will be built")
    parser.add_argument('-t', '--threads', type=int ,default=8,
                        help="Number of threads")
    parser.add_argument('-m', '--mem', action='store_true',
                        help="Allow kraken to run on memory (optionnal, otherwise will run kraken2 with the --memory-mapping option that avoid loading kraken databese into the memory but require significantly higher computation time. If you want to load kraken2 database into the memory best to make sure you allocated enough memory the PFPplus kraken2 database require 300G and the nt database 750G to run properly for example).")
    parser.add_argument('-n', '--name', required=True,
                        help="suffix given to every output file for this run")
    parser.add_argument('-k', '--kraken_db', required=True,
                        help="path to the kraken database folder")
    parser.add_argument('-c', '--contaminants', required=False ,default=None,
                        help="Path to a text file containing the taxon that should be systematically removed from the dataset (one taxid per line). Carefull it will eliminate all the children taxa as well! This file would ideally be generated from the previous steps of this workflow, but it can also be a 'blacklist' file with taxa the user want to exclude.")
    parser.add_argument('-j', '--keep_junk_contigs' ,action='store_true',
                        help="if the option is selected the contigs discarded by the assembly cleaning will be stored and analyzed for possible mobiles elements")
    return parser.parse_args()


def remove_sample(id_sample, error, logger, output_missing_file):
    """ report the sample id of a missing sample and write the id and the corresponding error message to the missing file."""
    with missing_id_lock:
        current_missing_id_list.append(id_sample)
        with open(output_missing_file, 'a') as missf:
            missf.write(id_sample + '\t' + error + '\n')
    logger.error(f"Sample {id_sample} removed: {error}")


def check_sample(id_sample,list_missing_id,logger):
    if id_sample in list_missing_id:
        logger.warning(f"{id_sample} has been found missing in previous steps and will been ignored, exiting.\n")
        return False
    else:
        logger.info(f"processing on preliminary assembly on {id_sample}\n")
        return True


def check_reads(path_input, id_sample, logger, output_missing_file):
    """ from an input path, a sample id, the number of paired reads (either 1 or 2), test the path then check if read files are present 
        and returns the file name if validated. If not, False is returned to cut the decontamination workflow for the sample. """
    message = "" #note: it is best that the sample id is separated from any suffix in the reads by '_', otherwise finding the right read can be tricky.    
    input_R1=f"{path_input}/{id_sample}_*R1.fastq*"
    input_R2=f"{path_input}/{id_sample}_*R2.fastq*"
    R1_list=glob.glob(input_R1)
    logger.info(f"R1 list: {R1_list}")
    R2_list=glob.glob(input_R2)
    if os.path.exists(path_input):
        logger.info(f"Files in input path: {glob.glob(path_input)}")
        if len(R1_list) < 1:
            message = f"Error : R1 reads file not found for {id_sample}"
            remove_sample(id_sample, message, logger, output_missing_file)
            return False, False
        elif len(R2_list) < 1:
            message = f"Error : R2 reads file not found for {id_sample}"
            remove_sample(id_sample, message, logger, output_missing_file)
            return False, False
        else:
            R1_file=R1_list[0]
            R2_file=R2_list[0]
            logger.info(f"Both reads files found for sample {id_sample}")
            return R1_file, R2_file
    else:
        logger.error("Invalid input path specified in the index file. Please check if the path truly exists.")
        sys.exit(1)


def launch_kraken(R1_file, R2_file, id_sample, output, kraken_db, suffix, threads, mem, logger, sample_dir):
    """
    create the result directories and launch the kraken command on paired reads
    """
    try: #creation of two directories for reports and results from kraken2
        os.makedirs(f"{output}/kraken_results")
    except:
        f"Warning : sample directory {output}/kraken_results already existing. Overwriting..."
    try:
        os.makedirs(f"{output}/kraken_reports")
    except:
        f"Warning : sample directory {output}/kraken_reports already existing. Overwriting..." 
    kr_out=f"{output}/kraken_results/{id_sample}_{suffix}_kr_results.txt"
    kr_rep=f"{output}/kraken_reports/{id_sample}_{suffix}_kr_reports.txt"
    cmd_mem = "--memory-mapping" if not mem else ""
    command_kr=f"kraken2 --threads {threads} --db {kraken_db} --paired {R1_file} {R2_file} -output {kr_out} --use-names --report {kr_rep} {cmd_mem}"
    logger.info(f"Kraken command for reads: {command_kr}")
    log_kr=f"{sample_dir}{id_sample}_{suffix}_kraken2.log"
    with open(log_kr,'w') as log_k:
        subprocess.run(command_kr, stdout=log_k, stderr=subprocess.STDOUT, shell=True)
    return kr_out, kr_rep


def blacklist_filtering(R1_file, R2_file, id_sample, output, kr_filter_out, suffix, contaminants, kr_filter_reports, logger, sample_dir):
    """ launch reads decontamination from a list of contaminants using krakentools."""
    output_dir=f"{output}/reads"
    try:                                
        os.makedirs(f"{output_dir}")
    except:
        f"Warning : sample directory {output_dir} already existing. Overwriting..."

    filter_contam=f"python extract_kraken_reads.py -1 {R1_file} -2 {R2_file} -k {kr_filter_out} -o {output_dir}/{id_sample}_{suffix}_R1.fastq -o2 {output_dir}/{id_sample}_{suffix}_R2.fastq -t {contaminants} --include-children --exclude --report {kr_filter_reports} --fastq-output"
    logger.info(f"Command for filtering reads: \n{filter_contam}")
    log_filter=f"{sample_dir}{id_sample}_{suffix}_filter_reads.log"
    with open(log_filter,'w') as log_f:
        subprocess.run(filter_contam, stdout=log_f, stderr=subprocess.STDOUT, shell=True)
    return f"{output_dir}/{id_sample}_{suffix}_R1.fastq", f"{output_dir}/{id_sample}_{suffix}_R2.fastq"


def genome_assembly(R1_file, R2_file, id_sample, suffix, assembly_dir, threads, logger, sample_out_dir):
    """
    perform the assembly on the decontaminated reads, is supposed to be run with fixed parameters, also filter out contigs shorter than 500bp.
    """
    try:
        os.makedirs(assembly_dir) #create a specific assembly dir, it is intended, at this step to be in the final_assembly folder.
    except:
        print(f"Warning : sample directory {assembly_dir} already existing. Overwriting...")
    
    log_assembly=f"{sample_out_dir}{id_sample}_{suffix}_assembly.log"
    cmd_spades=f"spades.py --sc --careful -t {threads} -1 {R1_file} -2 {R2_file} -o {assembly_dir}/{id_sample}_{suffix}_spades" #run spade in a specific temporary folder that is supposed to be erased.
    logger.info(f"running spades with: {cmd_spades}")
    cmd_reformat=f"reformat.sh minlength=500 in={assembly_dir}/{id_sample}_{suffix}_spades/contigs.fasta  out={assembly_dir}/{id_sample}_{suffix}_raw_assembly.fasta overwrite=true" #the trimmed assembly is saved directly in the assembly folder, outside of spades folder
    with open(log_assembly,'w') as log_a:
        subprocess.run(cmd_spades, stdout=log_a,stderr=subprocess.STDOUT, shell=True)
        subprocess.run(cmd_reformat, stdout=log_a,stderr=subprocess.STDOUT, shell=True) #the filtered assembly is writen directly in the assembly_dir
        if os.path.isdir(f"{assembly_dir}/{id_sample}_{suffix}_spades") :
            shutil.rmtree(f"{assembly_dir}/{id_sample}_{suffix}_spades") 
    logger.info(f"The preliminary assembly can be found in : {assembly_dir}/{id_sample}_{suffix}_prelim_assembly.fasta")
    return f"{assembly_dir}/{id_sample}_{suffix}_raw_assembly.fasta"
def junk_contig_analysis(final_assembly, contaminants, id_sample, suffix, threads, logger, sample_out_dir, mem, sample_dir)
    """extract the junk contigs for analysis"""
    junk_dir=f"{sample_out_dir}cleaned_assembly" #create the clean directory which will contain the final output of each sample.
    try:
        os.makedirs(junk_dir)
    except:
        f"Warning : sample directory {junk_dir} already existing. Overwriting..."
    log_junk=f"{sample_dir}{id_sample}_{suffix}_discarded_contigs_analysis.log"
    junk_contigs=f"{junk_dir}/{id_sample}_{suffix}_junk_contigs.fasta"
    cmd_junk_extract=cmd_krtools=f"python extract_kraken_reads.py -s {junk_contigs} -k {kr_p_out}/{id_sample}_{suffix}_kr_results.txt -o {junk_assembly} -t {contaminants} --include-children --report {kr_p_rep}/{id_sample}_{suffix}_kr_reports.txt"
    cmd_quast=f"quast -t {threads} --pe1 {decontam_reads_R1} --pe2 {decontam_reads_R2} -o {cleaned_dir}/quast/{id_sample}_cleaned_assembly {cleaned_assembly}"
    with open(log_junk,'w') as log_j:
        logger.info(f"running junk contig extraction... krakentools command: {cmd_junk_extract}")
        subprocess.run(cmd_junk_extract, stdout=log_j,stderr=subprocess.STDOUT, shell=True)
        logger.info(f"running  quast on discraded contigs: {cmd_quast}")
        subprocess.run(cmd_quast, stdout=log_j,stderr=subprocess.STDOUT, shell=True)

def contigs_cleaning(decontam_reads_R1, decontam_reads_R2, final_assembly, contaminants, id_sample, suffix, threads, logger, sample_out_dir, mem, sample_dir):
    """perform taxonomic classification on contigs and decontamination based on the the contaminants list defined at the previous stage of the pipeline, intended to be optionnal, it will be implemented later on."""
    #specific variables
    final_dir=f"{sample_out_dir}final_assembly"
    cleaned_dir=f"{sample_out_dir}cleaned_assembly" #create the clean directory which will contain the final output of each sample.
    try:
        os.makedirs(cleaned_dir)
    except:
        f"Warning : sample directory {cleaned_dir} already existing. Overwriting..."
    log_clean=f"{sample_dir}{id_sample}_{suffix}_contigs_clean.log"
    kr_p_rep=f"{final_dir}kraken_raw_assembly_reports"
    try:
        os.makedirs(kr_p_rep)
    except:
        f"Warning : sample directory {kr_p_rep} already existing. Overwriting..."
    kr_p_out=f"{final_dir}kraken_raw_assembly_results"
    try:
        os.makedirs(kr_p_out)
    except:
        f"Warning : sample directory {kr_p_out} already existing. Overwriting..."
    cmd_mem = "--memory-mapping" if not mem else ""
    cmd_kraken=f"kraken2 --threads {threads} --db {kraken_db} {final_assembly} -output {kr_p_out}/{id_sample}_{suffix}_kr_results.txt --use-names --report {kr_p_rep}/{id_sample}_{suffix}_kr_reports.txt {cmd_mem}"
    cleaned_assembly=f"{cleaned_dir}/{id_sample}_{suffix}_cleaned_assembly.fasta"
    cmd_krtools=f"python extract_kraken_reads.py -s {final_assembly} -k {kr_p_out}/{id_sample}_{suffix}_kr_results.txt -o {cleaned_assembly} -t {contaminants} --include-children --exclude --report {kr_p_rep}/{id_sample}_{suffix}_kr_reports.txt"
    with open(log_clean,'w') as log_c:
        subprocess.run(cmd_kraken, stdout=log_c,stderr=subprocess.STDOUT, shell=True)
        logger.info(f"krakentools command: {cmd_krtools}")
        subprocess.run(cmd_krtools, stdout=log_c,stderr=subprocess.STDOUT, shell=True)
        kr_c_rep=f"{cleaned_dir}/kraken_clean_assembly_reports"
        try:
            os.makedirs(kr_c_rep)
        except:
             print(f"Warning : sample directory {kr_c_rep} already existing. Overwriting...")
        kr_c_out=f"{cleaned_dir}/kraken_clean_assembly_results"
        try:
            os.makedirs(kr_c_out)
        except:
            print(f"Warning : sample directory {kr_c_out} already existing. Overwriting...")
        cmd_kr_clean=f"kraken2 --threads {threads} --db {kraken_db} {cleaned_assembly} -output {kr_c_out}/{id_sample}_{suffix}_kr_results.txt --use-names --report {kr_c_rep}/{id_sample}_{suffix}_kr_reports.txt {cmd_mem}"
        subprocess.run(cmd_kr_clean, stdout=log_c,stderr=subprocess.STDOUT, shell=True)
        cmd_quast=f"quast -t {threads} --pe1 {decontam_reads_R1} --pe2 {decontam_reads_R2} -o {cleaned_dir}/quast/{id_sample}_cleaned_assembly {cleaned_assembly}"
        subprocess.run(cmd_quast, stdout=log_c,stderr=subprocess.STDOUT, shell=True)
    return cleaned_assembly
#variable initiation

def process_sample(id_sample, path_input, suffix, kraken_db, run_dir, nb_threads, mem, previous_missing_id_list, output_missing_file, contaminants):
    ###check missing sample
    sample_dir=f"{path_input}/{id_sample}/"
    sample_out_dir=f"{run_dir}/{id_sample}/"
    taxon_dir= f"{path_input}/meta_analysis/"
    filter_dir=f"{sample_dir}/filtered_reads"
    decontam_dir=f"{sample_out_dir}/decontamination"
    
    log_file = f"{sample_out_dir}decontamination_{id_sample}_{suffix}.log"
    print(f"All prints will now be saved in {log_file} file.")
    
    if os.path.exists(log_file):    # generates a new file
        os.remove(log_file)
    
    logger = logging.getLogger(id_sample)
    logger.setLevel(logging.INFO)
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    logger.info(f"read filtering begin for sample {id_sample}")
    if not check_sample(id_sample, previous_missing_id_list, logger):
        return
    
    if os.path.exists(sample_dir):
        logger.info(f"Running the outlier detection for {id_sample} in {sample_dir} directory")
        logger.info(f"Successfully found {id_sample}, proceding to read analysis.")
    else:
        if id_sample in current_missing_id_list:
            logger.warning("the sample was removed from previous step of the pipeline, ignoring the sample")
        else:
            logger.error("ERROR, the sample directory do not exist, please check the path")
        
        return
    reads_dir=f"{filter_dir}/reads"
    filter_reads_R1, filter_reads_R2 = check_reads(reads_dir, id_sample, logger, output_missing_file)
    if not filter_reads_R1 or not filter_reads_R2:
        return
    #final_assembly dir creation
    try:
        os.makedirs(f"{sample_dir}final_assembly")
    except:
        print(f"Warning : sample directory {sample_dir}final_assembly already existing. Overwriting...")


    #reads_filtering after splashome decontamination:
    filter_suffix=f"{suffix}_filter_reads"
    kr_filter_out=f"{filter_dir}/kraken_results/{id_sample}_{filter_suffix}_kr_results.txt"
    kr_filter_rep=f"{filter_dir}/kraken_reports/{id_sample}_{filter_suffix}_kr_reports.txt"

    #print(kr_filter_rep) #check


    #creation of the contaminants list
    contaminants_list = []
    contam_file=contaminants
    try :
        with open(contam_file, 'r') as bl: 
            lines = bl.readlines()
            for i, line in enumerate(lines) :
                line = line.strip()
                contaminants_list.append(line)
                message = f"taxid {line} has been added as a contaminant\n"
                print(message)
    except:
        print(f"Error : the contaminant file cannot be found \n")
        sys.exit()
    contaminants_str=' '.join(f"{str(e)}" for e in contaminants_list)
    print(f"taxon selected for filtering: {contaminants_str}\nNow running reads decontamination\n")

    #variable initialization
    final_dir=f"{sample_dir}final_assembly"

    #read decontamination
    decontam_reads_R1, decontam_reads_R2 = blacklist_filtering(filter_reads_R1, filter_reads_R2, id_sample, final_dir, kr_filter_out, filter_suffix, contaminants_str, kr_filter_rep, logger, sample_out_dir)
    #taxonimic assignement on the reads, will help to assess reads decontamination, maybe not very usefull in the long run. to discuss.
    decontam_suffix=f"{suffix}_deconta_reads" #est ce utile pour le moment? a verifier
    kr_deconta_out,kr_deconta_reports = launch_kraken(decontam_reads_R1, decontam_reads_R2, id_sample,final_dir,kraken_db,suffix,nb_threads, mem, logger, sample_out_dir)
    #genome assembly variable
    output_dir=f"{decontam_dir}/reads"
    assembly_dir=f"{final_dir}/assembly"

    #run the final assembly on decontaminated reads.
    print(f"running final assembly in {assembly_dir}")
    deconta_assembly=genome_assembly(decontam_reads_R1, decontam_reads_R2, id_sample, suffix, assembly_dir, nb_threads, logger, sample_out_dir)

    #pas utile...?
    cleaned_dir=f"{output}/cleaned_assembly"
    try:
        os.makedirs(cleaned_dir)
    except:
        f"Warning : sample directory {cleaned_dir} already existing. Overwriting..."
    #deconta_assembly=f"{assembly_dir}/{id_sample}_{suffix}_raw_assembly.fasta"

    #def contigs_cleaning(R1_file, R2_file, assembly, id_sample, suffix, output, threads):
    print(f"running assembly cleaning in {cleaned_dir}")
    cleaned_assembly=contigs_cleaning(decontam_reads_R1, decontam_reads_R2, deconta_assembly, contaminants_str, id_sample, suffix, nb_threads, logger, sample_out_dir, mem, sample_dir)

    #ajouter une boucle if pour junk contigs
    if keep_junk == True:
        logger.info("keep discarded option activated.")
        junk_contig_analysis(decontam_reads_R1, decontam_reads_R2, deconta_assembly, contaminants_str, id_sample, suffix, nb_threads, logger, sample_out_dir, mem, sample_dir)

    message = f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} : decontamination script end for sample {id_sample}\n" #print the end time of the script, in order to have easy access to the total time the script ran.
    print(message)

if __name__ == "__main__":
    arguments = handle_program_options()
    list_missing_id = []
    path_input = arguments.input
    id_sample = arguments.id
    output = arguments.output           
    suffix = arguments.name               ### only one output path and suffix as meta-analysis is performed
    kraken_db = arguments.kraken_db
    mem = arguments.mem
    nb_threads = arguments.threads
    contaminants=arguments.contaminants
    keep_junk=arguments.keep_junk_contigs

    #implement path and verification step
    run_dir=f"{output}/"
    sample_dir=f"{run_dir}{id_sample}/"
    message = f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} : read filtering begin for sample {id_sample}\n"
    print(message)
    #folder verification
    taxon_dir= f"{run_dir}meta_analysis/"
    decontam_dir=f"{sample_dir}decontamination"

    #check if the relevant directories exist, in order to stop the script asap if the needed folder do not exist, as this is the 5th part of the workflow 
    if os.path.exists(run_dir):
        print("Now running reads decontamination and final assembly, this might take a while\n")
    else:
        print(f"the directory {run_dir} does not exist, it is intended to be created with the script Sincere-data_reads_filtering.py as it is the first step of this pipeline, please check your path or run the previous step of the pipeline on the samples first before running this script again.")
        sys.exit()

    if os.path.exists(taxon_dir):
        print(f"sucessfully found the taxon directory: {taxon_dir}\n")
    else:
        print(f"WARNING: the directory {taxon_dir} does not exist, it is likely that the script sincere-data_meta_analysis.py have not been run on the samples prior to running this Script. As the meta analysis is required for this script to work properly is an important part of the pipeline you should run it before running this script again.")
        sys.exit()

    if os.path.exists(decontam_dir):
            print(f"sucessfully found the taxon directory: {decontam_dir}\n")
    else:
        print(f"WARNING: the directory {decontam_dir} does not exist, it is likely that the script Sincere-data_outlier_predictor.py have not been run on the samples prior to running this Script. As the meta analysis is required for this script to work properly is an important part of the pipeline you should run it before running this script again.")
        sys.exit()
    input_missing_file=f"{path_input}/missing_id.txt"
    output_missing_file=f"{output}/missing_id.txt"
    try:
        with open(input_missing_file, 'r') as miss_f:
            previous_missing_id_list = [line.strip() for line in miss_f]
            with open(output_missing_file, 'w') as out_miss_f:
                for missing_sample in previous_missing_id_list:
                    out_miss_f.write(missing_sample + '\n')
    except:
        print("there is not missing file yet, the analysis will run on the complete dataset")
        previous_missing_id_list=[]
    print(f"Successfully found {id_sample}, proceding to reads decontamination.\n")
    if os.path.exists(sample_dir):
        print(f"Running the decontamination for {id_sample} in {sample_dir} directory\n")
    else:
        if sample_dir in list_missing_id:
            print("the sample was removed fro previous step of the pipeline, ignoring the sample")
            sys.exits()
        else:
            print("ERROR, the sample directory do not exist, please check the path")
            sys.exit()
    #process_sample(id_sample, path_input, suffix, kraken_db, run_dir, nb_threads, mem, previous_missing_id_list, output_missing_file, contaminant)
    #print("done!")

    # Get the maximum number of logical cores and divide by nb_threads  
    max_workers = max(1, cpu_count() // nb_threads)  # Also ensure at least 1 worker
    print("CPU count: ", cpu_count(), ", running with ", max_workers, " workers.", flush=True)
    
    partial_process_sample = partial(process_sample, 
                                     path_input=path_input,
                                     suffix=suffix, 
                                     run_dir=run_dir,
                                     nb_threads=nb_threads,
                                     mem=mem,
                                     kraken_db=kraken_db,
                                     contaminants=contaminants
                                     previous_missing_id_list=previous_missing_id_list,
                                     output_missing_file=output_missing_file)
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        executor.map(partial_process_sample, sample_ids)
        
    # Remove duplicates in the output missing file
    try: 
        with open(output_missing_file, 'r') as out_miss_f:
            unique_missing_ids = set(out_miss_f)
        with open(output_missing_file, 'w') as out_miss_f:
            for id_sample in unique_missing_ids:
                out_miss_f.write(id_sample)
    except:
        print("The output missing file is empty.")