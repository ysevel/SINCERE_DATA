##!/usr/bin/python
import os,sys
import glob
import time
from datetime import datetime
from pathlib import Path
from textwrap import dedent as twdd
#from Bio import SeqIO
import argparse
import subprocess
import logging

#functions
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
    parser.add_argument('-n', '--name', required=True,
                            help="Name of the run and main directory")
    parser.add_argument('-k', '--kraken_db', required=True,
                        help="path to the kraken database folder")
    parser.add_argument('-b', '--blacklist', required=False ,default=None,
                        help="path to a blacklist file containing the taxon that should be systematically removed from the dataset (one taxid per line). Carefull it will eliminate all the children taxa as well!")
    parser.add_argument('-t', '--threads', type=int ,default=1,
                        help="Number of threads")
    parser.add_argument('-m', '--mem', action='store_true',
                        help="allow kraken to run on memory (optionnal, otherwise will run kraken2 with the --memory-mapping option that avoid loading kraken databese into the memory but require significantly higher computation time. If you want to load kraken2 database into the memory best to make sure you allocated enough memory the PFPplus kraken2 database require 300G and the nt database 750G to run properly for example)")
    return parser.parse_args()

def remove_sample(id_sample,error):
    """ report the sample id of a missing sample and write the id and the corresponding error message to the missing file."""
    with open(missing_file, 'a') as missf:
        missf.write(id_sample + '\t' + error)
        missing_id.append(id_sample)
        print(id_sample)


def check_reads(path_input, id_sample):
    """ from an input path, a sample id, the number of paired reads (either 1 or 2), test the path then check if read files are present 
        and returns the file name if validated. If not, False is returned to cut the decontamination workflow for the sample. """
    message = "" #note: it is best that the sample id is separated from any suffix in the reads by '_', otherwise finding the right read can be tricky.
    input_R1=f"{path_input}/{id_sample}_*R1.fastq*"
    input_R2=f"{path_input}/{id_sample}_*R2.fastq*"
    R1_list=glob.glob(input_R1)
    print(R1_list)
    R2_list=glob.glob(input_R2)
    if os.path.exists(path_input):
        print(glob.glob(path_input))
        if len(R1_list) < 1:
            message = message + f"Error : R1 reads file not found for {id_sample} \n"
            remove_sample(id_sample,message)
            print(message)
            sys.exit()
        elif len(R2_list) < 1:
            message = message + f"Error : R2 reads file not found for {id_sample} \n"
            remove_sample(id_sample,message)
            print(message)
            sys.exit()
        else :
            R1_file=R1_list[0]
            R2_file=R2_list[0]
            message = message + f"both reads file found for sample {id_sample}\n"
            print(message)
            return R1_file, R2_file
    else :
        print ("Invalid input path specified in the index file. Please check if the path truly exist.")
        sys.exit()

def launch_kraken(R1_file, R2_file, id_sample,output,kraken_db, suffix, threads):
    """
    create the result directories and launch the krakern command on paired reads
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
    if mem == True:
        cmd_mem="" 
    else:
        cmd_mem="--memory-mapping" #allow to reduce memory usage if none specified, warning, kraken become painfully low, use smaller databases.
    command_kr=f"kraken2 --threads {threads} --db {kraken_db} --paired {R1_file} {R2_file} -output {kr_out} --use-names --report {kr_rep} {cmd_mem}"
    print("kraken will run on reads using :\n" + command_kr)
    log_kr=f"{sample_dir}{id_sample}_{suffix}_kraken2.log"
    with open(log_kr,'w')as log_k:
        subprocess.run(command_kr, stdout=log_k,stderr=subprocess.STDOUT, shell=True)
    return kr_out, kr_rep

def blacklist_filtering(R1_file, R2_file, id_sample,output,kr_filter_out, suffix, contaminants, kr_filter_reports):
    """ launche reads decontamination from an list of contaminants using krakentools."""
    output_dir=f"{output}/reads" #try to create the folder for the filtered reads
    try:                                
        os.makedirs(f"{output_dir}")
    except:
        f"Warning : sample directory {output_dir}/ already existing. Overwriting..."

    filter_contam=f"python /scratch/ysevellec/script/KrakenTools-master/extract_kraken_reads.py -1 {R1_file} -2 {R2_file} -k {kr_filter_out} -o {output_dir}/{id_sample}_{suffix}_R1.fastq -o2 {output_dir}/{id_sample}_{suffix}_R2.fastq -t {contaminants} --include-children --exclude --report {kr_filter_reports} --fastq-output"
    print("read filtering command :\n" + filter_contam) #specific log for the reads filtering operation
    log_filter=f"{sample_dir}{id_sample}_{suffix}_filter_reads.log"
    with open(log_filter,'w')as log_f:
        subprocess.run(filter_contam, stdout=log_f,stderr=subprocess.STDOUT, shell=True)
    return f"{output_dir}/{id_sample}_{suffix}_R1.fastq",f"{output_dir}/{id_sample}_{suffix}_R2.fastq"
### general parameters 
arguments = handle_program_options()

path_input = arguments.input
id_sample = arguments.id
output = arguments.output           
suffix = arguments.name               ### only one output path and suffix as meta-analysis is performed on the whole dataset, you should preserve it for all the steps of this workflow.
path_blacklist = arguments.blacklist
kraken_db = arguments.kraken_db
mem = arguments.mem
nb_thread = arguments.threads

### building run & sample directory
run_dir=f"{output}/"
try :
    os.makedirs(run_dir)
except :
    print(f"Warning : run directory {run_dir} already exist. Overwriting...")

sample_dir=f"{run_dir}{id_sample}/"
try :
    os.makedirs(sample_dir)
except :
    print(f"Warning : sample directory {run_dir} already exist. Overwriting...")

#redirect all prints to a general log file
general_log = f"{sample_dir}reads_filtering_{id_sample}_{suffix}.log"
print(f"All prints will now be saved in {general_log} file.")
if os.path.exists(general_log):    # generates a new file
    os.remove(general_log)
sys.stdout = open(general_log, 'a')

#same for errors
err_log = f"{sample_dir}reads_filtering_{id_sample}_{suffix}.err"
print(f"{err_log} file generated for possible errors.")
if os.path.exists(err_log):  
    os.remove(err_log)
sys.stderr = open(err_log, 'a')

#create file the subfolders
print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} : read filtering start")
try :
    os.makedirs(f"{sample_dir}trimmed_reads")
    os.makedirs(f"{sample_dir}filtered_reads")
    os.makedirs(f"{sample_dir}decontamination")
    os.makedirs(f"{sample_dir}final_assembly")
except :
    print(f"Warning : directories already existing for {id_sample}. Overwriting...\n")

#reinitialization of the missing file list
missing_file = f"{run_dir}missing_id.txt"
print(f"The missing sample will be written in {missing_file}.")
if os.path.exists(missing_file):    # generates a new file
    os.remove(missing_file)

missing_id=[]

#create the missing sample file:
reads_R1,reads_R2 = check_reads(path_input,id_sample)
print(reads_R1)
print(reads_R2)

#analyse the trimmed reads
print("Now running preliminary taxon classification on reads\n")
trim_suffix=f"{suffix}_trim_reads"
kr_trim_out,kr_trim_reports = launch_kraken(reads_R1, reads_R2, id_sample,f"{sample_dir}trimmed_reads",kraken_db,suffix,nb_thread)
#check the blacklist and create a contaminants list to filter out with krakentools
contaminants = []
if path_blacklist !=None:
    print(path_blacklist)
    try :
     with open(path_blacklist, 'r') as bl: 
         lines = bl.readlines()
         for i, line in enumerate(lines) :
             line = line.strip()
             contaminants.append(line)
             message = f"taxid {line} has been added as a contaminant from the blacklist\n"
             print(message)
    except:
        print(f"Error : the blacklist file cannot be found \n") #issue an error if the path to the blacklist is incorrect
        sys.exit()

output_filter=f"{sample_dir}filtered_reads"
filter_suffix=f"{suffix}_filter_reads"
contaminants_str=' '.join(f"{str(e)}" for e in contaminants)
print(f"taxon selected for filtering: {contaminants_str} from the blacklist \nNow running reads filtering\n") #checck on the blacklist
filter_reads_R1, filter_reads_R2 = blacklist_filtering(reads_R1, reads_R2, id_sample, output_filter, kr_trim_out, filter_suffix, contaminants_str, kr_trim_reports) #filter out the contaminated reads

kr_filter_out,kr_filter_reports = launch_kraken(filter_reads_R1, filter_reads_R2, id_sample,output_filter,kraken_db,filter_suffix,nb_thread) #launch kraken on the filtered reads in order to use the reports in the following steps of the workflow
print(f"script successfully ended check {filter_reads_R1} and {kr_filter_reports} now closing and waiting for main taxon meta analysis")
print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} : read filtering ends")