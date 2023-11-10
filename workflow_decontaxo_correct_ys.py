##!/usr/bin/python

import os,sys
import glob
import time
import matplotlib.pyplot as plt
import math
from datetime import datetime
import numpy as np 
from collections import OrderedDict
import pandas as pd
from pathlib import Path
from Bio import Entrez
from textwrap import dedent as twdd
import argparse
import convert_krep as krep

######################################################################################################################
##################################################### Functions ######################################################

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
    parser.add_argument('-o', '--output', required=True,
                            help="Output path where directories will be built")
    parser.add_argument('-n', '--name', required=True,
                            help="Name of the run and main directory")
    parser.add_argument('-f', '--fasta', default=True,
                        help="Give False as an argument to stop generation \
                            of fasta sequences of outliers found")
    parser.add_argument('-ts', '--thres_samples', type=float, default=0.1,
                        help="Allowed proportion for a taxa to be found in samples \
                            and considered as contaminant when less than the float value")
    parser.add_argument('-to', '--thres_outliers', type=float, default=0.9,
                        help="Minimal proportion for a taxa found in outliers to be \
                            considered as contaminant when upper than the float value")
    parser.add_argument('-b', '--blacklist', required=False ,default=None,
                        help="path to a blacklist file containing the taxon that should be systematically removed from the dataset (one taxid per line). Carefull it will eliminate all the chikdren taxa as well!")
    return parser.parse_args()

def convert_list(py_list):
    """ takes a python list as an input and converts it to a bash list, with the format ("elem1" "elem2" ...)"""
    sh_list=' '.join(f"{str(e)} " for e in py_list)
    return f"\"{sh_list}\""

def check_reads(path_input, id_sample, number, path_output):
    """ from an input path, a sample id, the number of paired reads (either 1 or 2), test the path then check if read files are present 
        and single-copied, and returns the file name if validated. If not, False is returned to cut the decontamination workflow for the sample. """
    message = ""
    if os.path.exists(path_input):
        file_list = glob.glob(path_input + f"{id_sample}*R{number}.fastq*")    # find files corresponding to regex
        if len(file_list) < 1:
            message = f"Error : no file found for {id_sample} R{number}:\nSample ejected from the workflow.\n"
            return False, message 
        elif len(file_list) > 1:
            message = f"Error : ambigous read files for sample {id_sample} R{number}:\n"
            for read_files in file_list:
                message = message + f"{str(read_files)}\n"
            message = message + f"Sample ejected from the workflow.\n"
            return False, message 
        else :
            reads = file_list[0]
            message = message + f"reads R{number} found for sample {id_sample}:\n{reads}\n"
            return reads, message 
    else :
        print ("Invalid input path specified in the index file. Please provide existing pathes.")
        sys.exit()

def remove_sample(i, nb_files):
    """ updates converted bash lists and python lists (global variables) using the index of a sample from python lists to remove """
    list_path_input_sh.replace(list_path_input[i], "")
    list_path_input.pop(i)

    list_ids_sh.replace(list_ids[i], "")
    list_ids.pop(i)

    nb_files -= 1
    return nb_files

def launch_sarray(name, nb_files, max_array, args, script, suffix):
    """
    Input :
        name : name of sarray (needs to be words separated by '_' for generation of log files) 
        nb_files : number of files to adjust sarray
        max_array : number max of arrays to run simultaneously (limited for Kraken, needing a lot of memory)
        args : string of arguments to give to bash scripts
        script : name of bash script
    Generates, prints and launches a sarray command and waits for the end of the job using *wait_sarray* function
    """
    tmp_file=f"tmp_{name}_{suffix}.log"
    if os.path.exists(tmp_file):    # generates a new file
        os.remove(tmp_file)
    parameters=f"--array=0-{nb_files-1}%{max_array} -o {logs_dir}{name}_%a.log -e {logs_dir}{name}_%a.log --job-name={name}_{suffix}"
    print(f"command {name} :\nsbatch {parameters} ./{script} {args} {tmp_file}")
    os.system(f"sbatch {parameters} ./{script} {args} {tmp_file}")
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} : {name.replace('_', ' ')} launched, logs are available at {logs_dir}.\nWaiting for the end of sarray...")
    wait_sarray(tmp_file, nb_files)
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} : {name.replace('_', ' ')} completed, temporary file {tmp_file} destroyed.\n")
    
def wait_sarray(tmp_file, nb_files):
    """ all sarray write the number of the iteration in a tmp file once complete. This function checks every minute if a file is generated, 
        once confirmed it checks if the number of iterations is reached, and breaks if it's the case """
    while not os.path.exists(tmp_file):
        time.sleep(30)
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} : Job has started")
    lines=[]
    while len(lines) < nb_files:
        with open (tmp_file,"r") as tmp:
            lines = tmp.readlines()
        time.sleep(30)
    os.remove(tmp_file)
    return 

def get_avg(list):
    return(round(sum(list) / len(list),2))
    """ à remplacer par la fonction de math qui donne la moyenne d'une liste"""

def get_visu_cov(cover, threshold, fig_name, id_file):
    """ Draw and save a jpg file named [[fig_name]] that shows values of coverage accross the genome and the value of the threshold
        defined depending on z-score (2) ; returns True if successful """
    try :
        fig, ax = plt.subplots()
        ax.plot(cover, label='cover')
        ax.axhline(y=threshold, color='r', label=f'threshold : {threshold}')
        ax.set_ylabel('coverage')
        ax.set_xlabel(f'position on assembly')
        ax.set_title(f"Visualization of purif outliers on {id_file}, berudes run")
        ax.legend()
        # print(f"Param : \ncover:\n{cover}\nthreshold : {threshold}, id_file : {id_file}\nfig_name :{fig_name}\n\n")
        ## adjust graph size
        fig.set_figwidth(fig.get_figwidth() * 2.5)
        fig.set_size_inches(fig.get_size_inches() * 1.2)
        plt.savefig(fig_name)
        plt.close()
        return True
    except :
        return False

def write_bedfile(name_bedfile, name_bedgraph, threshold):
    """
    Input : 
        name_bedfile : String, name of the file to generate
        name_bedgraph : String, name of the file to read
        threshold : Int, value of the threshold to write the bedfile
    Performs a bedgraph analysis using a theshold value to determine the outliers of a genome. Outliers 
    are merged when too close to another one ('gap' value).
    Output : 
        none, but bedfile is written
    """
    gap = 10
    min_length = 30
    start = None
    nb_fasta=0       # number of regions with high cover depth
    passed=False     # by default, passed is false when the threshold hasn't been reached
    with open(name_bedfile, 'w') as bedfile: 
        with open(name_bedgraph, 'r') as bedgraph: 
            nb_fasta=0       # number of regions with high cover depth
            lines = bedgraph.readlines()
            passed=False        # by default, passed is false when the threshold hasn't been reached
            waiting = False     # by default, there's no zone saved in memory waiting for a neighboring check

            for line in lines :
                line = line.strip()
                line = line.split("\t")

                #print(int(line[3]))
                if (int(line[3])>threshold and not passed):         # if a new exceeding region is starting...
                    if waiting:                                         # if the previous regions is waiting to merge... 
                        if (int(line[1]) - int(last_end) < gap and line[0]==last_contig) :    # and if one finished less than the authorized gap before...
                            passed = True                           # ...threshold is count as exceeded,...
                            start=last_start                        # ...previous zone's beginning is used, i.e they're merged...
                            last_start,last_end,last_contig=0,0,""  # ...last values are reset (is it enough ?)...
                            waiting = False                         # ...and there's no longer a zone waiting for check    
                        else :                          # else (i.e it's really a new exceeding region)
                            if (int(last_end) - int(last_start) > min_length):    # ... if the outlier is long enough...
                                nb_fasta+=1                                 # ...the previous one is written...
                                bedfile.write(f"{last_contig}\t{last_start}\t{last_end}\tregion_{nb_fasta}\n")
                            waiting = False                     # ...there's no longer a zone waiting for check
                            passed = True                       # ...passed is set as True,...
                            the_contig=line[0]                  # ...contig's name is saved,...
                            start=line[1]                       # ...and starting position is saved
                    else :  
                        passed = True                   # ...passed is set as True,...
                        the_contig=line[0]              # ...contig's name is saved,...
                        start=line[1]                   # ... and starting position is saved
                        
                if (int(line[3])>threshold and line[0]==the_contig and passed):   # if cover continues to exceed the threshold on the same contig...
                    end=line[2]                                             # ...ending position is updated
                    #print(f"continue :  {line} : {start}-{end}")

                elif(int(line[3])>threshold and line[0]!=the_contig and passed) : # if cover continues to exceed the threshold on another contig...
                    passed = False                                          # ...region is ended,...
                    if (int(end) - int(start) > min_length):    # ... if the outlier is long enough...
                        nb_fasta+=1                                             # ...counted as a new one...
                        bedfile.write(f"{the_contig}\t{start}\t{end}\tregion_{nb_fasta}\n") # ...and added to the file

                if (int(line[3])<=threshold and passed):  # if cover ends exceeding the threshold...
                    passed = False              # ...region is ended...
                    last_start=start                ##
                    last_end=end                     #
                    last_contig=the_contig           # 
                    waiting=True                     ## .... and saved and marked as waiting for next test 
                    #print(f"finished :  {line}\n")

            if waiting :    
                if (int(last_end) - int(last_start) > min_length):    # ... if the outlier is long enough...    
                    nb_fasta += 1                                                                   ##
                    bedfile.write(f"{last_contig}\t{last_start}\t{last_end}\tregion_{nb_fasta}\n")   ## if last outlier is saved in memory, print it

def count_taxa_outliers(list_outliers):
    """ From a list of dictionaries (1 per sample, key=taxa : value=occurrency), counts occurrencies of taxa 
        accross all samples. Returns a dictionary of taxas and their occurrency accross samples """
    taxas_found = OrderedDict()             # reset dictionnary 
    for dico_outliers in list_outliers:
        for taxa in dico_outliers.keys():
            if taxa in taxas_found.keys():
                taxas_found[taxa] += 1
            else :
                taxas_found[taxa] = 1
    return taxas_found

def taxName2taxID(list_names):
    """ as outliers' main taxas may be linked to samples' main taxas through going up in the taxonomy, taxID needs to be recovered using 
        Entrez module. 
        Input : 
            list of taxonomic names  
        Output : 
            dictionnary {taxName : taxID}, list of taxNames which failed their identification
    """
    Entrez.email = 'yann.sevellec@univ-rennes1.fr'
    dico_convert = {}
    error = []
    for i, taxName in enumerate(list_names) :
        if (taxName != ""):
            handle = Entrez.esearch(db='Taxonomy', term=taxName)
            record = Entrez.read(handle)
            print(f"{i} out of {len(list_names)} taxas to convert")
            
            if record['Count'] == '1':
                taxID = record['IdList'][0]
                dico_convert[taxName] = taxID
            elif record['Count'] == '0':
                error.append(taxName)
                dico_convert[taxName] = "Failed"
            else:
                error.append(taxName+"*")
                taxID = record['IdList']
                taxID = convert_list(taxID)
                dico_convert[taxName] = taxID
    
    return dico_convert, error

def meta_analysis(taxas_outliers, df_main_samples, dico_convert, suffix):
    """ 
    Input : 
        taxas_outliers : Ordered dictionnary of taxa names (keys) and their occurrency in samples' outliers (values)
        df_main_samples : dataframe of stats on main taxas from samples
        dico_convert : dictionnary of taxa names (keys) and their affiliated taxID (values)
    Performs meta-analysis on all samples to return a
    Output :
    df_meta, a 6-column dataframe : 
        "taxID"                     :   taxonomy ID of taxa
        "taxa"                      :   taxonomy name
        "occurrency_in_outliers"    :   number of samples having this taxa in their outliers
        "ratio_outliers"            :   occurrency_in_outliers divided by the number of samples
        "occurrency_in_samples"     :   number of samples having this taxa as main taxa from reads
        "ratio_samples"             :   occurrency_in_samples divided by the number of samples
    """
    df_meta = pd.DataFrame(columns=["taxID","taxName","occurrency_in_outliers","ratio_outliers","occurrency_in_samples","ratio_samples"])
    
    ### sort taxas from outliers by their occurrency  
    taxas_by_occurrency = OrderedDict(sorted(taxas_outliers.items(), key=lambda x: x[1], reverse=True))
    for taxa in taxas_by_occurrency.keys():

        ### get ratios from outliers
        ratio_outliers = taxas_by_occurrency[taxa]/nb_files

        ### get ratios from samples
        occurrency_in_main = df_main_samples['taxo'].apply(lambda x: taxa in x.values()).sum()
        ratio_samples = occurrency_in_main/len(df_main_samples)

        ### save results 
        line = [dico_convert[taxa],taxa,taxas_by_occurrency[taxa],ratio_outliers,occurrency_in_main,ratio_samples]
        df_meta.loc[len(df_meta)]=line

    return df_meta

######################################################################################################################
####################################################### Code #########################################################

### general parameters 
        
task = 10              # task authorized for standard sbatches 
task_for_kraken = 5    # task authorized for kraken sbatches (require 300Gb of memory each)

global list_path_input, list_ids
global list_path_input_sh, list_ids_sh
global nb_files

### assessing index file and others by getting the arguments 
arguments = handle_program_options()        

### assimilate them
path_index_file = arguments.input
output = arguments.output           ###
suffix = arguments.name               ### only one output path and suffix as meta-analysis is performed
get_Fasta = arguments.fasta         # option to generate an optional file of fasta sequences found in peaks (junk)
thres_outliers = arguments.thres_outliers       
thres_samples = arguments.thres_samples
path_blacklist = arguments.blacklist

list_path_input, list_ids = [], []

### read the index file and extract pathes for input and names of samples 
try :
    with open(path_index_file, 'r') as f: 
        lines = f.readlines()
        for i, line in enumerate(lines) :
            line = line.strip()
            line = line.split(",")
            try : 
                list_path_input.append(line[0])
                list_ids.append(line[1])
            except:
                print(f"Error : file does not contain 2 columns at line {i}, please make sure it's comma-separated and that the last line is filled")
                sys.exit()
except :
    print(f"Error : index file not found at {path_index_file}")
    sys.exit()

### checking the lecture of index file 
if list_path_input == [] or list_ids == []:
    print(f"Error : index file {path_index_file} is empty")
    sys.exit()


### preparing data for bash scripts
list_path_input_sh=convert_list(list_path_input)
list_ids_sh=convert_list(list_ids)

nb_files=len(list_ids)
print(f"\nSuccessfully found {nb_files} valid lines in index file. Processing directory building...\n")

################### Preparing environment : building output directories and checking the input #######################

### building run directory
run_dir=f"{output}{suffix}/"
try :
    os.makedirs(run_dir)
except :
    print(f"Warning : run directory {run_dir} already existing. Overwriting...")

### redirect all prints to a general log file
general_log = f"{run_dir}general_log_{suffix}.log"
print(f"All prints will now be saved in {general_log} file.")
if os.path.exists(general_log):    # generates a new file
    os.remove(general_log)
sys.stdout = open(general_log, 'a')

### same for errors
err_log = f"{run_dir}errors_log_{suffix}.log"
print(f"{err_log} file generated for possible errors.")
if os.path.exists(err_log):  
    os.remove(err_log)
sys.stderr = open(err_log, 'a')

i = 1
while i < nb_files:     # 'while' perform better than 'for' when in need to remove a sample from the set 
    
    ### building directories
    sample=list_ids[i]
    print(f"Sample {i} : {sample}...")
    sample_dir=f"{run_dir}{sample}"
    try :
        os.makedirs(sample_dir)
    except :
        print(f"Warning : main directory {sample_dir} already existing. Overwriting...")
    
    ## starting log file
    with open(f"{sample_dir}/{sample}.log","w") as log : 
        log.write(f"Welcome on log file of sample {i} of this run : {sample}\n\n")
        try :
            os.makedirs(f"{sample_dir}/decontamination")
            os.makedirs(f"{sample_dir}/final_assembly")
            print(f"subdirectories and log generated for {sample} available in :\n{sample_dir}\n")
        except :
            print(f"Warning : directories already existing for {sample}. Overwriting...\n")

        ### checking the reads 
        path_input=list_path_input[i]
        id_sample=list_ids[i]

        ### checking the arguments in input file and removing incomplete samples
        res_read1 = check_reads(path_input,id_sample,1,output)
        res_read2 = check_reads(path_input,id_sample,2,output)
        log.write(str(res_read1[1]))
        log.write(str(res_read2[1]))
        if(type(res_read1[0]) != str or type(res_read2[0]) != str):   
            nb_files = remove_sample(i, nb_files)
            print(f"ERROR : something went wrong while checking the read files for {sample}. \nPlease check {sample_dir}/{sample}.log for more information.\n")
            continue
    i += 1 

### create logs directory if not already existing
logs_dir=f"{run_dir}logs/"
try :
    os.makedirs(logs_dir)
except :
    pass

### prepare_decontaxo : assembly, alignment and bam generation
args = f"{list_path_input_sh} {list_ids_sh} {suffix} {output} "
launch_sarray("prepare_decontaxo", nb_files, task, args, "prepare_decontaxo.sh", suffix)  

### kraken_on_reads : get kraken reports of trimmed reads
launch_sarray("kraken_on_reads", nb_files, task_for_kraken, args, "kraken_on_reads.sh", suffix)      # let's not reach 3TB of memory for kraken sarrays (arrays limited to 5)  

### kraken_on_assembly : get kraken reports of contigs on trimmed 
launch_sarray("kraken_on_assembly", nb_files, task_for_kraken, args, "kraken_on_assembly.sh", suffix)     

## visu_cov and get_peaks : 
# df for overall stats 
df_summary = pd.DataFrame(columns=['file', 'mean', 'max', 'standard error', 'length', 'threshold'])
index = 0
while index < len(list_ids):        # needs a while instead of a for to reduce index if something failed

    print(f"\n{index+1} out of {len(list_ids)}")
    id_file=list_ids[index] 

    ### To define the threshold
    sample_dir=f"{run_dir}{id_file}"
    with open(f"{sample_dir}/{id_file}.log", "a") as log :   
        try :
            with open(f"{sample_dir}/decontamination/{id_file}_{suffix}_bedgraph.tsv", 'r') as f:
                contigs, pos, cover = [], [], []
                line = f.readline()
                line_nb = 1
                while line != "" :
                    line = line.split("\t")
                    contigs.append(line[0])
                    try:
                        begin=int(line[1])
                        end=int(line[2])
                    except:
                        log.write(f"position values not recognized in line {line_nb} :\n{line}\naborting\n")
                        print(f"position values not recognized in line {line_nb} :\n{line}\naborting\n")
                        nb_files = remove_sample(index, nb_files)
                        continue
                    while begin < end :
                        try : 
                            cover.append(int(line[3]))
                            begin+=1
                        except :
                            log.write(f"coverage column not found, or invalid data, at line {line_nb} :\n{line}\naborting\n")
                            print(f"position values not recognized in line {line_nb} :\n{line}\naborting\n")
                            nb_files = remove_sample(index, nb_files)
                            continue
                    line = f.readline()
                    line_nb += 1
        except :
            log.write(f"failed : bedgraph file could not be opened for {id_file}, check {logs_dir}align_{index}.log if bedgraph generation failed\n")  
            print(f"failed : bedgraph file could not be opened for {id_file}, check {logs_dir}align_{index}.log if bedgraph generation failed\n")
            nb_files = remove_sample(index, nb_files)
            continue

        try : 
            maxi=max(cover)     
        except :                # if file is empty....
            row=[id_file, None, None, None, None, None]
            df_summary.loc[len(df_summary)]=row
            message = f"\nfailed : {id_file}_{suffix}_bedgraph.tsv is empty"
            log.write(message+"\n") 
            print(message)
            nb_files = remove_sample(index, nb_files)
            continue            # ...abort file filling
        
        ### get values to define the threshold
        mean = get_avg(cover)
        std = np.std(cover)
        threshold = int(mean + 2*std)

        ## add values to summary file
        row=[id_file, mean, maxi, round(std,2), len(cover), threshold]
        df_summary.loc[len(df_summary)]=row

        ## draw coverage graph and save it
        fig_name=f"{sample_dir}/decontamination/{id_file}_{suffix}_visu_cov.jpg"
        if (get_visu_cov(cover, threshold, fig_name, id_file) == True):
            print (f"Visualization available at {fig_name}")
        else : 
            print (f"{fig_name} generation failed")

        ### writing bedfile
        bedfile = f"{sample_dir}/decontamination/{id_file}_{suffix}.bed"
        bedgraph = f"{sample_dir}/decontamination/{id_file}_{suffix}_bedgraph.tsv"
        write_bedfile(bedfile, bedgraph, threshold)

        ## writing and displaying results
        if get_Fasta : ### generates fasta files [option] 
            fasta_file=f"{sample_dir}/decontamination/{id_file}_{suffix}_peaks.fasta"
            os.system(f"bedtools getfasta -fi {sample_dir}/decontamination/{id_file}_{suffix}_ecrem_contigs.fasta -bed {bedfile} -name > {fasta_file}")
            message = f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} : {bedfile} and {fasta_file} generated\n"
        else :
            message = f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} : {bedfile} generated\n"
        log.write(message)
        print(message)

        index +=1

### print overall stats 
df_summary.to_csv(f"{logs_dir}summary_stats_{suffix}.tsv", sep="\t")

## peaks_to_fastq 
launch_sarray("peaks_to_fastq", nb_files, task, args, "peaks_to_fastq.sh", suffix)  

## kraken_on_peaks
args = f"{list_path_input_sh} {list_ids_sh} {suffix} {output} "
launch_sarray(f"kraken_on_peaks", nb_files, task_for_kraken, args, "kraken_on_peaks.sh", suffix)

## get main taxas
df_main_outliers = pd.DataFrame(columns=['sampleID','taxID', 'taxName', 'reads_nb', 'taxo']) 
df_main_outliers_unique_strains = pd.DataFrame(columns=['sampleID','taxID', 'taxName', 'reads_nb', 'taxo'])    # dataframe havins one example of each strain found as main in each outlier
df_main_samples = pd.DataFrame(columns=['sampleID','taxID', 'taxName', 'reads_nb', 'taxo'])
#creation outliers file
list_outliers = []
list_global_outliers=[]
outliers_dir =f"{output}{suffix}/outliers/"
try :
    os.makedirs(outliers_dir)
except :
    print(f"The run directory {outliers_dir} already exist")

for SAGs in list_ids:
    list_outliers=[]
    empty_reports=[]
    outlier_id=SAGs
    sample_dir=f"{run_dir}{SAGs}"
    reports_dir = f"{sample_dir}/decontamination/kraken_reports/"
    krep_csv =f"{sample_dir}/{SAGs}_converted_kraken_report.csv"
    krep_sample = f"{reports_dir}{SAGs}_{suffix}_trimmed_report.txt"
    df_taxas_main_sample = krep.convert_report(krep_sample)
    main_taxa_sample = krep.extract_main_taxa(df_taxas_main_sample)
    main_taxa_sample['taxo'] = krep.convert_tax(main_taxa_sample['taxo'])
    main_taxa_sample['sampleID']= SAGs
    df_main_samples.loc[len(df_main_samples)] = main_taxa_sample
    df_main_outliers = pd.DataFrame(columns=['sampleID','taxID', 'taxName', 'reads_nb', 'taxo'])
    valid_reports = len(glob.glob(f"{reports_dir}*junk*_report.txt"))
    for outlier in Path(reports_dir).glob('*junk*_report.txt') :
        krep_name_full = os.path.basename(outlier)             ## get full name of outlier's report 
        krep_name = krep_name_full.replace(suffix + "_","")
        krep_name = krep_name.replace("_report.txt","")
        df_taxas_outliers = krep.convert_report(outlier)
        if len(df_taxas_outliers) == 0:                             ##
            valid_reports -= 1                                       #
            empty_reports.append(krep_name_full)                     #
            continue                                                 ## if the dataframe is empty, its name is saved and next report is processed   
        else : 
            main_taxa_outlier = krep.extract_main_taxa(df_taxas_outliers)               ## ectraction of main taxa from outliers
            main_taxa_outlier['taxo'] = krep.convert_tax(main_taxa_outlier['taxo'])     ## conversion of taxa's taxonomy   
            main_taxa_outlier['sampleID']=SAGs
            df_main_outliers.loc[krep_name] = main_taxa_outlier
            #print(f"The main taxon for {krep_name} is {main_taxa_outlier['taxName']}")
            #print(df_main_outliers['sampleID']+ " "+ df_main_outliers["taxName"])
            if str(main_taxa_outlier['taxID']) not in list_outliers :    ## and adding it to the ddf counting taxas found in the sample if first time encountere
                list_outliers.append(str(main_taxa_outlier['taxID']))
                df_main_outliers_unique_strains.loc[krep_name] = main_taxa_outlier
                #print(f"added {df_main_outliers['taxName']} to the unique id list from {df_main_outliers['sampleID']}")
    taxas_outliers, nb_others = krep.sample_vs_outliers(main_taxa_sample, df_main_outliers)
    total = len(taxas_outliers)
    message = f"\nAmong {total} outliers, {nb_others} were found having a main taxa different than the one found in trimmed reads\n"
    message = message + f"Taxa \t\t occurrency_in_outliers\n"
    # print(f"For {id_sample}, {message}")
    #log.write(message+"\n")
    #for name in taxas_outliers.keys():
            #log.write(f"{name} \t\t {taxas_outliers[name]}\n")
            # print(f"{name} \t\t {taxas_outliers[name]}")
    list_global_outliers.append(taxas_outliers)
    #print(f"The main taxon for {main_taxa_sample['sampleID']} is {main_taxa_sample['taxName']}")
    #print(taxas_outliers)
    df_main_outliers.to_csv(f"{sample_dir}/{SAGs}_outliers.csv", index=True)
    if len(empty_reports) != 0:
        nb_outliers=len(glob.glob(f"{reports_dir}*junk*_report.txt"))
        message = f"{valid_reports} reports out of {nb_outliers} were not empty or not filled with unassigned reads.\n"
        message += "Please check :\n" + "\n".join(rep for rep in empty_reports)
        message += f"\nFrom {reports_dir} if any doubt"
        print(f"For {id_sample}, {message}")
        with open(f"{sample_dir}/{id_file}.log", "a") as log :
            log.write(message+"\n")
#print(df_main_outliers_unique_strains)
df_main_samples.to_csv(f"{output}Main_taxons_samples_{suffix}.csv", header=True)


#parse main taxon for anomalies
occurrency_in_main = df_main_samples['taxName'].value_counts()
if occurrency_in_main.iloc[0] > 0.7 * len(list_ids):     # if a taxa if found as main in more than 90% of samples... 
    conta_name = occurrency_in_main.index[0]
    dicoID, error = taxName2taxID([conta_name])
    conta_ID = dicoID[conta_name]
    message = f"WARNING : {conta_name} taxa of taxID {conta_ID} is found as main in more than 70% of the samples.\n"
    message += "Please consider decontaminating this taxa before running this pipeline\nTaxa occurrency in main :\n"
    #print(message,occurrency_in_main)
    sys.exit()
nb_files=len(list_ids)
df_main_outliers_unique_strains.to_csv(f"{output}Main_taxons_outliers_{suffix}.csv", header=True)
taxas_found = count_taxa_outliers(list_global_outliers)

#print(list_global_outliers)
#print(taxas_found)
#print(nb_files)
log = open(general_log, 'a')
dico_convert, err = taxName2taxID(taxas_found.keys())
if err :
    print(f"WARNING : conversion(s) of taxName failed for the following taxas :")
    print(" ".join(e for e in err))
    print("taxas with a * have several taxas IDs referenced that will be used for decontamination. Others don't have any ID proposed")
    # print(dico_convert)
#creation du rapport
df_meta = meta_analysis(taxas_found, df_main_samples, dico_convert, suffix)
meta_file = f"{output}results_meta_analysis_{suffix}.csv"
df_meta.to_csv(meta_file, header=True)
log.write(f"Results of meta-analysis available in file : \n{meta_file}")

### extracting contaminants into a list 
contaminants = []
for index, row in df_meta.iterrows():
    if row['ratio_outliers']>thres_outliers and row['ratio_samples']<thres_samples:
        contaminants.append(row['taxID'])
        message = f"{row['taxName']} of ID {row['taxID']} has been added as contaminant"
        log.write(message)

if len(contaminants) == 0:
    message = f"ERROR : no contaminant having an occurrency >{thres_outliers} in outliers and <{thres_samples} in samples in this run.\n"
    message = message + f"Please consider change settings of values 'thres_outliers' and 'thres_samples' while launching the code \n"
    log.write(message)
    sys.exit()
if path_blacklist !=None:
    try :
     with open(path_blacklist, 'r') as f: 
         lines = f.readlines()
         for i, line in enumerate(lines) :
             line = line.strip()
             contaminants.append(line)
    except:
        print(f"Error : the blacklist file cannot be found")
        sys.exit()

### convert the list and run krakentools to remove reads assigned to contaminants
contaminants_sh = convert_list(contaminants)
args_krakentools = f"{list_path_input_sh} {list_ids_sh} {suffix} {output} {contaminants_sh}"
launch_sarray(f"krakentools", nb_files, task, args_krakentools, "krakentools.sh", suffix)  
log.write("decontamination ended")
### remove krakentools logs as they take A LOT of space
for krakentool_log in Path(logs_dir).glob('krakentools*.log') : 
    os.remove(krakentool_log)

# ### assemble, ecremate and fastQC/multiQC 
#launch_sarray("final_assembly", nb_files, task, args, "final_assembly.sh", suffix) 

### get kraken reports of purified reads and assemblies
#args_final = f"{run_dir} {list_ids_sh} {suffix} {output}"
#launch_sarray("kraken_final", nb_files, task_for_kraken, args_final, "kraken_final.sh", suffix) 

### get overall stats on assemblies
#launch_sarray("stats_final", nb_files, task, args_final, "stats_final.sh", suffix) 

print(f"workflow is over and you can find assemblies and stats of each sample in {run_dir}analysis_on_final/")
log.close()
sys.exit()