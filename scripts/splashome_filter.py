##!/usr/bin/python
import os,sys
import glob
import time
from datetime import datetime
from textwrap import dedent as twdd
from Bio import SeqIO
import pandas as pd
import argparse
import convert_krep as krep
import subprocess

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
                        help="path to index file of samples to process.")
    parser.add_argument('-o', '--output', required=True,
                        help="Output path where directories will be built.")
    parser.add_argument('-n', '--name', required=True,
                        help="suffix added to each output file.")
    parser.add_argument('-b', '--blacklist', required=False ,default=None,
                        help="path to a blacklist file containing the taxon that should be systematically removed from the dataset (one taxid per line). Carefull it will eliminate all the children taxa as well!")
    parser.add_argument('-ta', '--thres_anomaly', type=float ,default=1,
                        help="must be a value between 0 and 1, indicate the proportion of the collection allowed to present the same main taxa (default= 1: for no alert for overrepresented main taxa)")
    return parser.parse_args()

def remove_sample(id_sample,error):
    """ report the sample id of a missing sample and write the id and the corresponding error message to the missing file."""
    with open(missing_file, 'a') as missf:
        missf.write(id_sample + '\t' + error)
        missing_id.append(id_sample)
        print(id_sample)

def check_missing_samples(path_index_file,path_missing_file):
    """ check for a sample if a sample have already have been reported as missing in a previous step of the pipeline. stop the script directly if the sample is already flagged as missing."""
    try :
        with open(path_index_file, 'r') as id_f: #retrieve sample id from a text file with a list, one id per line
            lines = id_f.readlines()
            for line in lines:
                line = line.strip()
                list_id.append(line)
    except :
        print(f"Error : index file not found at {path_index_file}") #if the file does not exist, stop the script
        sys.exit()
    try: 
        with open(path_missing_file,'r') as miss_f: #chek if a missing file exit and creata a list of the file found missing at previous step of the pipeline
            miss_lines = miss_f.readlines()
            for line in miss_lines:
                line=line.strip()
                list_missing_id.append(line)
    except:
        print("there is not missing file yet, the analysis will run on the complete dataset")
        return
    for sample_id in list_id:
        if sample_id in list_missing_id: #exclude sample id found in the missing list from the sample list.
            list_id.remove(sample_id)
            print(f"{sample_id} has been found missing in previos steps and will been removed from the samples.")

#variable initiation
arguments = handle_program_options()
list_id=[]
list_missing_id = []
path_index_file = arguments.input
output = arguments.output           
suffix = arguments.name               ### only one output path and suffix as meta-analysis is performed
path_blacklist = arguments.blacklist
main_abnormal_thresold= arguments.thres_anomaly

run_dir=run_dir=f"{output}/"

###check missing sample
missing_file=f"{run_dir}missing_id.txt"
check_missing_samples(path_index_file, missing_file)
nb_files=len(list_id)
print(f"\nSuccessfully found {nb_files} valid lines in index file.\n")

### checking the lecture of index file 
if list_id == []:
    print(f"Error : index file {path_index_file} is empty")
    sys.exit()

df_main_samples = pd.DataFrame(columns=['sampleID','taxID', 'taxName', 'reads_nb', 'taxo'])
#creation outliers file
taxon_dir= run_dir + "meta_analysis/"
try :
    os.makedirs(taxon_dir)
except :
    print(f"The run directory {taxon_dir} already exist")

#check the contaminant list
contaminants = []
if path_blacklist !=None:
    print(path_blacklist)
    try :
     with open(path_blacklist, 'r') as f: #create a list of taxid from the blacklist 
         lines = f.readlines()
         for i, line in enumerate(lines) :
             line = line.strip()
             contaminants.append(line)
             message = f"taxid {line} has been added as a contaminant from the blacklist\n"
             print(message)
    except:
        print(f"Error : the blacklist file cannot be found \n") #return an error and stop the script if the blacklist is missing
        sys.exit()

## get main taxas
df_main_samples = pd.DataFrame(columns=['sampleID','taxID', 'taxName', 'reads_nb', 'taxo'])
#creation outliers file
taxon_dir= run_dir + "meta_analysis/"
try :
    os.makedirs(taxon_dir)
except :
    print(f"The run directory {taxon_dir} already exist")

#run taxonomic analysis
for SAGs in list_id: 
    empty_reports=[]
    sample_dir=f"{run_dir}{SAGs}"
    out_filter=f"{sample_dir}/filtered_reads"
    reports_dir = f"{out_filter}/kraken_reports/"
    krep_csv =f"{sample_dir}/{SAGs}_converted_kraken_report.csv"
    krep_sample = f"{reports_dir}{SAGs}_{suffix}_filter_reads_kr_reports.txt" #find kraken report for each sample in the id list
    df_taxas_main_sample = krep.convert_report(krep_sample) #convert report to a readable taxonomy
    main_taxa_sample = krep.extract_main_taxa(df_taxas_main_sample) #extract the main taxon for the sample
    main_taxa_sample['taxo'] = krep.convert_tax(main_taxa_sample['taxo']) #
    main_taxa_sample['sampleID']= SAGs
    df_main_samples.loc[len(df_main_samples)] = main_taxa_sample

df_main_samples.to_csv(f"{taxon_dir}Main_taxons_samples_{suffix}.csv", header=True)
#parse main taxon for anomalies
if main_abnormal_thresold !=1:
    occurrency_in_main = df_main_samples['taxID'].value_counts()
    if occurrency_in_main.iloc[0] > main_abnormal_thresold * len(list_id):     # if a taxa if found as main in more samples than the thresold...
        conta_name = occurrency_in_main.index[0]
        print(conta_name)
        message = f"WARNING : {conta_name} taxa is found as main in more than {main_abnormal_thresold} of the samples it main be a sign of contamination.\n"
        message += "This taxa was added to the contaminant list and removed from the reads, it can be avoided by adding the corresponding taxid in the blacklist and set thres_anomaly to 1 and relaunch the pipeline:\n"
        print(message)
        contaminants.append(conta_name) ###... this abnormal main taxon is added to the blacklist to be removed at the next stap of the pipeline...
        print(f"the new blacklist is {contaminants}, add the '--splash' option in the next step in order to run the read filtering again, it is needed to properly take those result into account.")


#new blacklist edition
blacklist_filter=f"{taxon_dir}Blacklist_filter_update.txt"
try:
	new_blacklist=open(blacklist_filter,'x')
except FileExistsError:
	print('the updated blacklist file already exist, overwritting it')
	new_blacklist=open(blacklist_filter,'w')
for tax_id in contaminants:
    new_blacklist.write(tax_id+"\n")
new_blacklist.close()
print(f"the new updated blacklist is written in {blacklist_filter} it should be used for downstream analysis.\n")

