##!/usr/bin/python
import os,sys
import glob
import time
import math
from datetime import datetime
#import numpy as np 
from collections import OrderedDict
import pandas as pd
from pathlib import Path
from Bio import Entrez
from textwrap import dedent as twdd
import argparse
import convert_krep as krep

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
    parser.add_argument('-o', '--output', required=True,
                        help="Output path where directories will be built")
    parser.add_argument('-n', '--name', required=True,
                        help="suffix added to each output file")
    parser.add_argument('-ts', '--thres_samples', type=float, default=0.1,
                        help="Allowed proportion for a taxa to be found in samples \
                            and considered as contaminant when less than the float value, must be a value between 0 and 1")
    parser.add_argument('-to', '--thres_outliers', type=float, default=0.2,
                        help="Minimal proportion for a taxa found in outliers to be \
                            considered as contaminant when upper than the float value, must be a value between 0 and 1")
    parser.add_argument('-b', '--blacklist', required=False ,default=None,
                        help="path to a blacklist file containing the taxon that should be systematically removed from the dataset (one taxid per line). Carefull it will eliminate all the children taxa as well!")
    parser.add_argument('-e', '--email', required=True, help='an email adress to be provided to Entrez as it is required to use E-utilities.')
    return parser.parse_args()

def remove_sample(id_sample,error):
    with open(missing_file, 'a') as missf:
        missf.write(id_sample + '\t' + error)
        missing_id.append(id_sample)
        print(id_sample)

def check_missing_samples(path_index_file,path_missing_file):
    try :
        with open(path_index_file, 'r') as id_f: 
            lines = id_f.readlines()
            for line in lines:
                line = line.strip()
                list_id.append(line)
    except :
        print(f"Error : index file not found at {path_index_file}")
        #sys.exit(1)
    try: 
        with open(path_missing_file,'r') as miss_f:
            miss_lines = miss_f.readlines()
            for line in miss_lines:
                line=line.strip()
                list_missing_id.append(line)
    except:
        print("there is not missing file yet, the analysis will run on the complete dataset")
        return
    for sample_id in list_id:
        if sample_id in list_missing_id:
            list_id.remove(sample_id)
            print(f"{sample_id} has been found missing in previous steps and will been removed from the samples.")

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

def taxName2taxID(list_names,email):
    """ as outliers' main taxas may be linked to samples' main taxas through going up in the taxonomy, taxID needs to be recovered using 
        Entrez module. 
        Input : 
            list of taxonomic names  
        Output : 
            dictionnary {taxName : taxID}, list of taxNames which failed their identification
    """
    Entrez.email = email
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

### script
#variable initiation
arguments = handle_program_options()
list_id=[]
list_missing_id = []
path_index_file = arguments.input
output = arguments.output           
suffix = arguments.name               ### only one output path and suffix as meta-analysis is performed
path_blacklist = arguments.blacklist
thres_outliers = arguments.thres_outliers       
thres_samples = arguments.thres_samples
email = arguments.email


#implement path and verification step
run_dir=f"{output}/"
taxon_dir= f"{run_dir}meta_analysis/"

#redirect all prints to a general log file
general_log = f"{run_dir}meta_analysis_{suffix}.log"
print(f"All prints will now be saved in {general_log} file.")
if os.path.exists(general_log):    # generates a new file
    os.remove(general_log)
#sys.stdout = open(general_log, 'a')

# #same for errors
err_log = f"{run_dir}meta_analysis_{suffix}_{suffix}.err"
print(f"{err_log} file generated for possible errors.")
if os.path.exists(err_log):  
    os.remove(err_log)
#sys.stderr = open(err_log, 'a')

message = f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} : script begin for run {run_dir}\n"
print(message)

#folder verification
if os.path.exists(run_dir):
    print("Now run directory found...\n")
else:
    print(f"the directory {run_dir} does not exist, it is intended to be created with the script Sincere-data_reads_filtering.py as it is the first step of this pipeline, please check your path or run Sincere-data_reads_filtering.py and splashome_filter.py on the samples first before running this script again.")
    sys.exit()

if os.path.exists(taxon_dir):
    print(f"sucessfully found the taxon directory: {taxon_dir}\n")
else:
    print(f"WARNING: the directory {taxon_dir} does not exist, it is likely that the script splashome_filter.py have not been run on the samples prior to running this Script. As the splahome_filter.py is an important part of the pipeline you should consider running it before running this script again.")
    mkdir(taxon_dir)
missing_file=f"{run_dir}missing_id.txt"
check_missing_samples(path_index_file, missing_file)
nb_files=len(list_id)

### checking the lecture of index file 
if list_id == []:
    print(f"Error : index file {path_index_file} is empty")
    sys.exit()

for sample in list_id:
    sample_dir=f"{run_dir}{sample}"
    decontam_dir=f"{sample_dir}/decontamination"
    outliers_dir=f"{decontam_dir}/outliers"
if os.path.exists(outliers_dir):
    print(f"sucessfully found the taxon directory for {sample}.\n")
else:
    print(f"the directory {outliers_dir} does not exist, for {sample} it is intended to be created with the script Sincere-data_outlier_predictor.py as it is the previous step of this pipeline, please check your path or run Sincere-data_reads_filtering.py, splashome_filter.py and Sincere-data_outlier_predictor.py on the samples first before running this script again.")





#check the contaminant list
contaminants = []
if path_blacklist !=None:
    print(path_blacklist)
    try :
     with open(path_blacklist, 'r') as f: 
         lines = f.readlines()
         for i, line in enumerate(lines) :
             line = line.strip()
             contaminants.append(line)
             message = f"taxid {line} has been added as a contaminant from the blacklist\n"
             print(message)
    except:
        print(f"Error : the blacklist file cannot be found \n")
        sys.exit()

###meta-analyse
#initialization
df_main_samples = pd.DataFrame(columns=['sampleID','taxID', 'taxName', 'reads_nb', 'taxo'])
df_main_outliers_unique_strains = pd.DataFrame(columns=['sampleID','taxID', 'taxName', 'reads_nb', 'taxo'])    # dataframe havins one example of each strain found as main in each outlier
list_global_outliers=[]


###### work in progress ci dessous/////

print(list_id)
#extract and store main sample taxon and outliers main taxons
for SAGs in list_id:
    list_outliers=[]
    empty_reports=[]
    outliers_id=SAGs
    sample_dir=f"{run_dir}{SAGs}"
    decontam_dir=f"{sample_dir}/decontamination"
    outliers_dir=f"{decontam_dir}/outliers"
    reports_main = f"{decontam_dir}/kraken_reports/"
    reports_outlier=f"{outliers_dir}/kraken_reports/"
    krep_sample = f"{reports_main}{SAGs}_{suffix}_filter_reads_kr_reports.txt"
    df_taxas_main_sample = krep.convert_report(krep_sample)
    main_taxa_sample = krep.extract_main_taxa(df_taxas_main_sample)
    main_taxa_sample['taxo'] = krep.convert_tax(main_taxa_sample['taxo'])
    main_taxa_sample['sampleID']= SAGs
    df_main_samples.loc[len(df_main_samples)] = main_taxa_sample
    df_main_outliers = pd.DataFrame(columns=['sampleID','taxID', 'taxName', 'reads_nb', 'taxo'])
    valid_reports = len(glob.glob(f"{reports_outlier}{SAGs}*kr_reports.txt"))
    print(f"there is {valid_reports} kraken reports for {SAGs}.")
    for outlier in Path(reports_outlier).glob('*kr_reports.txt') :
        krep_name_full = os.path.basename(outlier)             ## get full name of outlier's report 
        krep_name = krep_name_full.replace(suffix + "_","")
        krep_name = krep_name.replace("_report.txt","")
        df_taxas_outliers = krep.convert_report(outlier)
        if len(df_taxas_outliers) == 0:                             ##
            valid_reports -= 1                                       #
            empty_reports.append(krep_name_full)                     #
            continue                                                 ## if the dataframe is empty, its name is saved and next report is processed   
        else : 
            main_taxa_outlier = krep.extract_main_taxa(df_taxas_outliers)               ## extraction of main taxa from outliers
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
    message = f"\nAmong {total} outliers, {nb_others} were found having a main taxa different than the one found in the main sample reads\n"
    message = message + f"Taxa \t\t occurrency_in_outliers\n"
    if len(empty_reports) != 0:
        nb_outliers=len(glob.glob(f"{reports_outlier}{SAGs}*kr_reports.txt"))
        message = f"{valid_reports} reports out of {nb_outliers} outlier were significant and will be exploited further.\n"
        message += "Please check :\n" + "\n".join(rep for rep in empty_reports)
        message += f"\nFrom {reports_outlier} if any doubt"
        print(f"For {SAGs}, {message}")
    list_global_outliers.append(taxas_outliers)
    df_main_outliers.to_csv(f"{outliers_dir}/{SAGs}_outliers.csv", index=True)
#print(df_main_outliers_unique_strains)
df_main_samples.to_csv(f"{taxon_dir}Main_taxons_samples_{suffix}.csv", header=True)
df_main_outliers_unique_strains.to_csv(f"{taxon_dir}Main_taxons_outliers_{suffix}.csv", header=True)
print(empty_reports)

taxas_found = count_taxa_outliers(list_global_outliers)
taxas_id, err = taxName2taxID(taxas_found.keys(),email)
if err :
    print(f"WARNING : conversion(s) of taxName failed for the following taxas :")
    print(" ".join(e for e in err))
    print("taxas with a * have several taxas IDs referenced that will be used for decontamination. Others don't have any ID proposed \n")
    # print(dico_convert)
#creation of metanalyse report
df_meta = meta_analysis(taxas_found, df_main_samples, taxas_id, suffix)
meta_file = f"{taxon_dir}results_meta_analysis_{suffix}.csv"
df_meta.to_csv(meta_file, header=True)
print(f"Results of meta-analysis available in file : \n{meta_file} \n")

for index, row in df_meta.iterrows():
    if row['ratio_outliers']>thres_outliers and row['ratio_samples']<thres_samples:
        contaminants.append(row['taxID'])
        message = f"{row['taxName']} of ID {row['taxID']} has been added as contaminant. \n"
        print(message)

if len(contaminants) == 0:
    message = f"ERROR : no contaminant having an occurrency >{thres_outliers} in outliers and <{thres_samples} in samples in this run.\n"
    message = message + f"Please consider change settings of values 'thres_outliers' and 'thres_samples' while launching the code, the script will run on the blacklist alone \n"
    print(message)


#write the final blacklist
#new blacklist edition
contam_file=f"{taxon_dir}contaminants_{suffix}.txt"
try:
	contams_final=open(contam_file,'x')
except FileExistsError:
	print('the updated blacklist file already exist, overwritting it')
	contams_final=open(contam_file,'w')
for tax_id in contaminants:
    contams_final.write(tax_id+"\n")
contams_final.close()
print(f"the new updated blacklist is written in {contam_file} it should be used for downstream analysis.\n")