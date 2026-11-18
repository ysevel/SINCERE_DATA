##!/usr/bin/python
import os,sys
import glob
import time
import matplotlib.pyplot as plt
import math
import shutil
from datetime import datetime
import numpy as np 
from collections import OrderedDict
import pandas as pd
from pathlib import Path
from Bio import Entrez
from textwrap import dedent as twdd
import argparse
import convert_krep as krep
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
                        help="Path to the SINCERE-DATA analysis folder used in the previous step(s).")
    parser.add_argument('-d', '--id', required=True,
                        help="ID for the sample, will be used to point to the correct files.")
    parser.add_argument('-o', '--output', required=True,
                        help="Output path where directories will be built.")
    parser.add_argument('-t', '--threads', type=int ,default=1,
                        help="Number of threads")
    parser.add_argument('-m', '--mem', action='store_true',
                        help="Allow kraken to run on memory (optionnal, otherwise will run kraken2 with the --memory-mapping option that avoid loading kraken database into the memory but require significantly higher computation time. If you want to load kraken2 database into the memory best to make sure you allocated enough memory the PFPplus kraken2 database require 300G and the nt database 750G to run properly for example).")
    parser.add_argument('-n', '--name', required=True,
                        help="Name of the run and main directory.")
    parser.add_argument('-f', '--fasta', action='store_true',
                        help="Add this argument if you want to generate a file of fasta sequences for the outliers found at this step.")
    parser.add_argument('-k', '--kraken_db', required=True,
                        help="Path to the kraken database folder.")
    parser.add_argument('-b', '--blacklist', required=False ,default=None,
                        help="Path to a blacklist file containing the taxon that should be systematically removed from the dataset (one taxid per line). Carefull it will eliminate all the children taxa as well!")
    parser.add_argument('--splash', action='store_true',
                        help='Use this option if you performed the step 2 of the pipeline and if it returned a splashome event: use this dataset as your own risk as it is likely too heavily contaminated to be of any value, run at your own risks.')
    parser.add_argument('-cs', '--confidence_score', required=False, default=0,
                        help="Kraken2 confidence score. See https://github.com/DerrickWood/kraken2/wiki/Manual#confidence-scoring for more details.")
    return parser.parse_args()


def remove_sample(id_sample, error, logger, missing_file, missing_id):
    """ report the sample id of a missing sample and write the id and the corresponding error message to the missing file."""
    with open(missing_file, 'a') as missf:
        missf.write(id_sample + '\t' + error + '\n')
        missing_id.append(id_sample)
    logger.error(f"Sample {id_sample} removed: {error}")


def check_sample(id_sample,missing_file,logger):
    try: 
        with open(missing_file,'r') as miss_f:
            miss_lines = miss_f.readlines()
            for line in miss_lines:
                line=line.strip()
                list_missing_id.append(line)
    except:
        logger.info("there is not missing file yet, the analysis will run on the complete dataset")
        return
    if id_sample in list_missing_id:
        logger.warning(f"{id_sample} has been found missing in previous steps and will been ignored, exiting.\n")
        sys.exit()
    else:
        logger.info(f"processing on preliminary assembly on {id_sample}\n")


def check_reads(path_input, id_sample, logger, missing_file, missing_id):
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
            remove_sample(id_sample, message, logger, missing_file, missing_id)
            return False, False
        elif len(R2_list) < 1:
            message = f"Error : R2 reads file not found for {id_sample}"
            remove_sample(id_sample, message, logger, missing_file, missing_id)
            return False, False
        else:
            R1_file=R1_list[0]
            R2_file=R2_list[0]
            logger.info(f"Both reads files found for sample {id_sample}")
            return R1_file, R2_file
    else:
        logger.error("Invalid input path specified in the index file. Please check if the path truly exists.")
        sys.exit(1)


def launch_kraken(R1_file, R2_file, id_sample, output, kraken_db, suffix, threads, mem, logger, sample_dir, confidence_score):
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
    command_kr=f"kraken2 --threads {threads} --db {kraken_db} --paired {R1_file} {R2_file} -output {kr_out} --use-names --report {kr_rep} {cmd_mem} --confidence {confidence_score}"
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


def genome_assembly(R1_file, R2_file, id_sample, suffix, output, threads, logger, sample_dir):
    """
    perform the assembly on the decontaminated reads, is supposed to be run with fixed parameters, also filter out contigs shorter than 500bp.
    """
    assembly_dir=f"{output}/assembly"
    try:
        os.makedirs(assembly_dir) #create a specific assembly dir, it is intended, at this step to be in the decontamination folder.
    except:
        print(f"Warning : sample directory {assembly_dir} already existing. Overwriting...")
    log_assembly=f"{sample_dir}{id_sample}_{suffix}_assembly.log"
    cmd_spades=f"spades.py --sc --careful -t {threads} -1 {R1_file} -2 {R2_file} -o {assembly_dir}/{id_sample}_{suffix}_spades" #run spade in a specific temporary folder that is supposed to be erased.
    logger.info(f"running spades with: {cmd_spades}")
    cmd_reformat=f"reformat.sh minlength=500 in={assembly_dir}/{id_sample}_{suffix}_spades/contigs.fasta  out={assembly_dir}/{id_sample}_{suffix}_prelim_assembly.fasta overwrite=true" #the trimmed assembly is saved directly in the assembly folder, outside of spades folder
    with open(log_assembly,'w') as log_a:
        subprocess.run(cmd_spades, stdout=log_a,stderr=subprocess.STDOUT, shell=True)
        subprocess.run(cmd_reformat, stdout=log_a,stderr=subprocess.STDOUT, shell=True) #the filtered assembly is writen directly in the assembly_dir
        if os.path.isdir(f"{assembly_dir}/{id_sample}_{suffix}_spades") :
            shutil.rmtree(f"{assembly_dir}/{id_sample}_{suffix}_spades") 
    logger.info(f"The preliminary assembly can be found in : {assembly_dir}/{id_sample}_{suffix}_prelim_assembly.fasta")
    return f"{assembly_dir}/{id_sample}_{suffix}_prelim_assembly.fasta"


def genome_alignement(R1_file, R2_file, assembly,id_sample, suffix, output, threads,logger, sample_dir):
    """This function perform reads alignement on the genome assembly, the final output is a bed file that can be used to determine the outliers regions"""
    align_dir=f"{output}/align_reads"
    try:
        os.makedirs(align_dir)
    except:
        f"Warning : sample directory {align_dir} already existing. Overwriting..."
    log_align=f"{sample_dir}{id_sample}_{suffix}_alignement.log"
    index=f"{output}/assembly/{id_sample}_{suffix}_prelim_assembly_index"
    cmd_index=f"bowtie2-build -f {assembly} {index}"
    sam_file=f"{align_dir}/{id_sample}_{suffix}.sam"
    cmd_bowtie=f"bowtie2 --threads {threads} -1 {R1_file} -2 {R2_file} -x {index} --sensitive -S {sam_file}"
    logger.info(f"Running bowtie2: {cmd_bowtie}")
    bam_file=f"{align_dir}/{id_sample}_{suffix}.bam"
    cmd_bam=f"samtools view -bS -o {bam_file} {sam_file}"
    bam_sorted=f"{align_dir}/{id_sample}_{suffix}_sorted.bam"
    logger.info(f"The sorted bam file for the sample in available at : {bam_sorted}")
    bai_sorted=f"{align_dir}/{id_sample}_{suffix}_sorted.bai"
    cmd_sort=f"samtools sort -o {bam_sorted} {bam_file}"
    cmd_bai=f"samtools index {bam_sorted} {bai_sorted}"
    bedgraph_file=f"{align_dir}/{id_sample}_{suffix}_bedgraph.tsv"
    cmd_bed=f"bedtools genomecov -ibam {bam_sorted} -bg > {bedgraph_file}"
    with open(log_align,'w') as log_al: #run all the commands on the same log file. ##note: should delete the alignement file the script will never use again 
        subprocess.run(cmd_index, stdout=log_al,stderr=subprocess.STDOUT, shell=True)
        subprocess.run(cmd_bowtie, stdout=log_al,stderr=subprocess.STDOUT, shell=True)
        subprocess.run(cmd_bam, stdout=log_al,stderr=subprocess.STDOUT, shell=True)
        os.remove(sam_file)
        subprocess.run(cmd_sort, stdout=log_al,stderr=subprocess.STDOUT, shell=True)
        os.remove(bam_file)
        subprocess.run(cmd_bai, stdout=log_al,stderr=subprocess.STDOUT, shell=True)
        subprocess.run(cmd_bed, stdout=log_al,stderr=subprocess.STDOUT, shell=True)
    logging.info(f"the bedgraph file is available in: {bedgraph_file}")
    return bedgraph_file #the bedgraph is the final file, the bam sorted file are used down the line as well


def get_avg(list):
    return(round(sum(list) / len(list),2))
    """ can be removed probably, math should perform this operation easily."""
#write_bed_file


def outlier_predictor(name_bedfile, name_bedgraph, threshold, align_d, outlier_d, sample_id, out_suffix, kraken_db, logger, sample_dir, nb_threads, mem, confidence_score):
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
    gap = 10 #maximal gap allowed between two outlier region to be merged
    min_length = 30 #minimal length of an outlier region to be taken ito consideration
    start = None
    nb_outlier=0       # number of regions with high cover depth
    passed=False     # by default, passed is false when the threshold hasn't been reached
    #variables a tester(optionnel?)
    outlier_dir=outlier_d
    align_dir=align_d
    id_sample=sample_id
    suffix=out_suffix
    with open(name_bedfile, 'w') as bedfile: 
        with open(name_bedgraph, 'r') as bedgraph: 
            nb_outlier=0       # number of regions with high cover depth
            lines = bedgraph.readlines()
            passed=False        # by default, passed is false when the threshold hasn't been reached
            waiting = False     # by default, there's no zone saved in memory waiting for a neighboring check
            outlier_region=""

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
                                outlier_region=f"{the_contig}:{start}-{end}"
                                process_outliers(align_dir, outlier_dir, sample_id, suffix, outlier_region, nb_outlier, kraken_db, logger, sample_dir, nb_threads, mem, confidence_score)
                                nb_outlier+=1                                 # ...the previous one is written...
                                bedfile.write(f"{last_contig}\t{last_start}\t{last_end}\tregion_{nb_outlier}\n")
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
                    passed = False                                          # ...region ended,...
                    if (int(end) - int(start) > min_length):    # ... if the outlier is long enough...
                        outlier_region=f"{the_contig}:{start}-{end}"
                        process_outliers(align_dir, outlier_dir, sample_id, suffix, outlier_region, nb_outlier, kraken_db, logger, sample_dir, nb_threads, mem, confidence_score)
                        nb_outlier+=1                                             # ...counted as a new one...
                        bedfile.write(f"{the_contig}\t{start}\t{end}\tregion_{nb_outlier}\n") # ...and added to the file
                if (int(line[3])<=threshold and passed):  # if cover ends exceeding the threshold...
                    passed = False              # ...region ended...
                    last_start=start                ##
                    last_end=end                     #
                    last_contig=the_contig           # 
                    waiting=True                     ## .... and saved and marked as waiting for next test 
                    #print(f"finished :  {line}\n")

            if waiting :    
                if (int(last_end) - int(last_start) > min_length):    # ... if the outlier is long enough...    
                    outlier_region=f"{the_contig}:{last_start}-{last_end}"
                    process_outliers(align_dir, outlier_dir, sample_id, suffix, outlier_region, nb_outlier, kraken_db, logger, sample_dir, nb_threads, mem, confidence_score)
                    nb_outlier += 1                                                                   ##
                    bedfile.write(f"{last_contig}\t{last_start}\t{last_end}\tregion_{nb_outlier}\n")   ## if last outlier is saved in memory, print it


#process_outliers
def process_outliers(align_dir, outlier_dir, id_sample, suffix, outlier_region, nb_outlier, kraken_db, logger, sample_dir, nb_thread, mem, confidence_score):
    try:
        os.makedirs(f"{outlier_dir}/outlier_fastq")
    except:
        print(f"Warning : sample directory {outlier_dir}/outlier_fastq already existing. Overwriting...")
    try:
        os.makedirs(f"{outlier_dir}/kraken_results")
    except:
        print(f"Warning : sample directory {outlier_dir}/kraken_results already existing. Overwriting...")
    try:
        os.makedirs(f"{outlier_dir}/kraken_reports")
    except:
        print(f"Warning : sample directory {output}/kraken_reports already existing. Overwriting...")
    kr_out=f"{outlier_dir}/kraken_results/{id_sample}_{suffix}_{nb_outlier}_kr_results.txt"
    kr_rep=f"{outlier_dir}/kraken_reports/{id_sample}_{suffix}_{nb_outlier}_kr_reports.txt"
    cmd_mem = "--memory-mapping" if not mem else ""
    log_outlier=f"{sample_dir}{id_sample}_{suffix}_outliers.log"
    cmd_samtools=f"samtools view -b -h {align_dir}/{id_sample}_{suffix}_sorted.bam {outlier_region} -o {outlier_dir}/outlier_fastq/{id_sample}_{suffix}_outlier_{nb_outlier}.bam"
    #print(cmd_samtools)
    cmd_sort=f"samtools sort -n -o {outlier_dir}/outlier_fastq/{id_sample}_{suffix}_outlier_{nb_outlier}_sorted.bam {outlier_dir}/outlier_fastq/{id_sample}_{suffix}_outlier_{nb_outlier}.bam"
    cmd_bamtofastq=f"bedtools bamtofastq -i {outlier_dir}/outlier_fastq/{id_sample}_{suffix}_outlier_{nb_outlier}_sorted.bam -fq {outlier_dir}/outlier_fastq/{id_sample}_{suffix}_outlier_{nb_outlier}_R1.fastq -fq2 {outlier_dir}/outlier_fastq/{id_sample}_{suffix}_outlier_{nb_outlier}_R2.fastq"
    cmd_kr_outlier=f"kraken2 --threads {nb_thread} --db {kraken_db} --paired {outlier_dir}/outlier_fastq/{id_sample}_{suffix}_outlier_{nb_outlier}_R1.fastq {outlier_dir}/outlier_fastq/{id_sample}_{suffix}_outlier_{nb_outlier}_R2.fastq -output {kr_out} --use-names --report {kr_rep} {cmd_mem} --confidence {confidence_score}"
    with open(log_outlier,'w') as log_outl:
        subprocess.run(cmd_samtools, stdout=log_outl,stderr=subprocess.STDOUT, shell=True)
        subprocess.run(cmd_sort, stdout=log_outl,stderr=subprocess.STDOUT, shell=True)
        os.remove(f"{outlier_dir}/outlier_fastq/{id_sample}_{suffix}_outlier_{nb_outlier}.bam") #remove outliers bam file as this can be huge and is not usefull later
        subprocess.run(cmd_bamtofastq, stdout=log_outl,stderr=subprocess.STDOUT, shell=True)
        #os.remove(f"{outlier_dir}/outlier_fastq/{id_sample}_{suffix}_outlier_{nb_outlier}_sorted.bam")
        subprocess.run(cmd_kr_outlier, stdout=log_outl,stderr=subprocess.STDOUT, shell=True)


#get visu cov (old)
def get_visu_cov(cover, threshold, fig_name, id_file):
    """ Draw and save a jpg file named [[fig_name]] that shows values of coverage accross the genome and the value of the threshold
        defined depending on z-score (2) ; returns True if successful """
    try :
        fig, ax = plt.subplots()
        ax.plot(cover, label='cover')
        ax.axhline(y=threshold, color='r', label=f'threshold : {threshold}')
        ax.set_ylabel('coverage')
        ax.set_xlabel(f'position on assembly')
        ax.set_title(f"Visualization of outliers regions on {id_file}.")
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


#fonction threshold
def threshold_calculator (id_sample, suffix, align_dir, decontam_dir,logger, sample_dir, kraken_db, get_Fasta, nb_threads, mem, output_missing_file, confidence_score):
    """ This fonction calculate the z-score threshold on coverage in order to find outliers regions in the preliminary assembly. The mean coverage for the whole assembled is computed and Z-score calculaed, outlier regions are determined as region with a coverage equal or supperior to 2 std over the mean coverage (mean + 2*std)"""
    outliers_dir=f"{decontam_dir}/outliers"
    df_summary = pd.DataFrame(columns=['file', 'mean', 'max', 'standard error', 'length', 'threshold'])
    with open(f"{sample_dir}{id_sample}_{suffix}_outliers.log", "a") as log_t :
        try:
            os.makedirs(outliers_dir)
        except:
            print(f"Warning : sample directory {outliers_dir} already existing. Overwriting...")
        try:
            with open(f"{align_dir}/{id_sample}_{suffix}_bedgraph.tsv", 'r') as bedg:
                contigs, pos, cover = [], [], []
                line = bedg.readline()
                line_nb = 1
                while line != "" :
                    line = line.split("\t")
                    contigs.append(line[0])
                    try:
                        begin=int(line[1])
                        end=int(line[2])
                    except:
                        message=f"position values not recognized in line {line_nb} :\n{line}\naborting\n"
                        log_t.write(message)
                        remove_sample(id_sample, message, logger, output_missing_file, current_missing_id_list)
                    while begin < end :
                        try : 
                            cover.append(int(line[3]))
                            begin+=1
                        except :
                            error=(f"coverage column not found, or invalid data, at line {line}; aborting\n")
                            remove_sample(id_sample, error, logger, output_missing_file, current_missing_id_list)
                            continue
                    line = bedg.readline()
                    line_nb += 1
        except :
            message=f"ERROR : bedgraph file could not be opened for {id_sample}, check the log file if bedgraph generation failed\n"
            log_t.write(message)  
            remove_sample(id_sample, message, logger, output_missing_file, current_missing_id_list)

        try : 
            maxi=max(cover)     
        except :                # if file is empty....
            row=[id_sample, None, None, None, None, None]
            df_summary.loc[len(df_summary)]=row
            message=f"{id_sample}_{suffix}_bedgraph.tsv is empty"
            log_t.write(message)
            logger.error(message)
            remove_sample(id_sample, message, logger, output_missing_file, current_missing_id_list)# ...remove the sample from thecollectipon          
        
    ### get values to define the threshold
        mean = get_avg(cover)
        std = np.std(cover)
        threshold = int(mean + 2*std)

    ## add values to summary file
        row=[id_sample, mean, maxi, round(std,2), len(cover), threshold]
        df_summary.loc[len(df_summary)]=row

        ## draw coverage graph and save it
        fig_name=f"{decontam_dir}/{id_sample}_{suffix}_visu_cov.png"
        if (get_visu_cov(cover, threshold, fig_name, id_sample) == True):
            logger.info(f"Visualization available at: {fig_name}")
        else : 
            log_t.write(f"{fig_name} generation failed")

        ### writing bedfile
        outlier_bedfile = f"{outliers_dir}/{id_sample}_{suffix}_outliers.bed"
        sample_bedgraph = f"{align_dir}/{id_sample}_{suffix}_bedgraph.tsv"
        outlier_predictor(outlier_bedfile, sample_bedgraph, threshold, align_dir, outliers_dir, id_sample, suffix,kraken_db, logger, sample_dir, nb_threads, mem, confidence_score)

        ## writing and displaying results
        if get_Fasta == True: ### generates fasta files [option] 
            fasta_file=f"{decontam_dir}/assembly/{id_sample}_{suffix}_peaks.fasta"
            os.system(f"bedtools getfasta -fi {decontam_dir}/assembly/{id_sample}_{suffix}_prelim_assembly.fasta -bed {outlier_bedfile} -name > {fasta_file}")
            message = f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} : {outlier_bedfile} and {fasta_file} generated\n"
        else :
            message = f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} : {outlier_bedfile} generated\n"
        log_t.write(message)
        logger.info(f"{outlier_bedfile} generated")
        return threshold



def process_sample(id_sample, path_input, suffix, path_blacklist, kraken_db, mem, run_dir,nb_threads,get_Fasta,splash,confidence_score):
    ###check missing sample
    missing_file=f"{run_dir}missing_id.txt"
    sample_dir=f"{run_dir}{id_sample}/"
    taxon_dir= f"{run_dir}meta_analysis/"
    filter_dir=f"{sample_dir}filtered_reads"
    decontam_dir=f"{sample_dir}decontamination"
    #implement the log for step3:
    log_file = f"{sample_dir}outliers_prediction_{id_sample}_{suffix}.log"
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
    
    logger.info(f"Processing sample: {id_sample}")
    missing_file = f"{run_dir}missing_id.txt"
    check_sample(id_sample, missing_file, logger)
    #folder verification
    if os.path.exists(sample_dir):
        logger.info(f"Running the outlier detection for {id_sample} in {sample_dir} directory")
        logger.info(f"Successfully found {id_sample}, proceding to read analysis.")
    else:
        if id_sample in missing_id:
            logger.warning("the sample was removed from previous step of the pipeline, ignoring the sample")
            sys.exits()
        else:
            logger.error("ERROR, the sample directory do not exist, please check the path")
            sys.exit()

    #reads to use
    reads_dir=f"{sample_dir}filtered_reads/reads"
    reads_R1, reads_R2 = check_reads(reads_dir, id_sample, logger, missing_file, missing_id)
    if not reads_R1 or not reads_R2:
        return
    #variables for kraken2
    filter_suffix=f"{suffix}_filter_reads"
    kr_filter_out=f"{filter_dir}/kraken_results/{id_sample}_{filter_suffix}_kr_results.txt"
    kr_filter_rep=f"{filter_dir}/kraken_reports/{id_sample}_{filter_suffix}_kr_reports.txt"
    #second filtering stap if given the splash option, otherwise will ignore this step and copy the filtered reads from the first step into the decontam_dir
    if splash==True: #if the splash option is given a second round of reads filtering is run on the dataset, otherwise 
        contaminants = []
        if path_blacklist !=None:
            print(path_blacklist)
            try :
                with open(path_blacklist, 'r') as bl: 
                    lines = bl.readlines()
                    for i, line in enumerate(lines) :
                        line = line.strip()
                        contaminants.append(line)
                        logger.info(f"taxid {line} has been added as a contaminant from the blacklist\n")
            except: #return an error if no blacklist have been provided.
                logger.info(f"Error : the blacklist file cannot be found \n")
                return
        #filtering out the new blacklist's taxa and run kraken on the filtered reads
        contaminants_str=' '.join(f"{str(e)}" for e in contaminants)
        logger.info(f"taxon selected for filtering: {contaminants_str} from the blacklist \nNow running reads filtering\n")
        filter_reads_R1, filter_reads_R2 = blacklist_filtering(reads_R1, reads_R2, id_sample, decontam_dir, kr_filter_out, filter_suffix, contaminants_str, kr_filter_rep,logger, sample_dir)
        #run kraken on the newly filtered reads including the splashome taxon.
        kr_filter_out,kr_filter_reports = launch_kraken(filter_reads_R1, filter_reads_R2, id_sample,decontam_dir,kraken_db,filter_suffix,nb_thread, mem, logger, sample_dir,confidence_score)
    else: #if no --splash option is provided, just copy the already existing filtered reads into the decontam_dir and copy kraken reports and results in order to have them in the correct folder.
        filter_reads_R1=reads_R1
        filter_reads_R2=reads_R2
        try: #creation of two directories for reports and results from kraken2
            os.makedirs(f"{decontam_dir}/kraken_results")
        except:
            print(f"Warning : sample directory {decontam_dir}/kraken_results already existing. Overwriting...")
        try:
            os.makedirs(f"{decontam_dir}/kraken_reports")
        except:
            print(f"Warning : sample directory {decontam_dir}/kraken_reports already existing. Overwriting...")
        filter_reads_dir=f"{decontam_dir}/reads"
        try:
            os.makedirs(filter_reads_dir)
        except:
            print(f"Warning : sample directory {filter_reads_dir} already existing. Overwriting...")
        shutil.copy(kr_filter_out,f"{decontam_dir}/kraken_results")
        shutil.copy(kr_filter_rep,f"{decontam_dir}/kraken_reports")
        shutil.copy(filter_reads_R1, filter_reads_dir)
        shutil.copy(filter_reads_R2, filter_reads_dir)
        logger.info(f"Copying {filter_reads_R1} and {filter_reads_R2} from {reads_dir} to {filter_reads_dir}")

    #genome assembly
    assembly_dir=f"{decontam_dir}/assembly"
    align_dir=f"{decontam_dir}/align_reads"
    logger.info(f"running preliminary assembly for {id_sample}")
    prelim_assembly=genome_assembly(filter_reads_R1, filter_reads_R2, id_sample, suffix, decontam_dir, nb_thread,logger, sample_dir)
    #command for debug, allow to bybass time consuming steps once they have run for the first time
    #filter_reads_R1=f"{output_dir}/{id_sample}_{filter_suffix}_R1.fastq"
    #filter_reads_R2=f"{output_dir}/{id_sample}_{filter_suffix}_R2.fastq"
    #prelim_assembly=f"{assembly_dir}/{id_sample}_{suffix}_prelim_assembly.fasta"

    #reads alignements
    logger.info(f"running reads alignement for {id_sample}")
    genome_bed=genome_alignement(filter_reads_R1, filter_reads_R2,prelim_assembly ,id_sample, suffix, decontam_dir, nb_thread,logger, sample_dir)
    logger.info(f"calculating threshold for {id_sample}")
    threshold=threshold_calculator (id_sample, suffix, align_dir, decontam_dir,logger, sample_dir, kraken_db, get_Fasta, nb_threads, mem, missing_file, confidence_score)
    logger.info(f"the threshold for {id_sample} is {threshold}")
    logger.info(f"script end for sample {id_sample}") #print the end time of the script, in order to have easy access to the total time the script ran.

#main and variable initialisation
if __name__ == "__main__":
    arguments = handle_program_options()
    #parse the argument, most of them would correspond to previous step of this pipeline.
    path_input = arguments.input
    id_sample = arguments.id
    output = arguments.output           
    suffix = arguments.name               ### only one output path and suffix as meta-analysis is performed should be the same as in the previous step of the workflow
    path_blacklist = arguments.blacklist
    kraken_db = arguments.kraken_db
    mem = arguments.mem
    nb_thread = arguments.threads
    get_Fasta = arguments.fasta
    splash=arguments.splash
    confidence_score=arguments.confidence_score

    run_dir=f"{output}/"
    taxon_dir= f"{run_dir}meta_analysis/"
    if os.path.exists(run_dir):
        print("Now running preliminary assembly and outlier prediction, this might take a while\n")

    else:
        print(f"the directory {run_dir} does not exist, it is intended to be created with the script Sincere-data_reads_filtering.py as it is the first step of this pipeline, please check your path or run Sincere-data_reads_filtering.py and, optionnally, splashome_filter.py on the samples first before running this script again.")
        sys.exit()
    if os.path.exists(taxon_dir):
        print(f"sucessfully found the taxon directory: {taxon_dir}\n")
    else:
        print(f"WARNING: the directory {taxon_dir} does not exist, it is likely that the script splashome_filter.py have not been run on the samples prior to running this Script. As the splahome_filter.py is an important part of the pipeline you should consider running it before running this script again.")
    missing_id = []
    process_sample(id_sample, path_input,suffix, path_blacklist, kraken_db,mem, run_dir,nb_thread,get_Fasta,splash,confidence_score)
