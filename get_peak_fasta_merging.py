# -- coding: utf-8 --

import csv
import matplotlib.pyplot as plt
import numpy as np
import math
import os

def get_avg(list):
    """return the average of a list's values"""
    return(round(sum(list) / len(list),2))

def get_std(list, mean):
    """return the standard deviations of a list's values"""
    variance = sum((x - mean) ** 2 for x in list) / len(list)
    return(round(math.sqrt(variance),2))

samples=("AG-345-I14","AG-321-F02","AG-347-K23","AG-390-I16","AG-463-P23","AG-429-E20","AG-363-A05","AG-424-A14","AG-410-B18","AG-410-D23", "AG-337-G06", "AG-347-K17")

suffix="berudes_rd0"
assembly_path="/groups/ecogeno/DEV/2023_RobertNoe/bin/align_run6_berudes/rd_0/assemblies/assemblies_500bp_filtered/"
assembly_suffix="_berudes_rd0_ecrem.fasta"
input_dir="/groups/ecogeno/DEV/2023_RobertNoe/bin/align_run6_berudes/rd_0/"

for filename in samples:
    ### getting coverage
    #with open(f"{path_input}{id_file}/bedtools/bedgraphs/{id_file}_{suffix}_bedgraph.tsv", 'r') as f: 
    try : 
        with open(f"{input_dir}bedtools/bedgraphs/{filename}_{suffix}_bedgraph.tsv", 'r') as f: 
            contigs=[]
            pos = []
            cover=[]
            line = f.readline()
            line_nb = 1
            while line != "" :
                line = line.split("\t")
                contigs.append(line[0])
                try:
                    begin=int(line[1])
                    end=int(line[2])
                except:
                    print(f"position values not recognized in line {line_nb} :\n{line}\naborting")
                    break
                while begin < end :
                    try : 
                        cover.append(int(line[3]))
                        begin+=1
                    except :
                        print (f"coverage column not found, or invalid data, at line {line_nb} :\n{line}\naborting")
                        break
                    
                line = f.readline()
                line_nb += 1
    except :
        print(f"{filename} failed")
        continue

    try : 
        maxi=max(cover)     # if file is empty....
    except : 
        continue            # ...abort file filling
    
    mean = get_avg(cover)
    std = get_std(cover,mean)

    threshold = int(mean+2*std)
    start = None
    regions = []

    # ### Building the output directories
    # if not os.path.exists("bedtools/peaks/"):
    #     os.makedirs("bedtools/peaks/")
    # if not os.path.exists("bedtools/peaks_fasta/"):
    #     os.makedirs("bedtools/peaks_fasta/")

    # ### Peak detection using bedgraph data
    # with open(f"bedtools/peaks/{filename}_{suffix}.bed", 'w') as output: 
    print(f"{filename}_{suffix}.bed")
    with open(f"{filename}_{suffix}.bed", 'w') as bedfile:
        with open(f"{input_dir}bedtools/bedgraphs/{filename}_{suffix}_bedgraph.tsv", 'r') as bedgraph: 
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
                        if (int(line[1]) - int(last_end) < 10 and line[0]==last_contig) :    # and if one finished less than 10 nuc before...
                            passed = True                           # ...threshold is count as exceeded,...
                            start=last_start                        # ...previous zone's beginning is used, i.e they're merged...
                            last_start,last_end,last_contig=0,0,""  # ...last values are reset (is it enough ?)...
                            waiting = False                         # ...and there's no longer a zone waiting for check    
                        else :                          # else (i.e it's really a new exceeding region)
                            nb_fasta+=1                         # the previous one is written...
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
                    nb_fasta+=1                                             # ...counted as a new one...
                    bedfile.write(f"{the_contig}\t{start}\t{end}\tregion_{nb_fasta}\n") # ...and added to the file
                    #print(f"changed of contig: {line}")

                if (int(line[3])<=threshold and passed):  # if cover ends exceeding the threshold...
                    passed = False              # ...region is ended...
                    last_start=start                ##
                    last_end=end                     #
                    last_contig=the_contig           # 
                    waiting=True                     ## .... and saved and marked as waiting for next test 
                    #print(f"finished :  {line}\n")

            if waiting :
                nb_fasta += 1
                bedfile.write(f"{last_contig}\t{last_start}\t{last_end}\tregion_{nb_fasta}\n") # ...and added to the file

    # generate fasta files
    os.system(f"bedtools getfasta -fi {assembly_path}/{filename}{assembly_suffix} -bed bedtools/peaks/{filename}_{suffix}.bed -name > bedtools/peaks_fasta/{filename}_{suffix}_peaks.fasta")
    print(f"bedtools/peaks/{filename}_{suffix}.bed and bedtools/peaks_fasta/{filename}_{suffix}_peaks.fasta generated")