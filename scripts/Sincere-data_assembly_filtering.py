##!/usr/bin/python
import os,sys
import argparse
import subprocess
import shutil
from textwrap import dedent as twdd


def handle_program_options():
    descr = """\
    Complete pipeline of single-cell decontamination using \
    taxonomic assignation. 
    Program arguments
    -----------------"""

    parser = argparse.ArgumentParser(description=twdd(descr),
                           formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input', required=True,
                        help="Path to index file of samples to process.")
    parser.add_argument('-o', '--output', required=True,
                        help="Output path where directories will be built.")
    parser.add_argument('-n', '--name', required=True,
                        help="Suffix added to each output file.")
    parser.add_argument('-t', '--threads', type=int ,default=1,
                        help="Number of threads")
    parser.add_argument('-co', '--completeness_threshold', required=False, type=int, default=50,
                        help="Must be a integer between 0 and 100, indicate the minimum completeness percentage to keep the assembly evaluated by CheckM.")
    parser.add_argument('-ct', '--contamination_threshold', required=False, type=int, default=10,
                        help="Must be a integer between 0 and 100, indicate the maximum contamination percentage to keep the assembly evaluated by CheckM.")
    return parser.parse_args()


def cleaned_assembly_filtering(output, suffix, nb_threads, completeness_threshold, contamination_threshold):
    keep_list = []
    cleaned_dir = f"{output}/cleaned_assembly"
    checkm_dir = f"{cleaned_dir}/checkm"
    checkm_outfile = f"{checkm_dir}/report_checkm_{suffix}.tsv"
    assembly_dir = f"{cleaned_dir}/assembly_fasta"

    filtered_assembly_dir = f"{cleaned_dir}/filtered_assembly"
    try:
        os.makedirs(filtered_assembly_dir)
    except:
        f"Warning : assembly directory {filtered_assembly_dir} already existing. Overwriting..."
    if os.path.isdir(f"{cleaned_dir}/genes_faa"):
        filtered_genes_dir = f"{cleaned_dir}/filtered_genes"
        try:
            os.makedirs(filtered_genes_dir)
        except:
            f"Warning : assembly directory {filtered_genes_dir} already existing. Overwriting..." 

        filtered_genbanks_dir = f"{cleaned_dir}/filtered_genbanks"
        try:
            os.makedirs(filtered_genbanks_dir)
        except:
            f"Warning : assembly directory {filtered_genbanks_dir} already existing. Overwriting..." 

    filtered_checkm_outfile = f"{checkm_dir}/report_filtered_checkm_{suffix}.tsv"
    cmd_checkm = f"checkm lineage_wf -t {nb_threads} -x fasta --tab_table -f {checkm_outfile} {assembly_dir} {checkm_dir}"
    print("running CheckM on cleaned assembly: {cmd_checkm}")
    subprocess.run(cmd_checkm, stderr=subprocess.STDOUT, shell=True)
    print("Evaluating the assembly using its completness and contamination calculated by CheckM...")
    with open(checkm_outfile, "r") as cm_file:
        lines = cm_file.readlines()
        with open(filtered_checkm_outfile, "w") as fcm_file:
            for i in range(0, len(lines)):
                line = lines[i]
                if i == 0:
                    fcm_file.write(line) #Write the header of the checkM result file
                else:
                    tabline = line.split("\t")
                    id = tabline[0].replace(f"_{suffix}_cleaned_assembly", "")
                    comp = float(tabline[11])
                    conta = float(tabline[12])
                    if comp >= completeness_threshold:
                        print(f"Assembly {id} passed the completeness threshold")
                        if conta <= contamination_threshold:
                            print(f"Assembly {id} passed the contamination threshold")
                            fcm_file.write(line)
                            sample_cleaned_assembly_dir = f"{output}/{id}/cleaned_assembly"
                            sample_assembly_file = f"{sample_cleaned_assembly_dir}/{id}_{suffix}_cleaned_assembly.fasta"
                            shutil.copy(sample_assembly_file, filtered_assembly_dir)
                            if os.path.isdir(f"{cleaned_dir}/genes_faa"):
                                sample_genes_file = f"{sample_cleaned_assembly_dir}/annotations/{id}_{suffix}.faa"
                                try:
                                    shutil.copy(sample_genes_file, filtered_genes_dir)
                                except:
                                    f"Warning : no amino acid gene sequence were found for {id}."
                                sample_genbank_file = f"{sample_cleaned_assembly_dir}/annotations/{id}_{suffix}.gbff"
                                try:
                                    shutil.copy(sample_genbank_file, filtered_genbanks_dir)
                                except:
                                    f"Warning : no annotation file were found for {id}."
                            keep_list.append(id)
                        else:
                            print(f"Assembly {id} didn't pass the contamination threshold, filtering it out")
                    else:
                        print(f"Assembly {id} didn't pass the completeness threshold, filtering it out")
    filtered_sample_file = f"{cleaned_dir}/Filtered_samples_list_{suffix}.txt"
    with open(filtered_sample_file, "w") as fs_file:
        for sample in keep_list:
            fs_file.write(f"{sample}\n")




#variable initiation
arguments = handle_program_options()
path_index_file = arguments.input
output = arguments.output           
suffix = arguments.name
nb_threads = arguments.threads
completeness_threshold = arguments.completeness_threshold
contamination_threshold = arguments.contamination_threshold

cleaned_assembly_filtering(output, suffix, nb_threads, completeness_threshold, contamination_threshold)
