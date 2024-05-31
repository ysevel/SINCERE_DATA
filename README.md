# SINCERE DATA : SINgle-CEll REads Decontamination through Automatic Taxonomic Assignation

This is the version 0.1 of the sincere_data pipeline actually provided for the genouest cluster users. 

The script is optimized to work with conda environment with the required programs.

It will run sarray bash script but an updated version will be provided soon that rely solely on python scripts.


To install the tools, create a conda environment with all the requiered tools and drop all the script from the bash files and scripts folder into the directory you wish to run the program from.

## Requirement:

fastqc V. 0.11+

spades V. 3.15+

kraken2 V. 2.1+

bowtie2 V. 2.4+

samtools V. 1.9+

bedtools V. 2.30+

krakentools : https://github.com/jenniferlu717/KrakenTools

quast V. 5.0+

### Python library requirement:

panda

numpy

Biopython

matplotlib

## Overview: 

This tools is designed to work on already trimmed reads and will make a preliminary assembly based on spades in singlecell mode. 
An reads alignment is performed on the assembled genome using bowtie2 and samtools and the coverage depth is calculated using bedtools.
it will then evaluate the reads alignment depth and detect outlier region based on the Z-score (number of square deviation from the mean coverage) More specifically, outlier being defined as the assembly region of more than 30bp presenting a z-score >2. 

kraken2 is used to detect the main taxon from the assembly and the outlier regions to identify potential contaminants based on divergences in the taxonomy and to detect contamination that overlap with the other single cell assembled genomes (SAGs) from a collection. The potential contaminant are 

## Usage:

python workflow_decontaxo_SAGs.py -i [input datafile] -o [output directory] -n [output_suffix] -ts [sample ratio default = 0.1]  -to [|outlier ratio, default =0.2]-b [blacklist file] -k [kraken db path] -c [conda environment to use with the bash files (beta version only)]

 ### Arguments :

-i --input (required): input tvs file containing the sample names and their path

-o --output (required): path to the output folder

-n --name (required): suffix for the output file and run directory

-f --fasta (default =false): Give False as an argument to stop generation of fasta sequences of outliers found
-k –kraken_db (required): Path to the kraken db 

-c -conda_env (required): name of the conda environment the pipeline is installed in (beta version, soon to be removed)

-ts –thres_sample (default= 0.1): Allowed proportion for a taxa to be found in samples and considered as contaminant when less than the float value

-to –thres_outlier (default= 0.2): Minimal proportion for a taxa found in the collection to be considered as contaminant when upper than the float value

-b –blacklist (default=none) : path to a blacklist file containing the taxon that should be systematically removed from the dataset (one taxid per line). Careful it will eliminate all the children taxa as well!

-ta thres_anomaly (default=1) must be a value between 0 and 1, indicate the proportion of the sample allowed to present the same main taxa. For no verification of this parameter keep the default thresold at 1.

