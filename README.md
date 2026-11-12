# SINCERE DATA : SINgle-CEll REads Decontamination through Automatic Taxonomic Assignation
This is the version 1.6.0 of the Sincere-DATA pipeline. 

The script is optimized to work with conda environment with the required programs.

The pipeline is divided in six different steps that perform:

- reads filtering of undesired taxons or known contaminants (Sincere-data_reads_filtering.py)
- screening the whole dataset for splashome type of contamination (splashome_filter.py)
- research of outlier regions and outlier analysis (Sincere-data_outlier_predictor.py)
- meta-analysis of the outliers to compute a list of potential contaminants (Sincere-data_meta_analysis.py)
- decontamination of the reads and contigs based on this contaminant list (Sincere-data_decontamination.py)
- final assembly filtering through completness and contamination metrics (Sincere-data_assembly_filtering.py)


To install the tools, create a conda environment with all the required tools and drop all the script from the bash files and scripts folder into the directory you wish to run the program from. Then, download Krakentools from https://github.com/jenniferlu717/KrakenTools and drop all the files in the Sincere-DATA installation directory.

## Requirements
### Mandatory
bbmap V.39.09 + : https://github.com/bbushnell/BBTools

bedtools V.2.31.1 + : https://github.com/arq5x/bedtools2

biopython V.1.84 : https://github.com/biopython/biopython

bowtie2 V.2.5.4 + : https://github.com/BenLangmead/bowtie2

Kraken2 V.2.1.3 + : https://github.com/DerrickWood/kraken2

KrakenTools : https://github.com/jenniferlu717/KrakenTools

matplotlib-base V.3.9.2 + : https://github.com/matplotlib/matplotlib

numpy V.1.26 + : https://github.com/numpy/numpy

pandas V.2.2.2 + : https://github.com/pandas-dev/pandas

QUAST V.5.2.0 + : https://github.com/ablab/quast

SAMtools V.1.21 + : https://github.com/samtools/samtools

SPAdes V.4.0.0 + : https://github.com/ablab/spades

### Step or option dependant
#### To perform the automatic annotation in Step 4
bakta V.1.8.1 + : https://github.com/oschwengers/bakta 

#### If using the '--keep_junk_contigs' option in Step 4
Platon V.1.7 + : https://github.com/oschwengers/platon

#### Assembly filtering step
CheckM V.1.2.5 + : https://github.com/Ecogenomics/CheckM 

## Overview
SINCERE-DATA is a workflow divided in six steps (Four mandatory and two optional):

First, a step of reads filtering to remove known contaminant or undesirable taxons that are likely present in the dataset. The reads filtering step is performed using Kraken2, a widely used taxonomic classifier and contamination assessment tool. Kraken2 uses exact k-mer matching to predict the taxonomic classification of a short or long DNA sequence based on a taxonomic database. Kraken2 provide several well maintained and regularly updated databases depending on the type of data-set to decontaminate. The filtering is performed by using a blacklist of undesirable taxids and to screen them based on each reads classification through krakentools utility.

An optional step is conducted to screen the data-set for splashome contamination. As Single-Cell is especially vulnerable to cross contamination between samples. To do so, the main taxon is identified for each sample and prevalence analysis is conducted. An alert is returned when a taxon encompassing more than the user-defined proportion of the SAGs. As the proportion can vary with the data-set, as a starting point, we propose screening for a taxon represented in more than 70% of the SAGs in a diverse collection. This option can be turned off when working with pure strains or with a uniform data-set. If a cross contamination event is detected the contaminant taxon will be added to an updated blacklist and removed from the filtered reads.

In the second main step the workflow will look for abnormally high coverage region in a preliminary assembly and identify high coverage regions based on z-score analysis (deviation from the mean). A Z-score threshold of 2 standard deviation above the mean have been selected as abnormally high for this purpose. Those high coverage regions are extracted and the reads are aligned on those regions. Taxonomic analysis of those regions is performed to extract the main taxon.

The third step is a meta-analysis conducted for each sample in order to identify the outliers regions that present both a abnormally high read coverage and a main taxon that differ from the whole sample main taxon. Those outliers regions are compiled and the frequency of each outlier taxon is reported in order to identify foreign DNA contaminants. The selection of the contaminants to consider for cleaning is performed thank to two user defined thresholds. (1) the ‘outlier prevalence’ threshold, which determines the proportion of SAGs (i.e. as a percentage of the data-set.) containing a given outlier region and (2) the ‘sample taxonomy’ threshold, which defines the proportion of SAGs (i.e. the percentage of the data-set.) with a main taxon identical to a detected potential contaminant in the outlier regions. Outlier regions both above the ‘outlier prevalence threshold’ and below the ‘sample taxonomy threshold’ were validated as contaminants and their associated DNA sequences (i.e. the sequencing reads of the contaminant and associated inferior taxonomic ranks) were removed from the whole data-set. using KrakenTools.

After that, the read decontamination is performed in the fourth step of the workflow based on a contaminant list from both the blacklist defined by the user and the outliers taxa identified during the meta-analysis. Those contaminants are filtered out of the reads and a final assembly is conducted. The contaminant identified are also eliminated from the assembly by filtering out the contigs that are classified to a contaminant taxon. The workflow also perform quality assessment on this cleaned assembly using QUAST. During this step the contigs removed by Sincere-data can be kept for further analysis. An utility have been implemented to retrieve plasmid associated contigs from the contigs: Sincere-data use Platon to perform an analysis on the "junk" contig removed during the decontamiation stp in order to reintegrate them into the cleaned assembly. Annotation of the cleaned assembly can also be performed using bakta. 

The final optional step of the workflow can filtered out low quality assembly. First, it will launch a global CheckM evaluation of all assemblies. Based on this evaluation, assemblies will be filtered out using two chosen threshold for completeness and contamination. To be kept the assembly have to show a completeness higher (or equal) than the completeness threshold and a contamination lower (or equal) than the contamination threshold.

![workflow_sincere-data-final](https://github.com/ysevel/SINCERE_DATA/tree/main/Image/Workflow.png)


## Usage
### Step 1: Read filtering
python Sincere-data_reads_filtering.py -i [reads_dir] -d [sample_id] -o [output_dir] -n [suffix] -k [path_to_kraken_db] (-b [blacklist_file] -t [nb_thread] -cs [kraken2 confidence score, default = 0] -m)

### Optional step: Splashome filter
python splashome_filter.py -i [sample_id_list.txt] -o [output_dir] -n [suffix] (-b [blacklist_file.txt] -ta [threshold_anomaly, default=1])

### Step 2: Outlier prediction
python Sincere-data_outlier_predictor.py -i [input_dir] -d [sample_id] -o [output_dir] -n [suffix] -k [path_to_kraken_db] (-b [blacklist_file] -t [nb_thread] -cs [kraken2 confidence score, default = 0] -m -f --splash)

### Step 3: Meta-analysis
python Sincere-data_meta_analysis.py -i [sample_id_list.txt] -o [output_directory] -n [suffix] -e [email] (-b [blacklist_file.txt] -ts [sample ratio, default = 0.1] -to [outlier ratio, default = 0.2])

### Step 4: Decontamination
python Sincere-data_decontamination.py -i [input_dir] -d [sample_id] -o [output_dir] -n [suffix] -k [path_to_kraken_db] (-c [contaminants_list.txt] -t [nb_thread] -j -ba [path_to_bakta_db] -p [path_to_platon_db] -cs [kraken2 confidence score, default = 0] -m)

### Optional step: Assembly filtering
python Sincere-data_assembly_filtering.py -i [sample_id_list.txt] -o [output_dir] -n [suffix] (-t [nb_thread] -co [completeness threshold, default = 50] -ct [contamination threshold, default = 10])

## Parameters
### Required parameters
-i --input (required) : Path to the input folders or files. Path to the trimmed reads folder for step 1. Path to the sample id txt file for splashome filter, step3 and assembly filtering. Path to the SINCERE-DATA analysis folder for step 2 and 4.

-d --id_file (required) : ID of the sample to analyze. (Only for steps 1, 2 and 4)

-o --output (required) : Path to the output folder (folder of the SINCERE-DATA analysis).

-n --name (required) : Suffix for the output file and run directory.

-k --kraken_db (required) : Path to a kraken database. (Only for steps 1, 2 and 4)

### Optional parameters
-cs --confidence_score (default=0) : Kraken2 confidence score. See https://github.com/DerrickWood/kraken2/wiki/Manual#confidence-scoring for more details. (Only for steps 1, 2 and 4)

-b --blacklist (default=none) : Path to a blacklist file containing the taxon that should be systematically removed from the dataset (one taxid per line). Careful it will eliminate all the children taxa as well.

-m --mem : Allow kraken to run on memory (optionnal, otherwise will run kraken2 with the --memory-mapping option that avoid loading kraken database into the memory but require significantly higher computation time. If you want to load kraken2 database into the memory best to make sure you allocated enough memory the PFPplus kraken2 database require 300G and the nt database 750G to run properly for example) (see kraken2 manual for more informations). (Only for steps 1, 2 and 4)

-t --threads (default=1) : Number of threads to use. (Only for steps 1, 2, 4 and assembly filtering)

### Step-specific parameters
#### Splashome finder
-ta --thres_anomaly (default=1) : Must be a value between 0 and 1, indicate the proportion of the sample allowed to present the same main taxa. For no verification of this parameter keep the default threshold at 1.

#### Step 2
--splash : Use this option if you performed the splashome analysis functionnality of the pipeline and if it returned a splashome event: use this dataset as your own risk as it is likely too heavily contaminated to be of any value.

-f --fasta : Add this argument if you want to generate a file of fasta sequences for the outliers found at this step.

#### Step 3
-e --email (required) : An email adress to be provided to Entrez as it is required to use E-utilities.

-ts --thres_sample (default=0.1): Allowed proportion for a taxa to be found in samples and considered as contaminant when less than the float value.

-to --thres_outlier (default=0.2): Minimal proportion for a taxa found in the collection to be considered as contaminant when upper than the float value.

#### Step 4
-c --contaminants (default=None) : Path to a text file containing the tax ID of taxa that should be systematically removed from the dataset (one tax ID per line). Careful it will eliminate all the children taxa as well! This file would ideally be generated from the previous steps of this workflow, but it can also be a 'blacklist' file with taxa the user want to exclude.

-a --annotation : If this option is selected the annotation of the assembly will be performed with bakta.

-ba --bakta_db : The path to bakta database for genome annotation. (required and only needed if '--annotation' option is enabled)

-j --keep_junk_contigs : If this option is selected the contigs discarded by the assembly cleaning will be stored and analyzed for possible mobile elements.

-p --platon_db : The path to platon database for plasmid identification. (required and only needed if '--keep_junk_contigs' option is enabled)

#### Assembly filtering
-co --completeness_threshold (default=50) : Must be an integer between 0 and 100, indicate the minimum completeness percentage (evaluated by CheckM) to keep the assembly.

-ct --contamination_threshold (default=10) : Must be an integer between 0 and 100, indicate the maximum contamination percentage (evaluated by CheckM) to keep the assembly.