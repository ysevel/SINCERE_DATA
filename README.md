# SINCERE DATA : SINgle-CEll REads Decontamination through Automatic Taxonomic Assignation

This is the version 1.1.0 of the sincere_data pipeline actually provided for the genouest cluster users. 

The script is optimized to work with conda environment with the required programs.

The pipeline consiste in five different step that perform

- reads filtering of undesired taxons or know contaminants (Sincere-data_reads_filtering.py)
- screening the whole dataset for sphlashome type of contamination (splashome_filter.py)
- research of outlier regions and outlier analysis (Sincere-data_outlier_predictor.py)
- méta-analysis of the outliers to compute a list of potential contaminants (Sincere-data_meta_analysis.py)
- decontamination of the reads and contigs based on this contaminant list (Sincere-data_decontamination.py)


To install the tools, create a conda environment with all the requiered tools and drop all the script from the bash files and scripts folder into the directory you wish to run the program from. Then, download Krakentools from https://github.com/jenniferlu717/KrakenTools and drop all the files in the Sincere-DATA installation directory.

## Requirement:
bbmap V.39.09 +

bedtools V.2.31.1 +

biopython V.1.84

bowtie2 V.2.5.4 +

kraken2 V.2.1.3 +

krakentools : https://github.com/jenniferlu717/KrakenTools

matplotlib-base V.3.9.2 +

numpy V. 1.26 +

pandas V.2.2.2 +

platon V.1.7 +

quast V.5.2.0 +

samtools V.1.21 +

spades V.4.0.0 +

## Overview: 

SINCERE-DATA is a workflow divided in three steps:

First, a step of reads filtering to remove known contaminant or undesirable taxons that are likely present in the dataset. The reads filtering step is performed using Kraken2, a widely used taxonomic classifier and contamination assessment tool. Kraken2 use exact k-mer matching to predict the taxonomic classification of a short or long DNA sequence based on a taxonomic database. Kraken2 provide several well maintained and regularly updated databases depending on the type of data-set to decontaminate. The filtering is performed by using a blacklist of undesirable taxids and to screen them based on each reads classification through krakentools utility.

An optional step is conducted to screen the data-set for widespread contamination. As singlecell is especially vulnerable to cross contamination between samples. To do so, the main taxon is identified for each sample and prevalence analysis is conducted. An alert is returned when a taxon encompassing more than the user-defined proportion of the SAGs. As the proportion can vary with the data-set, as a starting point, we propose screening for a taxon represented in more than 70% of the SAGs in a diverse collection. This option can be turned off when working with pure strains or with a uniform data-set. If a cross contamination event is detected the contaminant taxon will be added to an updated blacklist and removed from the filtered reads.

In the second step the workflow will look for abnormally high coverage region in a preliminary assembly and identify high coverage regions based on z-score analysis (deviation from the mean). A Z-score threshold of 2 standard deviation above the mean have been selected as abnormally high for this purpose.  Those high coverage regions are extracted and the reads are aligned on those regions. Taxonomic analysis is performed of those regions to extract the main taxon.

A meta-analysis is conducted for each sample in order to identify the outliers regions that present both a abnormally high read coverage and a main taxon that differ from the whole sample main taxon. Those outliers regions are compiled and the frequency of each outlier taxon is reported in order to identify foreign DNA contaminants. The selection of the contaminants to consider for cleaning is performed thank to two user defined thresholds. (1) the ‘outlier prevalence’ threshold, which determines the proportion of SAGs (i.e. as a percentage of the data-set.) containing a given outlier region and (2) the ‘sample taxonomy’ threshold, which defines the proportion of SAGs (i.e. the percentage of the data-set.) with a main taxon identical to a detected potential contaminant in the outlier regions. Outlier regions both above the ‘outlier prevalence threshold’ and below the ‘sample taxonomy threshold’ were validated as contaminants and their associated DNA sequences (i.e. the sequencing reads of the contaminant and associated inferior taxonomic ranks) were removed from the whole data-set. using KrakenTools.

Finally, the read decontamination is performed in the third step of the workflow based on a contaminants list from both the blacklist defined by the user and the outliers taxa identified during the meta-analysis. Those contaminants are filtered out of the reads and a final assembly is conducted. The contaminant identified are also eliminated from the assembly by filtering out the contigs that are classified to a contaminant taxon. The workflow also perform quality assessment on this cleaned assembly using quast. During this step the contigs removed by Sincere-data can be kept for further analysis. an utility have been implemented to retrieve plasmid associated contigs from the contigs: Sincere-data use Platon to perform an analysis on the "junk" co,tig removed during the decontamiation stp in order to reintegrate them into the cleaned assembly.

## Usage:

python workflow_decontaxo_SAGs.py -i [input datafile] -o [output directory] -n [output_suffix] -ts [sample ratio, default = 0.1]  -to [|outlier ratio, default =0.2]-b [blacklist file] -k [kraken db path] -c [conda environment to use with the bash files (beta version only)]

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

