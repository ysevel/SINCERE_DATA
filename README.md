# SINCERE DATA : SINgle-CEll REads Decontamination through Automatic Taxonomic Assignation
This is the version 1.6.0 of the Sincere-DATA pipeline. 

The pipeline is divided in six different steps that perform:

- reads filtering of undesired taxons or known contaminants (Sincere-data_reads_filtering.py)
- screening the whole dataset for splashome type contamination (splashome_filter.py)
- research of outlier regions and outlier analysis (Sincere-data_outlier_predictor.py)
- meta-analysis of the outliers to compute a list of potential contaminants (Sincere-data_meta_analysis.py)
- decontamination of the reads and contigs based on this contaminant list (Sincere-data_decontamination.py)
- final assembly filtering through completeness and contamination metrics (Sincere-data_assembly_filtering.py)

## Contents
- [Overview](#overview)
- [Installation](#installation)
- [Dependencies](#dependencies)
- [Usage](#usage)
- [Parameters](#parameters)
- [Main output](#main-output)
- [Citation](#citation)

## Overview
SINCERE-DATA is a workflow divided in six steps (Four mandatory and two optional):

First, a step of reads filtering to remove known contaminant or undesirable taxons that are likely present in the dataset. The reads filtering step is performed using Kraken2, a widely used taxonomic classifier and contamination assessment tool. Kraken2 uses exact k-mer matching to predict the taxonomic classification of a short or long DNA sequence based on a taxonomic database. Kraken2 provide several well maintained and regularly updated databases depending on the type of data-set to decontaminate. The filtering is performed by using a blacklist of undesirable taxids and to screen them based on each reads classification through krakentools utility.

An optional step is conducted to screen the data-set for splashome contamination. As Single-Cell is especially vulnerable to cross contamination between samples. To do so, the main taxon is identified for each sample and prevalence analysis is conducted. An alert is returned when a taxon encompassing more than the user-defined proportion of the SAGs. As the proportion can vary with the data-set, as a starting point, we propose screening for a taxon represented in more than 70% of the SAGs in a diverse collection. This option can be turned off when working with pure strains or with an uniform data-set. If a cross contamination event is detected the contaminant taxon will be added to an updated blacklist and removed from the filtered reads.

In the second main step, the workflow will look for abnormally high coverage region in a preliminary assembly and identify high coverage regions based on z-score analysis (deviation from the mean). A Z-score threshold of 2 standard deviation above the mean have been selected as abnormally high for this purpose. Those high coverage regions are extracted and the reads are aligned on those regions. Taxonomic analysis of those regions is performed to extract the main taxon.

The third step is a meta-analysis conducted for each sample in order to identify the outliers regions that present both a abnormally high read coverage and a main taxon that differ from the whole sample main taxon. Those outliers regions are compiled and the frequency of each outlier taxon is reported in order to identify foreign DNA contaminants. The selection of the contaminants to consider for cleaning is performed thank to two user defined thresholds. (1) the 'outlier prevalence' threshold, which determines the proportion of SAGs (i.e. as a percentage of the data-set.) containing a given outlier region taxonomic identification and (2) the 'sample taxonomy' threshold, which defines the proportion of SAGs (i.e. the percentage of the data-set.) with a main taxon identical to a detected potential contaminant in the outlier regions. Outlier regions both above the 'outlier prevalence threshold' and below the 'sample taxonomy threshold' were validated as contaminants and their associated DNA sequences (i.e. the sequencing reads of the contaminant and associated inferior taxonomic ranks) were removed from the whole data-set using KrakenTools.

After that, the read decontamination is performed in the fourth step of the workflow based on a contaminant list from both the blacklist defined by the user and the outliers taxa identified during the meta-analysis. Those contaminants are filtered out of the reads and a final assembly is conducted. The contaminant identified are also eliminated from the assembly by filtering out the contigs that are classified to a contaminant taxon. The workflow also perform quality assessment on this cleaned assembly using QUAST. During this step the contigs removed by Sincere-data can be kept for further analysis. An utility have been implemented to retrieve plasmid associated contigs from this contigs: Sincere-data use Platon to perform an analysis on the "junk" contig removed during the decontamiation step in order to reintegrate them into the cleaned assembly. Annotation of the cleaned assembly can also be performed using bakta. 

The final optional step of the workflow can filtered out low quality assembly. First, it will launch a global CheckM evaluation of all assemblies. Based on this evaluation, assemblies will be filtered out using two defined threshold for completeness and contamination. To be kept the assembly have to show a completeness higher (or equal) than the completeness threshold and a contamination lower (or equal) than the contamination threshold.

![workflow_sincere-data-final](https://github.com/ysevel/SINCERE_DATA/tree/main/Image/Workflow.png)

## Installation
SINCERE-DATA doesn't need any installation. Scripts can be run as described in Usage section. However, dependencies need to be installed either in a conda environment or in your path (See list below). For KrakenTools, you only need to download the extract_kraken_reads.py script and put it in `SINCERE_DATA/scripts`. 

## Dependencies
### Mandatory
- bbmap V.39.09 + : <https://github.com/bbushnell/BBTools>
- bedtools V.2.31.1 + : <https://github.com/arq5x/bedtools2>
- biopython V.1.84 : <https://github.com/biopython/biopython>
- bowtie2 V.2.5.4 + : <https://github.com/BenLangmead/bowtie2>
- Kraken2 V.2.1.3 + : <https://github.com/DerrickWood/kraken2>
- KrakenTools : <https://github.com/jenniferlu717/KrakenTools>
- matplotlib-base V.3.9.2 + : <https://github.com/matplotlib/matplotlib>
- numpy V.1.26 + : <https://github.com/numpy/numpy>
- pandas V.2.2.2 + : <https://github.com/pandas-dev/pandas>
- QUAST V.5.2.0 + : <https://github.com/ablab/quast>
- SAMtools V.1.21 + : <https://github.com/samtools/samtools>
- SPAdes V.4.0.0 + : <https://github.com/ablab/spades>

### Step or option dependant
#### If using the `-a --annotation` in Step 4
- bakta V.1.8.1 + : <https://github.com/oschwengers/bakta> 

#### If using the `-j --keep_junk_contigs` option in Step 4
- Platon V.1.7 + : <https://github.com/oschwengers/platon>

#### For the optional assembly filtering step
- CheckM V.1.2.5 + : <https://github.com/Ecogenomics/CheckM>


## Usage
### Classical usage (step by step)
#### Step 1: Read filtering
```bash
python Sincere-data_reads_filtering.py -i [reads_dir] -d [sample_id] -o [output_dir] -n [suffix] -k [path_to_kraken_db] (-b [blacklist_file] -t [nb_thread] -cs [kraken2 confidence score, default = 0] -m)
```

#### Optional step: Splashome filter
```bash
python splashome_filter.py -i [sample_id_list.txt] -o [output_dir] -n [suffix] (-b [blacklist_file.txt] -ta [threshold_anomaly, default=1])
```

#### Step 2: Outlier prediction
```bash
python Sincere-data_outlier_predictor.py -i [input_dir] -d [sample_id] -o [output_dir] -n [suffix] -k [path_to_kraken_db] (-b [blacklist_file] -t [nb_thread] -cs [kraken2 confidence score, default = 0] -m -f --splash)
```

#### Step 3: Meta-analysis
```bash
python Sincere-data_meta_analysis.py -i [sample_id_list.txt] -o [output_directory] -n [suffix] -e [email] (-b [blacklist_file.txt] -ts [sample ratio, default = 0.1] -to [outlier ratio, default = 0.2])
```

#### Step 4: Decontamination
```bash
python Sincere-data_decontamination.py -i [input_dir] -d [sample_id] -o [output_dir] -n [suffix] -k [path_to_kraken_db] (-c [contaminants_list.txt] -t [nb_thread] -j -ba [path_to_bakta_db] -p [path_to_platon_db] -cs [kraken2 confidence score, default = 0] -m)
```

#### Optional step: Assembly filtering
```bash
python Sincere-data_assembly_filtering.py -i [sample_id_list.txt] -o [output_dir] -n [suffix] (-t [nb_thread] -co [completeness threshold, default = 50] -ct [contamination threshold, default = 10])
```
### Slurm cluster specific usage
We have developped scripts that permit a more easy usage for slurm cluster users in one command line. First, you need to go into the Slurm_scripts directory :
```bash
cd Slurm_scripts
```
Then, you can modify the Configuration file present in the directory. This file contains all the parameters described below and some others option like the activation of optional step or the memory you want to use on the cluster. After the configuration file is modified, you just need to run the `Wrapper_sincere.sh` scripts with sbatch : 
```bash
sbatch Wrapper.sh
```
This script was tested on the GenOuest cluster maybe it will need modifications to work on other cluster. For example, the line 69 in `Launch_sincere.sh` file used to activate the conda environment is dependent on how your cluster works. 

## Parameters
### Required parameters
`-i --input` Path to the input directory or files dependent of the step. Path to the trimmed reads directory for step 1. Path to the sample id txt file for splashome filter step, step 3 and assembly filtering step. Path to the SINCERE-DATA analysis directory for step 2 and 4.

`-d --id_file` ID of the sample to analyze. (Only for steps 1, 2 and 4)

`-o --output` Path to the output directory (directory of the SINCERE-DATA analysis).

`-n --name` Suffix for the output files and run directory.

`-k --kraken_db` Path to a kraken database. (Only for steps 1, 2 and 4)

### Optional parameters
`-cs --confidence_score` (default=0) Kraken2 confidence score. See <https://github.com/DerrickWood/kraken2/wiki/Manual#confidence-scoring> for more details. (Only for steps 1, 2 and 4)

`-b --blacklist` (default=none) Path to a blacklist file containing the taxon that should be systematically removed from the dataset (one taxid per line). Careful it will eliminate all the children taxa as well.

`-m --mem` Allow kraken to run on memory (optionnal, otherwise will run kraken2 with the `--memory-mapping` option that avoid loading kraken database into the memory but require significantly higher computation time. If you want to load kraken2 database into the memory best to make sure you allocated enough memory the PFPplus kraken2 database require 300G and the nt database 750G to run properly for example) (see kraken2 manual for more informations). (Only for steps 1, 2 and 4)

`-t --threads` (default=1) Number of threads to use. (Only for steps 1, 2, 4 and assembly filtering)

### Step-specific parameters
#### Splashome finder
`-ta --thres_anomaly` (default=1) Must be a value between 0 and 1, indicate the proportion of the sample allowed to present the same main taxa. For no verification of this parameter keep the default threshold at 1.

#### Step 2
`--splash` Use this option if you performed the splashome analysis functionnality of the pipeline and if it returned a splashome event: use this dataset as your own risk as it is likely too heavily contaminated to be of any value.

`-f --fasta` Add this argument if you want to generate a file of fasta sequences for the outliers found at this step.

#### Step 3
`-e --email` (required) An email adress to be provided to Entrez as it is required to use E-utilities.

`-ts --thres_sample` (default=0.1) Allowed proportion for a taxa to be found in samples and considered as contaminant when less than the float value.

`-to --thres_outlier` (default=0.2) Minimal proportion for a taxa found in the collection to be considered as contaminant when upper than the float value.

#### Step 4
`-c --contaminants` (default=None) Path to a text file containing the tax ID of taxa that should be systematically removed from the dataset (one tax ID per line). Careful it will eliminate all the children taxa as well! This file would ideally be generated from the previous steps of this workflow, but it can also be a 'blacklist' file with taxa the user want to exclude.

`-a --annotation` If this option is selected the annotation of the assembly will be performed with bakta.

`-ba --bakta_db` The path to bakta database for genome annotation (required if `-a --annotation` option is enabled).

`-j --keep_junk_contigs` If this option is selected the contigs discarded by the assembly cleaning will be stored and analyzed for possible mobile elements.

`-p --platon_db ` The path to platon database for plasmid identification (required if `-j --keep_junk_contigs` option is enabled).

#### Assembly filtering
`-co --completeness_threshold` (default=50) Must be an integer between 0 and 100, indicate the minimum completeness percentage (evaluated by CheckM) to keep the assembly.

`-ct --contamination_threshold` (default=10) Must be an integer between 0 and 100, indicate the maximum contamination percentage (evaluated by CheckM) to keep the assembly.

## Main output
After running SINCERE-DATA, you will find a directory for each sample submitted which contains multiple results for all analysis. You will also find two others directory : 
- `meta_analysis` which contains the splashome filter step (if used) and the step 3 results. 
- `cleaned_assembly` which contains final assemblies and annotations (if performed) and filtered assemblies if this step was done. 

The main output for each step are described below. If you are only interested in the final result of SINCERE-DATA they will be in the `cleaned assembly` directory. 

### Step 1 : Reads filtering
- `<sample_id>/trimmed_reads/kraken_results` : Directory containing the Kraken2 result file for the input trimmed reads in step 1.
- `<sample_id>/trimmed_reads/kraken_reports` : Directory containing the Kraken2 report file for the input trimmed reads in step 1.
- `<sample_id>/filtered_reads/reads` : Directory containing the filtered reads after step 1.
- `<sample_id>/filtered_reads/kraken_results` : Directory containing the Kraken2 result file for the filtered reads after step 1.
- `<sample_id>/filtered_reads/kraken_reports` : Directory containing the Kraken2 report file for the filtered reads after step 1.

### Optional splashome filter step
- `meta_analysis/Blacklist_filter_update.txt` : Updated blacklist after splashome filter analysis, should be used for next analysis.

### Step 2 : Outlier prediction
- `<sample_id>/decontamination/assembly/<sample_id>_<run_name>_prelim_assembly.fasta` : Preliminary assembly by SPAdes that will be use as a reference for the outlier detection. 
- `<sample_id>/decontamination/align_reads/<sample_id>_<run_name>_sorted.bam` : Bowtie2 reads mapping on the preliminaly assembly to identify outlier regions.
-  `<sample_id>/decontamination/outliers/<sample_id>_<run_name>_outliers.bed` : Identified outlier regions.
- `<sample_id>/decontamination/<sample_id>_<run_name>_visu_cov.png` : Visual representation of the reads alignment on the preliminary assembly. 

### Step 3 : Meta-analysis
- `meta_analysis/Main_taxons_samples_<run_name>.csv` : Table containing the main taxon identified for each samples. 
- `meta_analysis/Main_taxons_outliers_<run_name>.csv` : Table containing the main taxon identified for each significant outliers i.e. with a different identification from the main sample taxonomy. 
- `meta_analysis/results_meta_analysis_<run_name>.csv` : Table containing the statistics for each potential identified contaminants.
- `meta_analysis/contaminants_<run_name>.txt` : Final list of identified contaminants (passing both threshold) that can be used for the step 4 of the workflow.

### Step 4 : Decontamination
- `<sample_id>/final_assembly/assembly` : Directory containing the SPAdes assembly of the decontaminated reads. 
- `<sample_id>/final_assembly/kraken_raw_assembly_results` : Directory containing the Kraken2 result file for contigs from the decontaminated reads assembly.
- `<sample_id>/final_assembly/kraken_raw_assembly_reports` : Directory containing the Kraken2 report file for contigs from the decontaminated reads assembly.
- `<sample_id>/cleaned_assembly/<sample_id>_<run_name>_cleaned_assembly.fasta` : Final cleaned assembly after removal of contaminants from the previous assembly. 
- `<sample_id>/cleaned_assembly/quast` : Directory containing the QUAST analysis of the cleaned assembly. 
- `<sample_id>/cleaned_assembly/kraken_raw_assembly_results` : Directory containing the Kraken2 result file for contigs from the cleaned assembly.
- `<sample_id>/cleaned_assembly/kraken_raw_assembly_reports` : Directory containing the Kraken2 report file for contigs from the cleaned assembly.
- `<sample_id>/cleaned_assembly/annotations` : Directory containing the bakta annotation of the cleaned assembly (if performed).
- `<sample_id>/junk_contigs/<sample_id>_<run_name>_junk_contigs.fasta` : Junk contigs which will be analyzed (if the option is enabled).
- `<sample_id>/junk_contigs/quast` : Directory containing the QUAST analysis of the junk_contigs (if the option is enabled). 
- `<sample_id>/junk_contigs/mobile_elements` : Directory containing the Platon analysis of the junk_contigs to found mobile elements (if the option is enabled). 
- `cleaned_assembly/assembly_fasta` : Directory containing all cleaned assemblies.
- `cleaned_assembly/genes_faa` : Directory containing all cleaned assemblies aminoacid sequences of genes annotated by bakta.
- `cleaned_assembly/genbanks` : Directory containing all cleaned assemblies annotation by bakta in gbff format.

### Optional assembly filtering step
- `cleaned_assembly/checkm/report_checkm_<run_name>.tsv` : CheckM report of all cleaned assemblies. 
- `cleaned_assembly/checkm/report_filtered_checkm_<run_name>.tsv` : CheckM report of filtered assemblies. 
- `cleaned_assembly/checkm/Filtered_samples_list_<run_name>.tsv` : List of samples that passed the completeness and contamination threshold. 
- `cleaned_assembly/filtered_assembly` : Directory containing all filtered assemblies.
- `cleaned_assembly/filtered_genes` : Directory containing all filtered assemblies aminoacid sequences of genes annotated by bakta.
- `cleaned_assembly/filtered_genbanks` : Directory containing all filtered assemblies annotation by bakta in gbff format.


## Citation
If you use SINCERE-DATA in your research, please cite this paper :
> Mauger, S., Sevellec, Y., Carret, L., Robert, N., Monard, C., Thion, C., Jarassier, L., & Vandenkoornhuyse, P. (2025). Wet-lab preparation and automated dataset decontamination procedure for bacterial single-cell genomics. Methods in Ecology and Evolution, 16, 2819-2835. https://doi.org/10.1111/2041-210X.70170

Please, don't forget to also cite all the other softwares used in SINCERE-DATA workflow :
- bbmap : <https://github.com/bbushnell/BBTools>
- bedtools : <https://doi.org/10.1093/bioinformatics/btq033>
- bowtie2 : <https://doi.org/10.1038/nmeth.1923>
- Kraken2 : <https://doi.org/10.1186/s13059-019-1891-0>
- KrakenTools : <https://doi.org/10.1038/s41596-022-00738-y>
- QUAST : <https://github.com/ablab/quast>
- SAMtools : <https://doi.org/10.1093/gigascience/giab008>
- SPAdes : <https://doi.org/10.1002/cpbi.102>
- bakta : <https://doi.org/10.1099/mgen.0.000685> 
- Platon : <https://doi.org/10.1099/mgen.0.000398>
- CheckM : <https://doi.org/10.1101/gr.186072.114>