from collections import OrderedDict
import re
import sys 
import os.path as osp
import csv
from pathlib import Path
import pandas as pd
import os

ranks = ["D", "P", "C", "O", "F", "G", "S", "SS"]
field_names = ["pct_reads", "clade_reads", "taxon_reads", "rank", "ncbi_tax", "sci_name"]

maxi_rank="D"
mini_rank="S"   

def convert_report(krep_fp):
    """
    Input :
        krep_fp : path to kraken report
    Inspired from kraken-biom codes. Converts a kraken report into complete taxonomy visualisation, and returns it as 
    a 4-column dataframe : 
        "taxID"                     :   taxonomy ID of taxa
        "taxName"                   :   taxonomy name
        "reads_nb"                  :   number of reads affiliated to this taxa
        "taxo"                      :   complete taxonomy in a string    
    """
    if not osp.isfile(krep_fp):
        raise RuntimeError(f"ERROR: File '{krep_fp}' not found.")

    ### load kraken-report files and parse it
    id_sample, names, sample_counts, taxa = process_sample(krep_fp, max_rank=maxi_rank, min_rank=mini_rank)
    taxID = list(taxa.keys())

    taxa_str=[]
    for key, value in taxa.items():
        taxa_sample = '_'.join(level for level in value)
        taxa_sample = taxa_sample.replace(',', '').replace("'", '').replace('"', '')
        taxa_str.append(taxa_sample)

    ## remove spaces in excess from taxonomy names 
    names_truncated = [name.strip() for name in names.values()]
    df = pd.DataFrame({'sampleID':id_sample,'taxID': taxID, 'taxName': names_truncated, 'reads_nb': pd.Series(sample_counts), 'taxo': taxa_str})

    return df

def process_sample(kraken_report, max_rank, min_rank):
    """
    Parse all kraken-report data files into sample counts dict
    and store global taxon id -> taxonomy data
    """
    taxa = OrderedDict()
    sample_counts = OrderedDict()
    if not osp.isfile(kraken_report):
        raise RuntimeError("ERROR: File '{}' not found.".format(kraken_report))

    # use the kraken report filename as the sample ID
    sample_id = osp.splitext(osp.split(kraken_report)[1])[0]

    with open(kraken_report, "rt") as kf:
        try:
            kdr = csv.DictReader(kf, fieldnames=field_names, 
                                    delimiter="\t")
            kdata = [entry for entry in kdr][1:]
        except OSError as oe:
            raise RuntimeError("ERROR: {}".format(oe))

    names, scounts, staxa = parse_kraken_report(kdata, max_rank=max_rank, min_rank=min_rank)

    ### 24/07/23 : staxa is a dictonnary made to avoid transforming the taxa in string then to dictonnary agai

    # update master records
    taxa.update(staxa)

    return sample_id, names, scounts, taxa

def parse_kraken_report(kdata, max_rank, min_rank):
    """
    Parse a single output file from the kraken-report tool. Return a list
    of counts at each of the acceptable taxonomic levels, and a list of 
    NCBI IDs and a formatted string representing their taxonomic hierarchies.
    :type kdata: str
    :param kdata: Contents of the kraken report file.
    """
    # map between NCBI taxonomy IDs and the string rep. of the hierarchy
    taxa = OrderedDict()
    # the master collection of read counts (keyed on NCBI ID)
    counts = OrderedDict()
    # the names of taxas (added)
    names = OrderedDict()
    # current rank
    r = 0
    max_rank_idx = ranks.index(max_rank)
    min_rank_idx = ranks.index(min_rank)

    for entry in kdata:
        # update running tally of ranks
        tax_lvl = parse_tax_lvl(entry)

        erank = entry['rank'].strip()
        if 'SS' in tax_lvl:
            erank = 'SS'

        if erank in ranks:
            r = ranks.index(erank)     

        # record the reads assigned to this taxon level, and record the taxonomy string with the NCBI ID
        if erank in ranks and min_rank_idx >= ranks.index(erank) >= max_rank_idx:
            taxon_reads = int(entry["taxon_reads"])
            clade_reads = int(entry["clade_reads"])
            if taxon_reads > 0 or (clade_reads > 0 and erank == min_rank):
                taxa[entry['ncbi_tax']] = tax_fmt(tax_lvl, r)
                names[entry['ncbi_tax']] = entry["sci_name"]

                if erank == min_rank:
                    counts[entry['ncbi_tax']] = clade_reads
                else:
                    counts[entry['ncbi_tax']] = taxon_reads

    return names, counts, taxa

def parse_tax_lvl(entry, tax_lvl_depth=[]):
    """
    Parse a single kraken-report entry and return a dictionary of taxa for its
    named ranks.
    :type entry: dict
    :param entry: attributes of a single kraken-report row.
    :type tax_lvl_depth: list
    :param tax_lvl_depth: running record of taxon levels encountered in
    previous calls.
    """
    # How deep in the hierarchy are we currently?  Each two spaces of
    # indentation is one level deeper.  Also parse the scientific name at this
    # level.
    depth_and_name = re.match('^( *)(.*)', str(entry['sci_name']))
    depth = len(depth_and_name.group(1))//2
    name = depth_and_name.group(2)
    # Remove the previous levels so we're one higher than the level of the new
    # taxon.  (This also works if we're just starting out or are going deeper.)
    del tax_lvl_depth[depth:]
    # Append the new taxon.
    erank = entry['rank']
    if erank == '-' and depth > 8 and tax_lvl_depth[-1][0] == 'S':
        erank = 'SS'
    tax_lvl_depth.append((erank, name))

    # Create a tax_lvl dict for the named ranks.
    tax_lvl = {x[0]: x[1] for x in tax_lvl_depth if x[0] in ranks}
    return(tax_lvl)

def tax_fmt(tax_lvl, end):
    """
    Create a string representation of a taxonomic hierarchy (QIIME for now).
    :type tax_lvl: dict
    :param tax_lvl: Keyed on the entries in ranks
    :type end: int
    :param end: The end rank index (0-based indexing). If end == 3
                then the returned string will contain K, P, C, and O.
    
    >>> tax_fmt({"K": "Bacteria", "P": "Firmicutes", "C": "Negativicutes", 
    ... "O": "Selenomonadales", "F": "Veillonellaceae", "G": "Veillonella", 
    ... "S": "Veillonella parvula", "SS": "Veillonella parvula DSM 2008"}, 4)
    'k__Bacteria; p__Firmicutes; c__Negativicutes; o__Selenomonadales'
    """
    if "S" in tax_lvl:
        if "G" in tax_lvl and tax_lvl["S"].startswith(tax_lvl["G"]):
            tax_lvl["S"] = tax_lvl["S"][len(tax_lvl["G"])+1:]
    if "SS" in tax_lvl:
        if "S" in tax_lvl and tax_lvl["SS"].startswith(tax_lvl["S"]):
            tax_lvl["SS"] = tax_lvl["SS"][len(tax_lvl["S"])+1:]
    
    tax = ["{}__{}".format(r.lower(), tax_lvl[r] if r in tax_lvl else '') 
             for r in ranks[:end+1]]
    # add empty identifiers for ranks beyond end
    #TODO: remove the :-1 when SS support is added
    tax.extend(["{}__".format(r.lower()) for r in ranks[end+1:-1]])

    # even though Bacteria, Archea are now technically Domains/superkingdoms
    # GreenGenes and other databases still list the traditional 
    # kingdom/phylum/class/etc. So this is a hack to shoehorn the kraken-report
    # data into that format.
    if tax[0].startswith('d'):
        tax[0] = "k"+tax[0][1:]

    return tax

def extract_main_taxa(df):
    """From a df, extracts and returns the row having the biggest value in 'reads_nb' column"""
    df = df.sort_values(by = 'reads_nb', ascending=False)
    main_taxa = df.iloc[0,:]
    return main_taxa

def convert_tax(taxo):
    """
    Input : 
        taxo : string (example : "k__Bacteria_p__Cyanobacteria_c___o__Synechococcales_f__Prochlorococcaceae_g__Prochlorococcus_s__marinus")
    Converts the string in a dictionnary arranged for taxonomy 
    Output :
        dicotaxo : {'k': 'Bacteria', 
                    'p': 'Cyanobacteria', 
                    'c': '', 
                    'o': 'Synechococcales', 
                    'f': 'Prochlorococcaceae', 
                    'g': 'Prochlorococcus', 
                    's': 'marinus'}
    """
    phylo_levels = ["k", "p", "c", "o", "f", "g", "s"]

    dicotaxo = {}
    for i in range(len(phylo_levels)):          # for each letter of phylo_levels...
        level = phylo_levels[i]
        pattern = f"{level}__(.*?)" + (f"_{phylo_levels[i+1]}" if i < len(phylo_levels)-1 else "$")
        
        match = re.search(pattern, taxo)        # the string following the letter is saved in the dictionnary 
        if match:
            dicotaxo[level] = match.group(1)
        else:
            dicotaxo[level] = ""

    return dicotaxo

def sample_vs_outliers(taxa_sample, taxas_outliers):
    """ 
    Input : 
        taxa_sample : row of df of main taxa from sample with columns ['taxID', 'taxName', 'reads_nb', 'taxo']
        taxas_outliers : df with columns ['taxID', 'taxName', 'reads_nb', 'taxo']
    Compares for each row of the df the taxonomy with the row of main taxa from sample : first at family level, then at 
    order level if family is empty, then at class level if order is empty. Counts contaminants when not matching with the 
    taxa from main sample. 
    Output :
        taxas_found : dictionnary {taxa : occurrency}
        contaminants : int of contaminants found 
    """
    taxas_found = {}
    contaminants = 0
    rank = "f"      # first rank of reference is the family                                    
    for index, row in taxas_outliers.iterrows():            # for each outlier found...
        taxa = row['taxo'][rank]                            # ... the taxa of the same rank as the main taxa is checked
        if taxa != taxa_sample['taxo'][rank]:               # if the taxa searched doesn't match the main taxa...
            if taxa == '':                                      # ...if it's empty at family level...
                new_rank = "o"     
                new_taxa = row['taxo'][new_rank]                    # ...a match is checked at order level
                if new_taxa != taxa_sample['taxo'][new_rank]:       # if it still doesn't match...
                    if new_taxa== "":                                   # ... because it's empty at order level...
                        new_rank = "c"                                      
                        new_taxa = row['taxo'][new_rank]                    # ...a match is checked at class level
                        if new_taxa != taxa_sample['taxo'][new_rank]:       # if it still doesn't match...
                            if new_taxa== "":                                   # ... because it's empty at class level...
                                continue                                        # ... its taxonomy is too blurred to be counted
                            contaminants += 1                               
                            if new_taxa not in taxas_found.keys():          # and if it's the first encounter...
                                taxas_found[new_taxa] = 1                       # ... it's added in the dictionnary
                            else : 
                                taxas_found[new_taxa] += 1                  # else its occurrency is incremented             
                    elif new_taxa not in taxas_found.keys():          # and if it's the first encounter...
                        #print(taxas_found.keys())
                        taxas_found[new_taxa] = 1                       # .... it's added in the dictionnary
                        contaminants += 1 
                    else : 
                        taxas_found[new_taxa] += 1                  # else its occurrency is incremented
                        contaminants += 1 
            elif taxa not in taxas_found.keys():            # if family is different and it's the first encounter...
                contaminants += 1
                taxas_found[taxa] = 1                           # ...it's added in the dictionnary       
            else : 
                taxas_found[taxa] += 1                          # else its occurrency is incremented

    return taxas_found, contaminants
