from snakemake.utils import validate
import pandas as pd

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
# singularity: "docker://continuumio/miniconda3"
containerized: "docker://ezherman/find-defence-systems"

##### load config and sample sheets #####

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t")
SAMPLE_NAMES = samples["sample"]
PADLOC_DB_DIR = config["padloc_db_dir"]
validate(samples, schema="../schemas/samples.schema.yaml")


################ define helper functions
import pandas as pd
import os
import glob
import gzip, shutil
import gffutils
import roman
import wget #installed with pip
from Bio import Entrez, SeqIO

## functions and class to download sequence files

# class to create combinations of data type and file extension
class ncbi_combination():
    def __init__(self, data_type, extension):
        self.data_type = data_type
        self.extension = extension

allowed_ncbi_combinations = [
    # from https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/147/155/GCF_900147155.1_390/ 
    # as an example
    ncbi_combination('assembly_report', 'txt'),
    ncbi_combination('assembly_stats', 'txt'),
    ncbi_combination('cds_from_genomic', 'fna.gz'),
    ncbi_combination('feature_count', 'txt.gz'),
    ncbi_combination('feature_table', 'txt.gz'),
    ncbi_combination('genomic', 'fna.gz'),
    ncbi_combination('genomic', 'gbff.gz'),
    ncbi_combination('genomic', 'gff.gz'),
    ncbi_combination('genomic', 'gtf.gz'),
    ncbi_combination('protein', 'faa.gz'),
    ncbi_combination('protein', 'gpff.gz'),
    ncbi_combination('rna_from_genomic', 'fna.gz'),
    ncbi_combination('translated_cds', 'faa.gz'),
    ncbi_combination('wgsmaster', 'gbff.gz'),
]

def seq_file_downloader(refseq_accession: str, combination: ncbi_combination, outdir: str):
    """Download a file from RefSeq given a GCF accession, 
       combination of data type and extension and an 
       output directory."""

    refseq_id = Entrez.read(Entrez.esearch(db = 'assembly', term = refseq_accession))['IdList'][0]
    ftp_base_path = Entrez.read(Entrez.esummary(db="assembly", id=refseq_id))["DocumentSummarySet"][
        "DocumentSummary"
    ][0]["FtpPath_RefSeq"]
    ftp_accession = ftp_base_path.split('/')[-1]
    ftp_path = ftp_base_path + '/' + ftp_accession + '_' + combination.data_type + '.' + combination.extension

    wget.download(ftp_path, out = outdir + refseq_accession + '.' + combination.extension)


def download_seq_file(refseq_accession: str, data_type: str, extension: str, outdir: str):
    """Check if the combination of data_type and str is allowed.
       If it is, download the file with seq_file_downloader."""

    combination = ncbi_combination(data_type, extension)
    if vars(combination) in [vars(c) for c in allowed_ncbi_combinations]:
        # download file
        seq_file_downloader(refseq_accession, combination, outdir)
    else:
        raise ValueError(
            'The provided combination of data_type and extension is not supported.'
            )

# Define a custom function for generalisation of system names
# between padloc and defense_finder
def generalise_system_names(row, software):

    system = row.system.lower()

    if software == 'defense_finder':

        if 'cas_' in system:
            # remove class from name
            system = re.sub('class[0-9]+-subtype-', 'subtype_', system)
        

        # zorya name misses and underscore, e.g. zorya_typei instead of zorya_type_i
        if 'zorya' in system:
            system = '_'.join(['zorya_type', system.split('zorya_type')[-1]])
        
        # old is called old_exonuclease here and simply old in padloc
        if system == 'old_exonuclease':
            system = 'old'
        
        # retron subsystems are delimited with '_', e.g. retron_i_b
        # in padloc they're delimited with '-', e.g. retron_i-b
        # set defensefinder to padloc naming
        if 'retron' in system:
            system_split = system.split('_') 
            subtype = system_split[-1] 
            system = '_'.join(system_split[0:-1]) + '-' + subtype
    
    elif software == 'padloc':

        if 'cas_' in system:

            # if the naming matches 'type_[a-z]+-',
            # where the final '-' is key,
            # it's actually a cas subtype definition
            system = re.sub('type_[a-z]+-', 
                            'subtype_' + system.split('_')[-1].split('-')[0] + '-',
                            system)
            
            # for some subtypes, padloc has multiple models (e.g. i-f1, i-f2, i-f3)
            # the defensefinder definitions aren't as granular
            # generalise the padloc naming by removing the numerical suffix
            system = re.sub('([a-z]+-[a-z])[0-9]', '\\1', system)

        # dsr systems miss the '_' separator
        if 'dsr' in system and '_' not in system:
            system = '_'.join(['dsr', system.split('dsr')[-1]])

    # remove type from name
    if 'cas_' not in system:
        system = re.sub('_type', '', system)
    
    # if numbering is integer, convert to roman
    if system.split('_')[-1].isdigit():
        roman_numeral = roman.toRoman(int(system.split('_')[-1])).lower()
        system = '_'.join(system.split('_')[0:-1]) + '_' + roman_numeral

    return system

## merge padloc and defense_finder dataframes
## if padloc is empty because no results were found
## process only defense_finder
def create_subsystem_table(df, program):

    if program == "padloc":
        df = (
            df
            .groupby(['system.number','system'])                     # group by system
            .agg(lambda x: ';'.join(map(str, x)))    # collapse rows within systems
            .reset_index()                           # ungroup
            .rename(columns={"protein.name":"protein_names", "target.name":"protein_IDs", "locus_tag":"locus_tags"})) # rename columns
        
        df = df[["system", "protein_names", "protein_IDs", "locus_tags"]]   # choose columns

    if program == "defense_finder":
        df = (
            df
            .groupby(['sys_id', 'system'])                     # group by system
            .agg(lambda x: ';'.join(map(str, x)))    # collapse rows within systems
            .rename(columns={"gene_name":"protein_names", "hit_id":"protein_IDs", "locus_tag":"locus_tags"})
            .reset_index()
            )
        df = df[["system", "protein_names", "protein_IDs", "locus_tags"]] # choose columns

    return df

## merge padloc and defense_finder dataframes
## if padloc is empty because no results were found
## process only defense_finder
def merge_subsystem_tables(df_defense_finder, df_padloc = None):

    # if no df_padloc exists, final table is defense_finder only
    if df_padloc is None:

        df = pd.concat([df_defense_finder],
                       keys = ["defense_finder"],
                       names = ("program", "row"))

        # remove the row index, then make the program index into a column
        df = df.reset_index(level = 'row', drop = True).reset_index()

    # if df_padloc exists, merge padloc and defense_finder
    if isinstance(df_padloc, pd.DataFrame):

        df = pd.concat([df_padloc, df_defense_finder],
                   keys = ["padloc", "defense_finder"],
                   names = ("program", "row"))          # merge tables

        # remove the row index, then make the program index into a column
        df = df.reset_index(level = 'row', drop = True).reset_index()


    return df


def matrix_all_samples(files: list[str], binary: bool) -> pd.DataFrame:
    """Merge defence system tables
       Also reduce data from number of hits to presence-absence if binary = True
    """

    # concatenate all analysis output
    df = pd.concat((pd.read_csv(f) for f in files),
                   keys = SAMPLE_NAMES,
                   names = ("sample", "row"))

    #remove row index
    df = df.reset_index(level = "row", drop = True)


    # add samples without defense system hits as "None"
    # if a file in files is empty, the corresponding key
    # is skipped and not added to the keys of df
    no_hit_samples = [a for a in SAMPLE_NAMES if a not in df.index]
    new_entry = {'program': 'NA', 'system': 'none', 'protein_names': 'NA', 'protein_IDs' : 'NA'}
    for a in no_hit_samples:
        df.loc[a] = new_entry

    # create wide table with systems as columns and samples as rows
    table = pd.crosstab(df.index, df.system)

    if binary:
        # reduce data to presence-absence
        for col in table.columns:
            table[col] = [1 if v != 0 else 0 for v in table[col]]

    # rename index
    table = table.rename_axis('sample')

    return table


