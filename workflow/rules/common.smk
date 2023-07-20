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
import wget #installed with pip

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

def seq_file_downloader(gcf: str, combination: ncbi_combination, outdir: str):
    """Download a file from RefSeq given a GCF accession, 
       combination of data type and extension and an 
       output directory."""

    base = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/"
    folder = base + gcf[4:7] + "/" + gcf[7:10] + "/" + gcf[10:13] + "/" + gcf + "/" 
    file = folder + gcf + '_' + combination.data_type + '.' + combination.extension

    wget.download(file, out = outdir + gcf + '.' + combination.extension)


def download_seq_file(gcf: str, data_type: str, extension: str, outdir: str):
    """Check if the combination of data_type and str is allowed.
       If it is, download the file with seq_file_downloader."""

    combination = ncbi_combination(data_type, extension)
    if vars(combination) in [vars(c) for c in allowed_ncbi_combinations]:
        # download file
        seq_file_downloader(gcf, combination, outdir)
    else:
        raise ValueError(
            'The provided combination of data_type and extension is not supported.'
            )

## wrangle pandas dataframe for padloc and defense_finder
def create_subsystem_table(df, program):

    if program == "padloc":
        df = df.groupby(['system.number',
                         'system'])                     # group by system
        df = df.agg(lambda x: ';'.join(map(str, x)))    # collapse rows within systems
        df = df.reset_index()                           # ungroup
        df = df.rename(columns={"protein.name":"protein_names",
                                "target.name":"protein_IDs"}) # rename columns
        df = df[["system", "protein_names", "protein_IDs"]]   # choose columns

    if program == "defense_finder":
        df = df.groupby(['sys_id'])                     # group by system
        df = df.agg(lambda x: ';'.join(map(str, x)))    # collapse rows within systems
        df = df.reset_index()                           # ungroup

        #rename columns
        df = df.rename(columns={"sys_id":"system",
                                "name_of_profiles_in_sys":"protein_names",
                                "protein_in_syst":"protein_IDs"})
        df = df[["system", "protein_names", "protein_IDs"]] # choose columns

        #remove "UserReplicon_" from system names
        df["system"] = df["system"].str.replace('UserReplicon_', '')

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

        df = df.reset_index(level = "row", drop = True)

    # if df_padloc exists, merge padloc and defense_finder
    if isinstance(df_padloc, pd.DataFrame):

        df = pd.concat([df_padloc, df_defense_finder],
                   keys = ["padloc", "defense_finder"],
                   names = ("program", "row"))          # merge tables

        df = df.reset_index(level = "row", drop = True)     # drop row index


    return df


def matrix_all_samples(files: str, binary: bool) -> pd.DataFrame:
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


