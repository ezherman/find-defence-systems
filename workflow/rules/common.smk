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
validate(samples, schema="../schemas/samples.schema.yaml")


################ define helper functions
import pandas as pd
import hashlib
import os
import glob
import gzip, shutil
import wget #installed with pip

## function to download sequence files
## strain_to_sra created in obtain_genome_data rule

def download_seq_files(gcf):
    base = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/"
    folder = base + gcf[4:7] + "/" + gcf[7:10] + "/" + gcf[10:13] + "/" + gcf + "/" 
    gff = folder + gcf + "_genomic.gff.gz"
    faa = folder + gcf + "_protein.faa.gz"

    # download gff and faa files
    wget.download(gff, out = "data/annotation/" + gcf + ".gff.gz")
    wget.download(faa, out = "data/protein_seq/" + gcf + ".faa.gz")

## find conda environment hash, adapted from
## https://bioinformatics.stackexchange.com/a/9346
def find_conda_env_hash(yaml):
    md5hash = hashlib.md5()
    md5hash.update((os.getcwd() + "/.snakemake/conda").encode())
    f = open(yaml, 'rb')
    md5hash.update(f.read())
    f.close()
    h = md5hash.hexdigest()
    return h

## wrangle pandas dataframe for padloc and defense_finder
def create_subsystem_table(df, program):

    if program == "padloc":
        df = df.groupby(['system.number',
                         'system'])         # group by system
        df = df.agg(';'.join)               # collapse rows within systems
        df = df.reset_index()               # ungroup
        df = df.rename(columns={"protein.name":"protein_names",
                                "target.name":"protein_IDs"}) # rename columns
        df = df[["system", "protein_names", "protein_IDs"]]   # choose columns

    if program == "defense_finder":
        df = df.groupby(['sys_id'])         # group by system
        df = df.agg(';'.join)               # collapse rows within systems
        df = df.reset_index()               # ungroup

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


## merge system tables of all samples
## also reduce data from number of hits to presence-absence
def matrix_all_samples(files):

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

    # reduce data to presence-absence
    for col in table.columns:
        table[col] = [1 if v != 0 else 0 for v in table[col]]

    # rename index
    table = table.rename_axis('sample')

    return table


