# Update padloc database
# find_conda_env_hash is used because the location of the
# padloc env directory that Snakemake creates includes a hash key that is
# dependent on the user's machine
rule update_padloc_db:
    output: directory(PADLOC_DB_DIR + "/hmm")
    conda: "../envs/padloc.yaml"
    params: pdd = PADLOC_DB_DIR
    shell:
        "padloc --data {params.pdd} --db-update"
