# Update padloc database
# find_conda_env_hash is used because the location of the
# padloc env directory that Snakemake creates includes a hash key that is
# dependent on the user's machine
rule update_padloc_db:
    output: directory(".snakemake/conda/" + find_conda_env_hash("scripts/env/padloc.yaml") + "/data/hmm")
    conda: "env/padloc.yaml"
    shell:
        "padloc --db-update"
