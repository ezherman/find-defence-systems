# Update defense-finder models
# By default, the models are downloaded to $HOME/.macsyfinder/data
rule update_defense_finder_models:
    output: directory(os.path.expanduser('~') + "/.macsyfinder/data")
    conda: "../envs/defensefinder.yaml"
    shell:
        r"""mkdir -p {output}
            macsydata install -U --org  mdmparis -f defense-finder-models
        """