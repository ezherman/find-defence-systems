# Run defense finder
rule run_defense_finder:
    output:
        systems = "results/intermediate/defense_finder/defense_finder_{sample}/defense_finder_systems.tsv",
        genes   = "results/intermediate/defense_finder/defense_finder_{sample}/defense_finder_genes.tsv",
        index   = temp("tmp_unzipped_defense_finder_inputs/{sample}.faa.idx"),
        faa     = temp("tmp_unzipped_defense_finder_inputs/{sample}.faa")
    input: 
        data = "data/protein_seq/{sample}.faa.gz",
        models = os.path.expanduser('~') + "/.macsyfinder/data"
    conda: "../envs/defensefinder.yaml"
    params: outdir_prefix = "results/intermediate/defense_finder/defense_finder"
    threads: 8
    shell:
        r"""rm -rf {params.outdir_prefix}_{wildcards.sample} #if the outdir exists, remove it
            gunzip -c {input.data} > {output.faa}
            defense-finder run {output.faa} --out-dir {params.outdir_prefix}_{wildcards.sample} -w {threads}
         """
