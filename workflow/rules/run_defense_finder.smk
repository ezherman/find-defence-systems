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
    threads: 8
    shell:
        r"""outdir=$(dirname {output.systems})
            rm -rf $outdir #if the outdir exists, remove it
            gunzip -c {input.data} > {output.faa}
            defense-finder run {output.faa} --out-dir $outdir -w {threads}
            rename "{wildcards.sample}_" "" $outdir/* #defense finder includes sample in the filenames
         """
