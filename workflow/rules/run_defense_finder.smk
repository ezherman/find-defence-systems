# Run defense finder
# defense-finder creates a folder inside out-dir
# this rule moves files from that folder one directory up
# after running defense-finder
rule run_defense_finder:
    output:
        dir = directory("results/intermediate/defense_finder/defense_finder_{sample}"),
        tsv = "results/intermediate/defense_finder/defense_finder_{sample}/defense_finder_systems.tsv"
    input: 
        data = "data/protein_seq/{sample}.faa",
        models = os.path.expanduser('~') + "/.macsyfinder/data"
    conda: "../envs/defensefinder.yaml"
    shell:
        r"""defense-finder run {input.data} --out-dir {output.dir} -w 1
            rm {input.data}.idx #remove index file
         """
