# Run defense finder
# defense-finder creates a folder inside out-dir
# this rule moves files from that folder one directory up
# after running defense-finder
rule run_defense_finder:
    output:
        dir = directory("results/intermediate/defense_finder/defense_finder_{assembly}"),
        tsv = "results/intermediate/defense_finder/defense_finder_{assembly}/defense_finder_systems.tsv"
    input: "data/protein_seq/{assembly}.faa"
    conda: "../envs/defensefinder.yaml"
    resources:
        n_defense_finder = 1
    shell:
        r"""defense-finder run {input} --out-dir {output.dir}
            mv {output.dir}/*-defense-finder/* {output.dir}
            rm -r {output.dir}/*-defense-finder
            rm {input}.idx #remove index file
         """
