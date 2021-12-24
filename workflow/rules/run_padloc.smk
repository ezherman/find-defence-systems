# Run padloc
rule run_padloc:
    output:
        dir = directory("data/intermediate_data/padloc/padloc_{assembly}"),
        csv = "data/intermediate_data/padloc/padloc_{assembly}/{assembly}_padloc.csv"
    input:
        faa = "data/input_data/protein_sequences/{assembly}.faa",
        gff = "data/input_data/annotations/{assembly}.gff",
        hmm = ".snakemake/conda/" + find_conda_env_hash("scripts/env/padloc.yaml") + "/data/hmm"
    conda: "env/padloc.yaml"
    shell:
        r"""touch {output.csv} #empty output file in case padloc finds no results
                               #otherwise snakemake exits due to lack of output file
            padloc --force --faa {input.faa} --gff {input.gff} --outdir {output.dir}
         """
