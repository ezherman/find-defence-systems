# Run padloc
rule run_padloc:
    output:
        dir = directory("results/intermediate/padloc/padloc_{sample}"),
        csv = "results/intermediate/padloc/padloc_{sample}/{sample}_padloc.csv"
    input:
        faa = "data/protein_seq/{sample}.faa",
        gff = "data/annotation/{sample}.gff",
        hmm = ".snakemake/conda/" + find_conda_env_hash("workflow/envs/padloc.yaml") + "/data/hmm"
    conda: "../envs/padloc.yaml"
    shell:
        r"""touch {output.csv} #empty output file in case padloc finds no results
                               #otherwise snakemake exits due to lack of output file
            padloc --force --faa {input.faa} --gff {input.gff} --outdir {output.dir}
         """
