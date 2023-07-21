# Run padloc
rule run_padloc:
    output:
        dir = directory("results/intermediate/padloc/padloc_{sample}"),
        csv = "results/intermediate/padloc/padloc_{sample}/{sample}_padloc.csv"
    input:
        faa = "data/protein_seq/{sample}.faa",
        gff = "data/annotation/{sample}.gff",
        hmm = PADLOC_DB_DIR + "/hmm"
    conda: "../envs/padloc.yaml"
    threads: 8
    params: pdd = PADLOC_DB_DIR
    shell:
        r"""touch {output.csv} #empty output file in case padloc finds no results
                               #otherwise snakemake exits due to lack of output file
            padloc --force --data {params.pdd} --faa {input.faa} --gff {input.gff} --outdir {output.dir} --cpu {threads}
         """
