# Run padloc
rule run_padloc:
    output:
        dir     = directory("results/intermediate/padloc/padloc_{sample}"),
        csv     = "results/intermediate/padloc/padloc_{sample}/{sample}_padloc.csv",
        faa     = temp("data/protein_seq/{sample}.faa"),
        gff     = temp("data/annotation/{sample}.gff")
    input:
        faa_gz  = "data/protein_seq/{sample}.faa.gz",
        gff_gz  = "data/annotation/{sample}.gff.gz",
        hmm = PADLOC_DB_DIR + "/hmm"
    conda: "../envs/padloc.yaml"
    threads: 8
    params: pdd = PADLOC_DB_DIR
    shell:
        r"""touch {output.csv} #empty output file in case padloc finds no results
                               #otherwise snakemake exits due to lack of output file
            gunzip -c {input.faa_gz} > {output.faa}
            gunzip -c {input.gff_gz} > {output.gff}
            padloc --force --data {params.pdd} --faa {output.faa} --gff {output.gff} --outdir {output.dir} --cpu {threads}
         """
