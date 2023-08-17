# Run padloc
rule run_padloc:
    output:
        csv     = "results/intermediate/padloc/padloc_{sample}/{sample}_padloc.csv",
        faa     = temp("tmp_unzipped_padloc_inputs/{sample}.faa"),
        gff     = temp("tmp_unzipped_padloc_inputs/{sample}.gff")
    input:
        faa_gz  = "data/protein_seq/{sample}.faa.gz",
        gff_gz  = "data/annotation/{sample}.gff.gz",
        hmm = PADLOC_DB_DIR + "/hmm",
        empty_padloc_out = 'resources/padloc/empty_padloc_output.csv'
    conda: "../envs/padloc.yaml"
    threads: 8
    params: 
        pdd             = PADLOC_DB_DIR,
        outdir_prefix   = "results/intermediate/padloc/padloc"
    shell:
        r"""cp {input.empty_padloc_out} {output.csv} #empty output file in case padloc finds no results
                                                     #otherwise snakemake exits due to lack of output file
            gunzip -c {input.faa_gz} > {output.faa}
            gunzip -c {input.gff_gz} > {output.gff}
            padloc --force --data {params.pdd} --faa {output.faa} --gff {output.gff} --outdir {params.outdir_prefix}_{wildcards.sample} --cpu {threads}
         """
