# Run padloc
rule run_padloc:
    output:
        csv     = "results/intermediate/padloc/padloc_{sample}/{sample}_padloc.csv",
        faa     = temp("tmp_unzipped_padloc_inputs/{sample}.faa"),
        gff     = temp("tmp_unzipped_padloc_inputs/{sample}.gff")
    input:
        faa_gz  = "data/protein_seq/{sample}.faa.gz",
        gff_gz  = "data/annotation/{sample}.gff.gz",
        hmm = PADLOC_DB_DIR + "/hmm"
    conda: "../envs/padloc.yaml"
    threads: 8
    params: 
        pdd             = PADLOC_DB_DIR,
        outdir_prefix   = "results/intermediate/padloc/padloc"
    shell:
        r"""#empty output file in case padloc finds no results
            #otherwise snakemake exits due to lack of output file
            #and downstream processing fails
            echo "system.number,seqid,system,target.name,hmm.accession,hmm.name,protein.name,full.seq.E.value,domain.iE.value,target.coverage,hmm.coverage,start,end,strand,target.description,relative.position,contig.end,all.domains,best.hits" > {output.csv}

            gunzip -c {input.faa_gz} > {output.faa}
            gunzip -c {input.gff_gz} > {output.gff}
            padloc --force --data {params.pdd} --faa {output.faa} --gff {output.gff} --outdir {params.outdir_prefix}_{wildcards.sample} --cpu {threads}
         """
