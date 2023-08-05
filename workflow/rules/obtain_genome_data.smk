# Download genome data from NCBI
rule obtain_genome_data:
    output:
        faa_gz  = "data/protein_seq/{sample}.faa.gz",
        gff_gz  = "data/annotation/{sample}.gff.gz",
        faa     = temp("data/protein_seq/{sample}.faa"),
        gff     = temp("data/annotation/{sample}.gff")
    retries: 2
    resources: ncbi_connection = 1
    run:

        # obtain genome data
        for out in output:

            outdir = '/'.join(out.split('/')[:-1]) + '/'
            
            if 'protein_seq' in out:
                data_type = 'protein'
                extension = 'faa.gz'
            else:
                data_type = 'genomic'
                extension = 'gff.gz'
            
            download_seq_file(wildcards.sample, data_type, extension, outdir)

        
        
        # sometimes a corrupt file is obtained, which results in an unzip error in ppanggolin.
        # try to gunzip here. 
        # if bash throws an error, the rule is retried -> files are downloaded again.
        shell("gunzip -c {output.fna_gz} > {output.fna}")
        shell("gunzip -c {output.gff_gz} > {output.gff}")
            


        


        



