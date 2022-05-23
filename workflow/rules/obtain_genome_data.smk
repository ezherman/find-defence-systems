# Download genome data from NCBI
rule obtain_genome_data:
    output:
        prt_sq  = "data/protein_seq/{sample}.faa",
        annot   = "data/annotation/{sample}.gff"
    run:

        # obtain genome data for each strain
        download_seq_files(wildcards.sample)
        
        # unzip the files
        for f in [output.prt_sq]:
            with gzip.open(f + '.gz', 'r') as f_in, open(f, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        for f in [output.annot]:
            with gzip.open(f + '.gz', 'r') as f_in, open(f, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)


        


        



