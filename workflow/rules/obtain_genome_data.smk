# Download genome data from NCBI
rule obtain_genome_data:
    output:
        prt_sq  = "data/protein_seq/{sample}.faa",
        annot   = "data/annotation/{sample}.gff"
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

        
        
        # unzip the files
        for f in output:
            with gzip.open(f + '.gz', 'r') as f_in, open(f, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)


        


        



