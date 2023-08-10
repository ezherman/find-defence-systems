# Download genome data from NCBI
rule obtain_genome_data:
    output:
        gz_file = "data/{file_type}/{sample}.{ext}.gz",
    retries: 2
    resources: ncbi_connection = 1
    run:

        # in theory other data types can be obtained from ncbi
        # these are not supported at the moment here
        if wildcards.file_type == 'protein_seq':
            data_type = 'protein'
        elif wildcards.file_type in ['annotation', 'assembly']:
            data_type = 'genomic'
        else:
            ValueError('Unsuppoted file_type wildcard')

        outdir = '/'.join(output.gz_file.split('/')[:-1]) + '/'
            
        download_seq_file(wildcards.sample, data_type, wildcards.ext + '.gz', outdir)
        
        # sometimes a corrupt file is obtained, which results in downstream unzipping errors
        # try to gunzip here. 
        # if bash throws an error, the rule is retried -> files are downloaded again.
        shell("gunzip -c {output.gz_file} > tmp.{wildcards.sample}.{wildcards.ext}; rm -f tmp.{wildcards.sample}.{wildcards.ext}")
            


        


        



