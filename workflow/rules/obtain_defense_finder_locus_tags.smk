rule obtain_defense_finder_locus_tags:
    output:
        hits            = "results/intermediate/defense_finder/defense_finder_{sample}/defense_finder_genes_ltags.csv",
        single_line_faa = temp("single_line_{sample}.faa")
    input:
        defense_finder  = "results/intermediate/defense_finder/defense_finder_{sample}/defense_finder_genes.tsv",
        gbff_gz         = "data/annotation/{sample}.gbff.gz",
        faa_gz          = "data/protein_seq/{sample}.faa.gz"
    run:
        # check if defensefinder found any hits
        # if not, an empty spreadsheet can be returned
        hits = pd.read_table(input.defense_finder)
        if len(hits) == 0:
            hits['locus_tag'] = []
            hits.to_csv(output.hits, index = False)
            shell('touch {output.single_line_faa}') #snakemake requires this file to complete the job
            

        else:
            # check if annotation was performed with Bakta
            # if so, the target.name column is already the locus tag
            with gzip.open(input.gbff_gz, "r") as f:
                lines = f.readlines()
                bakta_true = "Bakta" in str(lines[8])

            if bakta_true:
                hits['locus_tag'] = hits['hit_id']
                hits.to_csv(output.hits, index = False)
                shell('touch {output.single_line_faa}') #snakemake requires this file to complete the job 
        
            else:
        
                #-------- load data
                
                # parse the gbff file into a list
                with gzip.open(input.gbff_gz, "rt") as f:
                    gbff_out = list(SeqIO.parse(f, format="genbank"))
    
                # convert the multiline fasta into single-line
                shell("gunzip -c {input.faa_gz} | seqtk seq -l0 > {output.single_line_faa}")
    
                # load the defense finder hits into a dictionary
                hits = pd.read_table(input.defense_finder)
                unique_hits = hits.drop_duplicates(subset = ['hit_id', 'hit_pos'])
    
                hit_dicts = [
                    {
                        "hit_id": hit_id,
                        "hit_pos": hit_pos,
                    }
                    for hit_id, hit_pos in zip(
                        unique_hits["hit_id"], unique_hits["hit_pos"]
                    )
                ]
    
                #-------- add locus tags to hit_dicts using the hit sequences
    
                # first construct a dict of sequences to locus tags
                locus_tag_dict = {}
                for contig in gbff_out:
                    for feature in contig.features:
                        # not all CDS features have an associated translation, so need to include
                        # that in the if statement
                        if feature.type == "CDS" and "translation" in feature.qualifiers.keys():
                            locus_tag_dict.update(
                                {
                                    feature.qualifiers["translation"][0]: feature.qualifiers[
                                        "locus_tag"
                                    ][0]
                                }
                            )
                
                # then add locus_tag values of interest by matching the protein sequence
                with open(output.single_line_faa, "r") as f:
                    sequences = [l for l in f.readlines()]
                    for entry in hit_dicts:
                        # hit_pos * 2 - 1 because there are two lines per hit: fasta header and sequence
                        # and the first entry has index 0 rather than 1
                        sequence = sequences[entry["hit_pos"] * 2 - 1].strip("\n")
                        locus_tag = locus_tag_dict[sequence]
                        entry.update({"locus_tag": locus_tag})
    
                # -------- merge locus tags onto the defensefinder output table
                pd_out = hits.merge(pd.DataFrame.from_records(hit_dicts), on = ['hit_id', 'hit_pos'])
                pd_out.to_csv(output.hits, index = False)

