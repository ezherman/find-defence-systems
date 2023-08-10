rule obtain_padloc_locus_tags:
    output:
        hits    = "results/intermediate/defence_system_genes/padloc/{sample}.csv",
        gff_db  = temp('gff_{sample}.db')
    input:
        padloc  = "results/intermediate/padloc/padloc_{sample}/{sample}_padloc_filtered_renamed.csv",
        gff_gz  = "data/annotation/{sample}.gff.gz"
    run:
        # -------- Load data
        hits = pd.read_csv(input.padloc)  # padloc data
        unique_hits = hits.drop_duplicates(subset = ['target.name', 'start', 'end'])

        # create list of dictionaries for downstream lookup in gff database
        hit_dicts = [
            {
                "target.name": target_name,
                "start": start,
                "end": end,
            }
            for target_name, start, end in zip(
                unique_hits["target.name"],
                unique_hits["start"],
                unique_hits["end"]
            )
        ]

        # create database from gff file
        gff_db = gffutils.create_db(
            input.gff_gz,
            output.gff_db,
            merge_strategy="create_unique",
            force=True,
        )

        # -------- add locus tags to hit_dicts using the hit start and end
        for hit in hit_dicts:
            # the gff file may have multiple entries for a target.name (i.e. protein ID)
            # in that case, 2nd, 3rd etc. IDs have the suffix _1, _2 etc.
            # the protein IDs are prefixed by 'cds-' in the gff_db
            base_ID = "cds-" + hit["target.name"]
            IDs = [e.id for e in gff_db.all_features(featuretype="CDS") if base_ID in e.id]

            for ID in IDs:
                gff_entry = gff_db[ID]

                # verify that this gff_entry matches the hit by its
                # start and end coordinates
                # then store the locus tag in the hit dictionary
                if [gff_entry.start, gff_entry.end] == [hit["start"], hit["end"]]:
                    locus_tag = gff_entry.attributes["locus_tag"][0]
                    hit.update({"locus_tag": locus_tag})

        # -------- merge locus tags onto the padloc output table
        pd_out = hits.merge(pd.DataFrame.from_records(hit_dicts), on = ['target.name', 'start', 'end'])
        pd_out.to_csv(output.hits, index=False)
