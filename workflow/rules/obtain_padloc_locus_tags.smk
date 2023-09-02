rule obtain_padloc_locus_tags:
    output:
        hits    = "results/intermediate/padloc/padloc_{sample}/{sample}_padloc_ltags.csv",
        gff_db  = temp('gff_{sample}.db')
    input:
        padloc  = "results/intermediate/padloc/padloc_{sample}/{sample}_padloc.csv",
        gff_gz  = "data/annotation/{sample}.gff.gz"
    run:
        #-------- if padloc did not find hits, create empty output files
        if os.stat(input.padloc).st_size == 0:
            shell('cp {input.padloc} {output.hits}')
            shell('touch {output.gff_db}')

        else:
            # -------- Load data 
            hits = pd.read_csv(input.padloc)  # padloc data

            # first check if annotation was performed with Bakta
            # if so, the target.name column is already the locus tag
            with gzip.open(input.gff_gz, "r") as f:
                lines = f.readlines()
                bakta_true = "Bakta" in str(lines[2])

            if bakta_true:
                hits["locus_tag"] = hits["target.name"]
                hits.to_csv(output.hits, index=False)

            # if not bakta (i.e. refseq), link target name to locus tag using the gff file
            else:
                unique_hits = hits.drop_duplicates(subset=["target.name", "start", "end"])

                # create list of dictionaries for downstream lookup in gff database
                hit_dicts = [
                    {
                        "target.name": target_name,
                        "start": start,
                        "end": end,
                    }
                    for target_name, start, end in zip(
                        unique_hits["target.name"], unique_hits["start"], unique_hits["end"]
                    )
                ]

                # create database from gff file
                gff_db = gffutils.create_db(
                    "pao1.gff.gz",
                    "pao1.db",
                    merge_strategy="create_unique",
                    force=True,
                )

                # -------- add locus tags to hit_dicts using the hit start and end
                for hit in hit_dicts:
                    # the gff file may have multiple entries for a target.name (i.e. protein ID)
                    # in that case, 2nd, 3rd etc. IDs have the suffix _1, _2 etc.
                    # the protein IDs are prefixed by 'cds-' in the gff_db
                    base_ID = "cds-" + hit["target.name"]
                    IDs = [
                        e.id
                        for e in gff_db.all_features(featuretype="CDS")
                        if any(n in e.id for n in [base_ID, hit["target.name"]])
                    ]

                    for ID in IDs:
                        gff_entry = gff_db[ID]

                        # verify that this gff_entry matches the hit by its
                        # start and end coordinates
                        # then store the locus tag in the hit dictionary
                        if [gff_entry.start, gff_entry.end] == [hit["start"], hit["end"]]:
                            locus_tag = gff_entry.attributes["locus_tag"][0]
                            hit.update({"locus_tag": locus_tag})

                # -------- merge locus tags onto the padloc output table
                pd_out = hits.merge(
                    pd.DataFrame.from_records(hit_dicts), on=["target.name", "start", "end"]
                )
                pd_out.to_csv(output.hits, index=False)
