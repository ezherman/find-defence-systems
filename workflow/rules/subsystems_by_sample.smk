# Create subsystems spreadsheet
# Either defense_finder only or
# Merge systems found by defensefinder and padloc
rule subsystems_by_sample:
    output:
        csv = "results/intermediate/subsystems_by_assembly/subsystems_{sample}.csv"
    input:
        dfinder = "results/intermediate/defense_finder/defense_finder_{sample}/defense_finder_genes_renamed.csv",
        padloc  = "results/intermediate/padloc/padloc_{sample}/{sample}_padloc.csv"
    run:

        # if padloc found hits
        if os.stat(input.padloc).st_size > 0:

            # import dataframes
            padloc = pd.read_csv(input.padloc)
            defense_finder = pd.read_csv(input.dfinder)

            # wrangle dataframes
            padloc = create_subsystem_table(padloc, "padloc")
            defense_finder = create_subsystem_table(defense_finder, "defense_finder")

            # merge dataframes and export
            merge_subsystem_tables(defense_finder, padloc).to_csv(output.csv)

        # if padloc did not find hits
        if os.stat(input.padloc).st_size == 0:

            # import and wrangle
            defense_finder = pd.read_csv(input.dfinder)
            defense_finder = create_subsystem_table(defense_finder, "defense_finder")

            # export
            merge_subsystem_tables(defense_finder).to_csv(output.csv)

