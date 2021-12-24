# Create final systems spreadsheet
# Either defense_finder only or
# Merge systems found by defensefinder and padloc
rule final_systems_spreadsheet:
    output:
        csv = "data/analysis_data/system_tables/system_table_{assembly}.csv"
    input:
        dfinder = "data/intermediate_data/defense_finder/defense_finder_{assembly}/defense_finder_systems.tsv",
        padloc  = "data/intermediate_data/padloc/padloc_{assembly}/{assembly}_padloc.csv"
    run:

        # if padloc found hits
        if os.stat(input.padloc).st_size > 0:

            # import dataframes
            padloc = pd.read_csv(input.padloc)
            defense_finder = pd.read_table(input.dfinder)

            # wrangle dataframes
            padloc = create_system_table(padloc, "padloc")
            defense_finder = create_system_table(defense_finder, "defense_finder")

            # merge dataframes and export
            final_system_table(defense_finder, padloc).to_csv(output.csv)

        # if padloc did not find hits
        if os.stat(input.padloc).st_size == 0:

            # import and wrangle
            defense_finder = pd.read_table(input.dfinder)
            defense_finder = create_system_table(defense_finder, "defense_finder")

            # export
            final_system_table(defense_finder).to_csv(output.csv)

