# Create final systems spreadsheet
# Either defense_finder only or
# Merge systems found by defensefinder and padloc
rule systems_by_assembly:
    output:
        csv = "results/intermediate/systems_by_assembly/systems_{assembly}.csv"
    input:
        dfinder = "results/intermediate/defense_finder/defense_finder_{assembly}/defense_finder_systems.tsv",
        padloc  = "results/intermediate/padloc/padloc_{assembly}/padloc_{assembly}.csv"
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

