# Create system tables at the level of parent systems (i.e. ignoring subsystems)
rule final_systems_spreadsheet_summary:
    output:
        csv_out = "data/analysis_data/system_tables_summarised/system_table_summarised_{assembly}.csv"
    input:
        csv_in = "data/analysis_data/system_tables/system_table_{assembly}.csv"
    run:

        df = pd.read_csv(input.csv_in)

        # if no systems were found, export the empty table
        if df.empty:
            df.to_csv(output.csv_out)

        else:

            # system-level results
            df.system = df.system.str.split("_", expand = True)[0]

            # remove capital letters for cross-compatibility between programs
            # e.g. Gabija in defense_finder and gabija in padloc
            df.system = df.system.str.lower()

            # merge abi systems, for which subsystems are not delimited with "_" in padloc
            # and some abi systems do not contain abi in their name,
            # see https://github.com/padlocbio/padloc-db/blob/master/sys_meta.txt
            abortive_systems = ("abi", "rexab", "lit", "bsta", "old", "tin", "prrc", "pifa")
            df.system[df.system.str.startswith(abortive_systems)] = "abi"

            # merge dsr systems, for which subsystems are not delimited with "_" in padloc
            df.system[df.system.str.startswith("dsr")] = "dsr"

            # merge dsr systems, for which subsystems are not delimited with "_" in defense_finder
            df.system[df.system.str.startswith("thoeris")] = "thoeris"

        # export
        df.to_csv(output.csv_out, index = False)
