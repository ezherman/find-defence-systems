# Create system tables (i.e. at the level of parent systems, ignoring subsystems)
rule systems_by_sample:
    output:
        csv_out = "results/intermediate/systems_by_assembly/systems_{sample}.csv"
    input:
        csv_in = "results/intermediate/subsystems_by_assembly/subsystems_{sample}.csv"
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

            # merge dsr systems, for which subsystems are not delimited with "_" in padloc
            df.system[df.system.str.startswith("dsr")] = "dsr"

            # merge dsr systems, for which subsystems are not delimited with "_" in defense_finder
            df.system[df.system.str.startswith("thoeris")] = "thoeris"

        # export
        df.to_csv(output.csv_out, index = False)
