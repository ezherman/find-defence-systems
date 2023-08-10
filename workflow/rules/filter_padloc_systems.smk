# sometimes padloc calls a system as *_other (e.g. DMS_Other)
# when the same genes have also been called as a known system (e.g. RM_Type_I)
# remove the *_other calls in those cases
rule filter_padloc_systems:
    output:
        padloc  = "results/intermediate/padloc/padloc_{sample}/{sample}_padloc_ltags_filtered.csv"
    input:
        padloc  = "results/intermediate/padloc/padloc_{sample}/{sample}_padloc_ltags.csv"
    run:
        #-------- import data
        padloc = pd.read_csv(input.padloc)

        #-------- filter data
        # find systems that have the same start and end positions for their genes
        all_positions = (
            padloc
            .groupby(['system.number', 'system'])
            .apply(lambda x: [list(x['start']), list(x['end'])])
            .apply(pd.Series)
            .rename(columns = {0:'start', 1:'end'})
        )

        # find rows with the same start and end positions
        # astype(str) is needed because the start and end columns contain lists
        all_redundant_positions = all_positions[all_positions.astype(str).duplicated(keep = False)] # filter rows which have non-unique start and end values

        # for each non-unique start and end positions, find the system numbers and systems
        all_redundant_systems = (
            all_redundant_positions
            .astype(str)
            .reset_index()
            .groupby(['start', 'end'])
            .apply(lambda x: [list(x['system.number']), list(x['system'])])
            .apply(pd.Series)
            .rename(columns = {0: 'system.number', 1: 'system'})
        )

        # isolate the system numbers associated with *_other systems that share
        # start and end positions with at least one non-*_other system
        system_numbers_to_remove = []
        for i in all_redundant_systems.index:
            row = all_redundant_systems.loc[i]
            other_systems = [s.endswith('_other') for s in row.system]

            # at least one but less than all systems are *_other
            if sum(other_systems) > 0 and sum(other_systems) < len(other_systems):
                for j in range(0, len(row['system.number'])):
                    if other_systems[j]:
                        system_numbers_to_remove.append(row['system.number'][j])

        # remove these hits
        padloc = padloc.loc[~padloc['system.number'].isin(system_numbers_to_remove)]        

        #-------- export the data
        padloc.to_csv(output.padloc, index = False)