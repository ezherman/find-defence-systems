# Create subsystems spreadsheet
# Either defense_finder only or
# Merge systems found by defensefinder and padloc
rule subsystems_by_sample:
    output:
        csv = "results/intermediate/subsystems_by_assembly/subsystems_{sample}.csv"
    input:
        dfinder = "results/intermediate/defense_finder/defense_finder_{sample}/defense_finder_genes_ltags_renamed.csv",
        padloc  = "results/intermediate/padloc/padloc_{sample}/{sample}_padloc_ltags_renamed.csv"
    run:

        #-------- merging and reindexing 

        # if padloc found hits
        if os.stat(input.padloc).st_size > 0:

            # import dataframes
            padloc = pd.read_csv(input.padloc)
            defense_finder = pd.read_csv(input.dfinder)

            # wrangle dataframes
            padloc = create_subsystem_table(padloc, "padloc")
            defense_finder = create_subsystem_table(defense_finder, "defense_finder")

            # merge dataframes
            df = merge_subsystem_tables(defense_finder, padloc)

        # if padloc did not find hits
        if os.stat(input.padloc).st_size == 0:

            # import and wrangle
            defense_finder = pd.read_csv(input.dfinder)
            defense_finder = create_subsystem_table(defense_finder, "defense_finder")

            # merge dataframes, still running this because the function
            # also resets the index
            df = merge_subsystem_tables(defense_finder)
        

        #-------- remove duplicates
        df['group_id'] = df.groupby('system').ngroup()
        not_singletons = df.groupby('group_id').filter(lambda x: len(x)>1)['group_id'].unique() # group_id with occurence > 1

        idx_to_exclude = []
        for one_group_id in not_singletons:
            group_data = df[df['group_id'] == one_group_id] # subset the rows belonging to this group
            group_data['protein_IDs'] = group_data['protein_IDs'].str.split(';') #turn protein_IDs into lists
            group_data['n_protein_IDs'] = group_data['protein_IDs'].str.len() #number of protein IDs per list
            
            # sort then reset index to simplify downstream looping
            group_data = group_data.sort_values(by='n_protein_IDs').reset_index(names='original_index') 

            # loop through all but the final row
            # for each row, check whether its protein_IDs are a subset of any downstream row protein IDs
            # if they are, add the row index / system ID to a list of rows to remove
            final_index = group_data.index[-1]
            for i in group_data.index[0:-1]:
                protein_IDs = group_data.loc[i]['protein_IDs']
                possible_supersets = group_data.loc[range(i+1, final_index+1)].protein_IDs.to_list()

                subset_boolean = [set(protein_IDs) <= set(e) for e in possible_supersets]
                
                if any(subset_boolean):
                    idx_to_exclude.append(group_data.loc[i]['original_index'])
            
        # exclude the entries which form a subset of another entry with the same system designation
        df = df.loc[~ df.index.isin(idx_to_exclude)].drop('group_id', axis = 1)

        #-------- remove *_other systems where their locus tags have a more informative system annotation

        # find systems that have a non-unique set of locus tags
        df['dup'] = df.duplicated(subset = 'locus_tags', keep = False)
        non_unique_locus_tags = df[df['dup'] == True]
        non_unique_locus_tags_groups = (
            non_unique_locus_tags
            .reset_index()
            .groupby('locus_tags')
            .apply(lambda x: [list(x['index']), list(x['system'])])
            .apply(pd.Series)
            .rename(columns = {0: 'indices', 1: 'systems'})
        )

        # isolate the indices associated with *_other systems that have
        # locus tags identical to at least one non-*_other system
        indices_to_remove = []
        for i in non_unique_locus_tags_groups.index:
            row = non_unique_locus_tags_groups.loc[i]
            other_systems = [s.endswith('_other') for s in row.systems]

            # at least one but less than all systems are *_other
            if sum(other_systems) > 0 and sum(other_systems) < len(other_systems):
                for j in range(0, len(row['indices'])):
                    if other_systems[j]:
                        indices_to_remove.append(row['indices'][j])

        # remove these hits
        df = df.loc[~df.index.isin(indices_to_remove)].drop('dup', axis = 1)

        #-------- export
        df.to_csv(output.csv, index = False)

