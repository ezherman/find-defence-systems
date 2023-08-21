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
        

        # -------- remove drt_class_* systems
        # these are broadly defined system classes specific to padloc
        # also tend to overlap with other systems
        df = df[~df["system"].str.contains("drt_class")]

        # -------- remove duplicates
        df["group_id"] = df.groupby("system").ngroup()
        not_singletons = (
            df.groupby("group_id").filter(lambda x: len(x) > 1)["group_id"].unique()
        )  # group_id with occurence > 1

        idx_to_exclude = []
        for one_group_id in not_singletons:
            group_data = df[
                df["group_id"] == one_group_id
            ]  # subset the rows belonging to this group
            group_data["protein_IDs"] = group_data["protein_IDs"].str.split(
                ";"
            )  # turn protein_IDs into lists
            group_data["n_protein_IDs"] = group_data[
                "protein_IDs"
            ].str.len()  # number of protein IDs per list

            # sort then reset index to simplify downstream looping
            group_data = group_data.sort_values(by="n_protein_IDs").reset_index(
                names="original_index"
            )

            # loop through all but the final row
            # for each row, check whether its protein_IDs are a subset of any downstream row protein IDs
            # if they are, add the row index / system ID to a list of rows to remove
            final_index = group_data.index[-1]
            for i in group_data.index[0:-1]:
                protein_IDs = group_data.loc[i]["protein_IDs"]
                possible_supersets = group_data.loc[
                    range(i + 1, final_index + 1)
                ].protein_IDs.to_list()

                subset_boolean = [set(protein_IDs) <= set(e) for e in possible_supersets]

                if any(subset_boolean):
                    idx_to_exclude.append(group_data.loc[i]["original_index"])

        # exclude the entries which form a subset of another entry with the same system designation
        df = df.loc[~df.index.isin(idx_to_exclude)].drop("group_id", axis=1)

        # -------- remove *_other systems where their locus tags have a more informative system annotation

        # turn locus tag strings into sets
        df["locus_tags_set"] = df["locus_tags"].str.split(";").apply(set)

        # separate data into _other systems and non _other
        df_other = df[["_other" in s for s in df["system"]]]
        df_non_other = df[[i not in df_other.index for i in df.index]]

        sys_indices_to_remove = []

        for i in df_other.index:
            one_other_system = df_other.loc[i]["locus_tags_set"]

            # check whether one_other_system is a subset of any non _other system
            subsets = [one_other_system <= ltags for ltags in df_non_other["locus_tags_set"]]

            if any(subsets):
                sys_indices_to_remove.append(i)

        # remove the _other systems which are subsets of non _other systems
        df = df.loc[~df.index.isin(sys_indices_to_remove)].drop("locus_tags_set", axis=1)

        # -------- remove cases where a family annotation falls within a type annotation
        # e.g. hachiman vs hachiman_i, only if the hachiman locus_tags fall within those of hachiman_i

        # find hits with systems not defined at type level, for which a type level annotation does exist
        # some family-level systems are defined with '-fam', e.g. lamassu-fam
        systems_to_check = {
            s
            for s in df.system
            if any(df["system"].str.contains("|".join([s + "_", s.split("-fam")[0] + "_"])))
        }

        for system in systems_to_check:
            if system.endswith("-fam"):
                system_clean = system.split("-fam")[0]
            else:
                system_clean = system

            # find locus tags associated with the systems defined at type level
            # system + '_' indicates there's a type definition
            systems_with_type = df["system"][
                df.system.str.contains(system_clean + "_")
            ].to_list()
            locus_tag_sets = df["locus_tags"].str.split(";").apply(set)
            locus_tag_sets_with_type = locus_tag_sets[
                df["system"].isin(systems_with_type)
            ].to_list()

            for index in df[df["system"] == system].index:
                locus_tag_set = set(df.loc[index]["locus_tags"].split(";"))

                if any([locus_tag_set <= ltags for ltags in locus_tag_sets_with_type]):
                    df = df.drop(index)

        #-------- export
        df.to_csv(output.csv, index = False)

