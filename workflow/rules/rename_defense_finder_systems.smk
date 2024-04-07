# Rename systems identified by defensefinder 
# so that they are comparable to the renamed padloc output
rule rename_defense_finder_systems:
    output:
        dfinder = "results/intermediate/defense_finder/defense_finder_{sample}/defense_finder_genes_ltags_renamed.csv"
    input:
        dfinder = "results/intermediate/defense_finder/defense_finder_{sample}/defense_finder_genes_ltags.csv",
    run:
        #-------- if defense finder did not find hits, create empty output files
        dfinder = pd.read_table(input.dfinder)
        if len(hits) == 0:
            dfinder['system'] = []
            dfinder['sys_id'] = []
            dfinder.to_csv(output.hits, index = False)

        else:   
            #-------- wrangle data, then rename the systems
            # system name is sys_id with the UserReplicon_ prefix and numerical suffix (e.g. _1) removed
            # the suffix increases (to e.g. _2) when multiple instances of the same system are found
            # within an isolate
            dfinder['system'] = ['_'.join(d.split('_')[1:-1]) for d in dfinder['sys_id']]
            dfinder['sys_id'] = dfinder.groupby(['sys_id']).ngroup() # simple numerical ID
    
            # rename the systems
            # generalise_system_names is defined in workflow/rules/common.smk
            dfinder['system'] = dfinder.apply(generalise_system_names, software="defense_finder", axis=1)
    
            #-------- export the data
            dfinder.to_csv(output.dfinder, index = False)
