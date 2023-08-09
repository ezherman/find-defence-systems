# Rename systems identified by defensefinder 
# so that they are comparable to the renamed padloc output
rule rename_padloc_systems:
    output:
        padloc  = "results/intermediate/padloc/padloc_{sample}/{sample}_padloc_renamed.csv"
    input:
        padloc  = "results/intermediate/padloc/padloc_{sample}/{sample}_padloc.csv"
    run:
        #-------- import data
        padloc = pd.read_table(input.padloc)

        #-------- rename the systems
        # generalise_system_names is defined in workflow/rules/common.smk
        padloc['system'] = padloc.apply(generalise_system_names, software="padloc", axis=1)

        #-------- export the data
        padloc.to_csv(output.padloc, index = False)