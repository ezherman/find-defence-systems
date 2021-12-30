# Create system and subsystem tables across all assemblies
rule systems_all_assemblies:
    output:
        sys = "results/main/systems_summarised.csv",
        sub = "results/main/systems.csv"
    input:
        summary_output  = expand("results/intermediate/systems_by_assembly_summarised/systems_summarised_{assembly}.csv", assembly = ASSEMBLIES),
        full_output     = expand("results/intermediate/systems_by_assembly/systems_{assembly}.csv", assembly = ASSEMBLIES)
    run:
        systems_all_assemblies(os.path.dirname(input.summary_output[0])).to_csv(output.sys)
        systems_all_assemblies(os.path.dirname(input.full_output[0])).to_csv(output.sub)
