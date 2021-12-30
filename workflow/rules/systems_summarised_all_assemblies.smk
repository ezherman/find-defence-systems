# Create system and subsystem tables across all assemblies
rule systems_summarised_all_assemblies:
    output:
        sys = "results/main/systems_summarised.csv",
        sub = "results/main/systems.csv"
    input:
        sys_dir         = "results/intermediate/systems_summarised_by_assembly",
        sub_dir         = "results/intermediate/systems_by_assembly",
        summary_output  = expand("results/intermediate/systems_summarised_by_assembly/systems_summarised_{assembly}.csv", assembly = ASSEMBLIES)
    run:
        summarise_all_assemblies(input.sys_dir).to_csv(output.sys)
        summarise_all_assemblies(input.sub_dir).to_csv(output.sub)
