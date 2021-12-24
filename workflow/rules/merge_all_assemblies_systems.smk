# Create system and subsystem tables across all assemblies
rule merge_all_assemblies_systems:
    output:
        sys = "data/analysis_data/system_tables_summarised/system_table_summarised_all_assemblies.csv",
        sub = "data/analysis_data/system_tables/system_table_all_assemblies.csv"
    input:
        sys_dir         = "data/analysis_data/system_tables_summarised",
        sub_dir         = "data/analysis_data/system_tables",
        summary_output  = expand("data/analysis_data/system_tables_summarised/system_table_summarised_{assembly}.csv", assembly = AS>    
    run:
        merge_assemblies_systems(input.sys_dir).to_csv(output.sys)
        merge_assemblies_systems(input.sub_dir).to_csv(output.sub)
