# Create system or subsystem tables across all assemblies
rule matrix_all_samples:
    output:
        mat = "results/main/matrix_{level}.csv"
    input:
        mats = expand("results/intermediate/{{level}}_by_assembly/{{level}}_{sample}.csv", sample = SAMPLE_NAMES)
    run:
        matrix_all_samples(files = input.mats, binary = True).to_csv(output.mat)
