rule heatmap:
    output: "results/main/heatmap_systems.png"
    input: "results/main/matrix_systems.csv"
    params:
        samples = config["samples"]
    conda: "../envs/r.yaml"
    script:
        "../scripts/heatmap.R"
