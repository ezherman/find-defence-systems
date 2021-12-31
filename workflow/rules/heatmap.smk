rule heatmap:
    output: "results/main/heatmap.png"
    input: "results/main/systems_summarised.csv"
    params:
        samples = config["samples"]
    conda: "../envs/r.yaml"
    script:
        "../scripts/heatmap.R"
