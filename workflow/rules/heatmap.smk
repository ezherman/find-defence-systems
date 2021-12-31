rule heatmap:
    output: "results/main/heatmap.png"
    input: "results/main/systems_summarised.csv"
    conda: "../envs/r.yaml"
    script:
        "../scripts/heatmap.R"
