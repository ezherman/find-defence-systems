conda create -c conda-forge -c bioconda -n r_env r-base r-readr r-dplyr r-circlize bioconductor-complexheatmap 

conda activate r_env

conda env export > ../envs/y.yaml
