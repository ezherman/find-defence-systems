## TO BE PASTED INTO THE COMMAND LINE
## FROM THE ROOT OF THE REPOSITORY

mamba create -c conda-forge -c bioconda -n r_env r-base r-readr r-dplyr r-circlize bioconductor-complexheatmap -y

mamba env export -n r_env > workflow/envs/r.yaml
