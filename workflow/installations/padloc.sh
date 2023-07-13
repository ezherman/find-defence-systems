## TO BE PASTED INTO THE COMMAND LINE
## FROM THE ROOT OF THE REPOSITORY

mamba create -n padloc -c conda-forge -c bioconda -c padlocbio padloc -y

mamba env export -n padloc > workflow/envs/padloc.yaml
