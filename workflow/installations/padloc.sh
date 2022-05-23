## TO BE PASTED INTO THE COMMAND LINE

conda create -n padloc -c conda-forge -c bioconda -c padlocbio padloc

conda activate padloc

conda env export > ../envs/padloc.yaml
