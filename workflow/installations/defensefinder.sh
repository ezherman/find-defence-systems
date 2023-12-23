## TO BE PASTED INTO THE COMMAND LINE
## FROM THE ROOT OF THE REPOSITORY

mamba create -c bioconda -c conda-forge -n defensefinder pip hmmer -y

mamba activate defensefinder

pip install mdmparis-defense-finder

mamba env export > workflow/envs/defensefinder.yaml
