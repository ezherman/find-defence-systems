## TO BE PASTED INTO THE COMMAND LINE

conda create -c bioconda -c conda-forge -n defensefinder pip=21.0 hmmer

conda activate defensefinder

pip install mdmparis-defense-finder

conda env export > ../envs/defensefinder.yaml
