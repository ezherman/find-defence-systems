## TO BE PASTED INTO THE COMMAND LINE
## FROM THE ROOT OF THE REPOSITORY

mamba create -c conda-forge -c bioconda -n snakemake pandas snakemake mamba pip biopython openpyxl -y

mamba activate snakemake
pip install wget

mamba env export  >  workflow/envs/snakemake.yaml
