## TO BE PASTED INTO THE COMMAND LINE
## FROM THE ROOT OF THE REPOSITORY

mamba create -c conda-forge -c bioconda -c bioinf-mcb -n snakemake pandas pyjanitor biopython snakemake mamba pip sra-downloader entrez-direct -y

mamba activate snakemake
pip install wget

mamba env export  >  workflow/envs/snakemake.yaml
