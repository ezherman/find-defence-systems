## TO BE PASTED INTO THE COMMAND LINE

conda create -n snakemake -c conda-forge -c bioconda -c bioinf-mcb pandas pyjanitor biopython snakemake mamba pip sra-downloader entrez-direct -y

conda activate snakemake
pip install wget

conda env export  >  ../envs/snakemake.yaml
