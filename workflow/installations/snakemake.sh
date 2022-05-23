## TO BE PASTED INTO THE COMMAND LINE

conda create -n snakemake -c conda-forge -c bioconda pandas pyjanitor biopython snakemake mamba pip -y

conda activate snakemake
pip install wget

conda env export  >  ../envs/snakemake.yaml
