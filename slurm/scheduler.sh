#!/bin/bash
#SBATCH --job-name=snakemake	         	# Job name
#SBATCH --mail-type=END,FAIL             	# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=elh605@york.ac.uk    	# Where to send mail  
#SBATCH --ntasks=1                       	# Run on a single CPU
#SBATCH --mem=4gb                        	# Job memory request
#SBATCH -c 1				 	# CPUs
#SBATCH --time=03:00:00                  	# Time limit hrs:min:sec
#SBATCH --output=snakemake_%j.log    		# Standard output and error log
#SBATCH --account=biol-evophage-2020     	# Project account

cd ../

mkdir -p logs_slurm
snakemake -F --use-conda --profile slurm/slurm_profile 
