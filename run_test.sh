#!/bin/bash

#SBATCH --job-name=glimpse
#SBATCH --partition=amd
#SBATCH --time=1-00:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

module load nextflow/22.04.3
module load any/singularity/3.11.1

nextflow -log logs/.nextflow.log run main.nf -profile tartu_hpc_test \
    --out_dir results.test \
    -resume
