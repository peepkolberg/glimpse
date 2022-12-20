#!/bin/bash

#SBATCH --job-name=glimpse
#SBATCH --partition=amd
#SBATCH --time=1-00:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G



module load any/jdk/1.8.0_265
module load any/singularity/3.5.3
module load nextflow
module load any/bcftools
module load any/glimpse/1.1.1

# Workflow follows the steps from this tutorial: https://odelaneau.github.io/GLIMPSE/docs/tutorials/getting_started/#imputation-using-1000gp-reference-panel
# maybe sample bam files need to be sorted. An example is given in data/samples.tsv. Each sample expects an index file to be present in the same directory as the bam file. Index file must have the same name as bam file + ".bai". This is hard-coded in main.nf workflow{} section.
# chromosomes is a comma-separated list of autosomal chromosomes. Sex chromosomes will probably break the code.
# Reference panel and genome expect index files to be present in the same directories as the references. Index files must have the same name as reference file + ".csi" for reference panel and + ".fai" for reference genome. These are hard-coded in main.nf workflow{} section.
# genetic map points to a directory. In the directory, expects a separate file for each chromosome. The filenames must follow this convention: chr{i}.gmap.gz, where {i} is replaced by chromosome number.
# To tweak: each process' memory, CPU, time limits
# TODO: For a small (maybe a few minutes) time save, refactor the workflow so that process sample_GLs runs in parallel with process chunk_chr.
nextflow run main.nf -profile tartu_hpc \
    --samples data/samples.tsv \
    --chromosomes 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 \
    --reference_panel /path/to/reference_panel.vcf.gz \
    --reference_genome /path/to/reference_genome.fa \
    --genetic_map /path/to/genetic_map \
#    -resume
