#!/bin/bash

#SBATCH --job-name=glimpse
#SBATCH --partition=amd
#SBATCH --time=1-00:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

module load nextflow/22.04.3
module load any/singularity/3.11.1

# Maybe sample bam files need to be sorted. An example is given in data/samples.tsv. Each sample expects an index file to be present in the same directory as the bam file. Index file must have the same name as bam file + ".bai".
# chromosomes is a comma-separated list of autosomal chromosomes. Sex chromosomes will probably break the code.
# Reference panel and genome expect index files to be present in the same directories as the references. Index files must have the same name as reference file + ".csi" for reference panel and + ".fai" for reference genome.
# genetic_map points to a directory. In the directory, expects a separate file for each chromosome. The filenames must follow this convention: chr{i}.gmap.gz, where {i} is replaced by chromosome number.
nextflow -log logs/.nextflow.log run main.nf -profile tartu_hpc \
    --samples /path/to/samples.tsv \
    --chromosomes 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 \
    --reference_panel /gpfs/space/projects/genomic_references/imputation/glimpse_reference_070223/reference_panel/1kGP_high_coverage_Illumina.all.filtered.SNV_INDEL_SV_phased_panel.no_chr.independent_inds.vcf.gz \
    --reference_genome /gpfs/space/projects/genomic_references/imputation/glimpse_reference_070223/reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --genetic_map /gpfs/space/projects/genomic_references/imputation/glimpse_reference_070223/genetic_map.b38 \
    --outdir ./results \
#    -resume
