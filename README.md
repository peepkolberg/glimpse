### Recent changes
* Add UT HPC test profile and paths to data
* Automatically uses a prebuilt container from quay.io  -  no longer depends on local installations of bcftools, tabix and glimpse

### How to run test?
The test will impute chromosomes 21 and 22 of two samples (individuals). Runs ca. 40 minutes.

Submit `run_test.sh` to HPC. No need to specify parameters. Imputed genotypes will appear in directory `results.test/glimpse/imputed.ligated`.

### How to run test?
Submit `run.sh` to HPC. Imputed genotypes will appear in directory `{--outdir}/glimpse/imputed.ligated`.

Here's the descriptions and requirements for input parameters:
* **--samples** TSV file with sample IDs and BAM files. Example in [data/test_samples.tsv](https://github.com/peepkolberg/glimpse/blob/main/data/test_samples.tsv). Each sample expects an index file to be present in the same directory as the BAM file. Index file must have the same name as BAM file + ".bai".
* **--chromosomes** comma-separated list of autosomal chromosomes. Sex chromosomes will probably break the code.
* **--reference_panel** VCF(.gz) file. Expects an index file to be present in the same directory as the panel. Index file must have the same name as panel file + ".csi".
* **--reference_genome** FA file. Expects an index file to be present in the same directory as the genome. Index file must have the same name as genome file + ".fai".
* **--genetic_map** Points to a directory. In the directory, expects a separate genetic map file for each chromosome. The filenames must follow this convention: `chr{i}.gmap.gz`, where `{i}` is replaced by chromosome number.
* **--outdir** Optional parameter to specify the name of the output directory. (default: `./results`)

Some additional notes that might be necessary to follow:
* Sort input BAM files

To-Do:
* Tweak each process' memory, CPU, time requirements
* Refactor the workflow so that process sample_GLs runs in parallel with process chunk_chr. Might save a few minutes.

Workflow follows the steps from this tutorial: https://odelaneau.github.io/GLIMPSE/docs/tutorials/getting_started/#imputation-using-1000gp-reference-panel

#### Good luck!
