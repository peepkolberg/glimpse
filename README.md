To run the workflow, submit `run.sh` to HPC. Imputation results will appear in directory `results/glimpse`.

Here's the descriptions and requirements for input parameters:
* **--samples** TSV file with sample IDs and BAM files. Example in [data/samples.tsv](https://github.com/peepkolberg/glimpse/blob/main/data/samples.tsv). Each sample expects an index file to be present in the same directory as the BAM file. Index file must have the same name as BAM file + ".bai". This is hard-coded in main.nf workflow{} section.
* **--chromosomes** comma-separated list of autosomal chromosomes. Sex chromosomes will probably break the code.
* **--reference_panel** VCF(.gz) file. Expects an index file to be present in the same directory as the panel. Index file must have the same name as panel file + ".csi". This is hard-coded in main.nf workflow{} section.
* **--reference_genome** FA file. Expects an index file to be present in the same directory as the genome. Index file must have the same name as genome file + ".fai". This is hard-coded in main.nf workflow{} section.
* **--genetic_map** Points to a directory. In the directory, expects a separate genetic map file for each chromosome. The filenames must follow this convention: `chr{i}.gmap.gz`, where `{i}` is replaced by chromosome number.

Some additional notes that might be necessary to follow:
* Sort input BAM files

To-Do:
* Tweak each process' memory, CPU, time requirements
* Refactor the workflow so that process sample_GLs runs in parallel with process chunk_chr. Might save a few minutes.

Workflow follows the steps from this tutorial: https://odelaneau.github.io/GLIMPSE/docs/tutorials/getting_started/#imputation-using-1000gp-reference-panel

#### Good luck!
