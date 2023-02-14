// Follows the tutorial at https://odelaneau.github.io/GLIMPSE/docs/tutorials/getting_started/#imputation-using-1000gp-reference-panel

nextflow.enable.dsl=2



samples_ch = Channel.fromPath(params.samples)
    .splitCsv(header: true, sep: "\t", strip: true)
    .map( row -> [row.sample, row.bam] )

chrs_ch = Channel.fromList(params.chromosomes.split(",") as List)

ref_pan_ch = Channel.fromPath(params.reference_panel)

ref_gen_ch = Channel.fromPath(params.reference_genome)

genetic_map_ch = Channel.fromPath(params.genetic_map, type: "dir")



process split_ref_pan {
    input:
        tuple val(chr), path(ref_pan), path(ref_pan_idx)

    output:
        tuple val(chr), path(ref_pan_chr), path(ref_pan_chr_idx)

    script:
        ref_pan_chr = "ref_pan.chr${chr}.vcf.gz"
        ref_pan_chr_idx = "${ref_pan_chr}.csi"

        """
        bcftools view -r ${chr} -m2 -M2 -v snps --threads 4 ${ref_pan} -Oz -o ${ref_pan_chr}
        bcftools index -f --threads 4 ${ref_pan_chr}
        """
}

process extract_ref_pan_sites {
    input:
        tuple val(chr), path(ref_pan_chr), path(ref_pan_chr_idx)

    output:
        tuple val(chr), path(ref_pan_chr), path(ref_pan_chr_idx), path(sites_vcf), path(sites_vcf_idx), path(sites_tsv), path(sites_tsv_idx)

    script:
        sites_vcf = "ref_pan.chr${chr}.sites.vcf.gz"
        sites_tsv = "ref_pan.chr${chr}.sites.tsv.gz"
        sites_vcf_idx = "${sites_vcf}.csi"
        sites_tsv_idx = "${sites_tsv}.tbi"

        """
        bcftools view --drop-genotypes --threads 4 ${ref_pan_chr} -Oz -o ${sites_vcf}
        bcftools index --threads 4 -f ${sites_vcf}
        bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' ${sites_vcf} | bgzip --threads 4 -c > ${sites_tsv}
        tabix -s1 -b2 -e2 ${sites_tsv}
        """
}

process chunk_chr {
    input:
        tuple val(chr), path(ref_pan_chr), path(ref_pan_chr_idx), path(sites_vcf), path(sites_vcf_idx), path(sites_tsv), path(sites_tsv_idx)

    output:
        tuple val(chr), path(ref_pan_chr), path(ref_pan_chr_idx), path(sites_vcf), path(sites_vcf_idx), path(sites_tsv), path(sites_tsv_idx), path(chunks)

    script:
        window_size = 2000000
        buffer_size = 200000
        chunks = "ref_pan.chr${chr}.chunks.txt"

        """
        GLIMPSE_chunk --input ${sites_vcf} --region ${chr} --window-size ${window_size} --buffer-size ${buffer_size} --thread 4 --output ${chunks}
        """
}

// TODO: Run in parallel with the previous process.
// TODO: To impute X-chromosome, give --ploidy-file to bcftools call. Use val(chr) to determine ploidy. 
process sample_GLs {
    input:
        tuple val(sample), path(bam), path(bam_idx), path(ref_gen), path(ref_gen_idx), val(chr), path(ref_pan_chr), path(ref_pan_chr_idx), path(sites_vcf), path(sites_vcf_idx), path(sites_tsv), path(sites_tsv_idx), path(chunks)

    output:
        tuple val(sample), val(chr), path(ref_pan_chr), path(ref_pan_chr_idx), path(chunks), path(sample_GLs), path(sample_GLs_idx)

    script:
        sample_GLs = "${sample}.chr${chr}.GLs.vcf.gz"
        sample_GLs_idx = "${sample_GLs}.csi"
        sample_name = "sample_name.txt"
        sample_GLs_temp = "${sample}.chr${chr}.GLs.new_samplename.vcf.gz"

        """
        bcftools mpileup -f ${ref_gen} --skip-indels --redo-BAQ --annotate 'FORMAT/DP' -T ${sites_vcf} --regions ${chr} ${bam} -Ou | bcftools call --threads 4 -Aim -C alleles -T ${sites_tsv} -Oz -o ${sample_GLs}
        echo $sample > $sample_name
        bcftools reheader --samples $sample_name --threads 4 $sample_GLs -o $sample_GLs_temp
        mv $sample_GLs_temp $sample_GLs
        bcftools index -f $sample_GLs
        """
}

process impute {
    input:
        tuple val(sample), val(chr), path(ref_pan_chr), path(ref_pan_chr_idx), path(chunks), path(sample_GLs), path(sample_GLs_idx), path(map)

    output:
        tuple val(sample), val(chr), path(ligated), path(ligated_idx)

    script:
        imputed_filenames = "imputed_filenames.txt"
        ligated = "${sample}.chr${chr}.imputed.ligated.bcf"
        ligated_idx = "${ligated}.csi"

        """
        while ifs="" read -r line || [ -n "\$line" ]; do
            printf -v id "%02d" \$(echo \$line | cut -d" " -f1)
            irg=\$(echo \$line | cut -d" " -f3)
            org=\$(echo \$line | cut -d" " -f4)
            imputed=${sample}.chr${chr}.imputed.\${id}.bcf
            GLIMPSE_phase --input ${sample_GLs} --reference ${ref_pan_chr} --map ${map} --input-region \$irg --output-region \$org --output \$imputed
            bcftools index --threads 4 -f \$imputed
        done < ${chunks}

        ls ${sample}.chr${chr}.imputed.*.bcf > ${imputed_filenames}
        GLIMPSE_ligate --thread 4 --input ${imputed_filenames} --output ${ligated}
        bcftools index --threads 4 -f ${ligated}
        """
}

process sample {
    publishDir "${projectDir}/results/glimpse/phased/${sample}/chr${chr}", mode: "move", pattern: "${phased}*"
    
    input:
        tuple val(sample), val(chr), path(ligated), path(ligated_idx)
    
    output:
        tuple val(sample), path(ligated), path(ligated_idx), path(phased), path(phased_idx)
    
    script:
        phased = "${sample}.chr${chr}.phased.bcf"
        phased_idx = "${phased}.csi"

        """
        GLIMPSE_sample --thread 4 --input ${ligated} --solve --output ${phased}
        bcftools index --threads 4 -f ${phased}
        """    
}

// Everything from this point on is my original creation.
process concat_chrs {
    input:
        tuple val(sample), path(chromosomes_bcfs), path(chromosomes_idxs)

    output:
        tuple path(sorted), path(sorted_idx)

    script:
        concat_filenames = "concat_filenames.txt"
        sorted = "${sample}.all_chrs.sorted.bcf"
        sorted_idx = "${sorted}.csi"

        """
        ls ${sample}.chr*.imputed.ligated.bcf > $concat_filenames
        bcftools concat --threads 4 --file-list $concat_filenames -Ob | bcftools sort -Ob -o $sorted
        bcftools index --threads 4 $sorted
        """
}

process merge_inds {
    publishDir "${projectDir}/results/glimpse/imputed.ligated", mode: "move"

    input:
        path(samples_files)

    output:
        tuple path(merged), path(merged_idx)

    script:
        merge_filenames = "merge_filenames.txt"
        merged = "merged.imputed.ligated.vcf.gz"
        merged_idx = "${merged}.csi"

        """
        ls *.all_chrs.sorted.bcf > $merge_filenames
        bcftools merge -m id --file-list $merge_filenames --threads 4 -Oz -o $merged
        bcftools index --threads 4 $merged
        """
}



workflow {
    //println "chr, ref_pan"
    ref_pan_chr_ch = chrs_ch
        .combine(ref_pan_ch)
        .map( it -> [it[0], it[1], it[1]+".csi"] )
    //ref_pan_chr_ch.view()

    //println "chr, ref_pan, sites, chunks"
    split_ref_pan_ch = split_ref_pan(ref_pan_chr_ch) | extract_ref_pan_sites | chunk_chr
    //split_ref_pan_ch.view()


    //println "sample, bam, ref_gen, chr, ref_pan, sites, chunks"
    job_params_ch = samples_ch
        .map( it -> [it[0], it[1], it[1]+".bai"] )
        .combine(ref_gen_ch)
        .map( it -> [it[0], it[1], it[2], it[3], it[3]+".fai"] )
        .combine(split_ref_pan_ch)
    //job_params_ch.view()

    //println "sample, chr, ref_pan, chunks, sample_GLs"
    job_params_ch = sample_GLs(job_params_ch)
    //job_params_ch.view()
        

    //println "sample, chr, ref_pan, chunks, sample_GLs, map"
    job_params_ch = job_params_ch
        .combine(genetic_map_ch)
        .map( it -> [it[0], it[1], it[2], it[3], it[4], it[5], it[6], it[7]+"/chr"+it[1]+".gmap.gz"] )
    //job_params_ch.view()

    //println "sample, ligated, phased"
    imputed_chromosomes_ch = impute(job_params_ch) | sample 
    //imputed_chromosomes_ch.view()

    //println "sample, [ligated]"
    imputed_chromosomes_grouped_ch = imputed_chromosomes_ch
        .map( it -> it[0..2] )  // phased files were published in the previous process. Continue only with ligated.
        .groupTuple()
    //imputed_chromosomes_grouped_ch.view()

    //println "merged"
    concat_chrs(imputed_chromosomes_grouped_ch) | \
    collect | \
    merge_inds //| view()
}
