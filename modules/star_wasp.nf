process STAR_WASP {

    tag "$sample"
    conda 'bioconda::star=2.7.11b conda-forge::gzip'

    input:
    tuple val(sample), path(r2), path(r1)

    path star_index
    path whitelist
    val chemistry

    tuple val(vcf_sample), path(heterozygous_vcf)


    output:
    tuple val(sample),
          path("${sample}_WASPAligned.sortedByCoord.out.bam"),
          path("${sample}_WASPSolo.out"),
          emit: wasp_output


    script:

    def r2_arg = r2.collect { it.name }.join(',')
    def r1_arg = r1.collect { it.name }.join(',')

    """
    STAR \
        --genomeDir ${star_index} \
        --readFilesIn ${r2_arg} ${r1_arg} \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --soloType Droplet \
        --soloCBwhitelist ${whitelist} \
        --outFileNamePrefix ${sample}_WASP \
        --varVCFfile ${heterozygous_vcf} \
        --waspOutputMode SAMtag \
        --outSAMattributes vA vG NH HI AS NM nM CB UB \
        --outSAMmapqUnique 50 \
        --outSAMstrandField intronMotif \
        --outSAMattrRGline 'ID:${sample} SM:${sample} LB:${sample} PL:Illumina' \
        --runThreadN ${task.cpus} \
        --soloCBstart ${chemistry.cb_start} \
        --soloCBlen ${chemistry.cb_len} \
        --soloUMIstart ${chemistry.umi_start} \
        --soloUMIlen ${chemistry.umi_len}
    """
}