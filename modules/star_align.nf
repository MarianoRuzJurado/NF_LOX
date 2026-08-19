process STAR_ALIGN {

    tag "$sample"

    input:
    tuple val(sample), path(r2), path(r1)

    path star_index
    path whitelist
    val chemistry

    output:
    tuple val(sample),
          path("${sample}Aligned.sortedByCoord.out.bam"),
          path("${sample}Solo.out"),
          emit: star_output

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
        --outFileNamePrefix ${sample} \
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