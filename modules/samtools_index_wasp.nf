process SAMTOOLS_INDEX_WASP {

    tag "$sample"

    input:

    tuple val(sample), path(bam), path(solo)


    output:

    tuple val(sample),
          path(bam),
          path("${bam}.bai"),
          path(solo),
          emit: indexed_wasp_bam


    script:

    """
    samtools index \
        -@ ${task.cpus} \
        ${bam}
    """
}