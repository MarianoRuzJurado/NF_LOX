process SAMTOOLS_INDEX {

    tag "$sample"

    input:

    /*
     * BAM and Solo output produced by STAR.
     */
    tuple val(sample), path(bam), path(solo)


    output:

    /*
     * Pass the BAM, its index, and the Solo directory onward.
     */
    tuple val(sample),
          path(bam),
          path("${bam}.bai"),
          path(solo),
          emit: indexed_bam


    script:

    """
    samtools index \
        -@ ${task.cpus} \
        ${bam}
    """
}