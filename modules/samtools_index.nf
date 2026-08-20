process SAMTOOLS_INDEX {

    tag "$sample"
    conda 'bioconda::samtools=1.22.1'

    input:

    /*
     * BAM and Solo output from STAR.
     */
    tuple val(sample), path(bam), path(solo)


    output:

    /*
     * BAM, its index, and the Solo directory needs to be passed on.
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