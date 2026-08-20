process CALL_VARIANTS {

    tag "$sample"
    conda 'bioconda::bcftools=1.22'

    input:

    /*
     * BAM produced by STAR and indexed by SAMtools.
     *
     * Solo directory need to be provided so it is
     * carried along for the workflow to use it
     * later.
     */
    tuple val(sample), path(bam), path(bai), path(solo)

    /*
     * Reference genome.
     */
    path genome_fasta


    output:

    tuple val(sample),
          path("${sample}variants.vcf.gz"),
          emit: variants


    script:

    """
    bcftools mpileup \
        -a AD \
        -f ${genome_fasta} \
        -q 20 \
        -d 1000000 \
        -r X \
        ${bam} |
    bcftools call \
        -mv \
        -Oz \
        -o ${sample}variants.vcf.gz
    """
}