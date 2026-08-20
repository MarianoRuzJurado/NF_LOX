process INTERSECT_XESCAPE_SNVS {

    tag "$sample"
    conda 'bioconda::bcftools=2.31.1'

    input:
    tuple val(sample), path(snv_bed), path(xescape_bed)

    output:
    tuple val(sample),
          path("${sample}_Xescape_SNVs.tsv"),
          emit: xescape_snvs

    script:
    """
    bedtools intersect \
        -a ${snv_bed} \
        -b ${xescape_bed} \
        -wa -wb \
        > ${sample}_Xescape_SNVs.tsv
    """
}