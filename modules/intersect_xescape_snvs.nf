process INTERSECT_XESCAPE_SNVS {

    tag "$sample"

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