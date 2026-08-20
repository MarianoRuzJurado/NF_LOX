process VCF_TO_X_SNV_BED {

    tag "$sample"

    input:
    tuple val(sample), path(vcf)

    output:
    tuple val(sample),
          path("${sample}_heterozygous_X_SNVs.bed"),
          emit: x_snv_bed

    script:
    """
    bcftools query \
        -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' \
        ${vcf} \
    | awk 'BEGIN {OFS="\\t"} {
        print \$1, \$2-1, \$2, \$3, \$4
    }' \
    > ${sample}_heterozygous_X_SNVs.bed
    """
}