process WASP_VCF {

    tag "$sample"

    input:

    tuple val(sample), path(isec_dir), path(vcf4wasp)


    output:

    tuple val(sample),
          path("${sample}heterozygous_SNVs.vcf"),
          emit: wasp_vcf


    script:

    """
    perl ${vcf4wasp} \
        ${isec_dir}/0002.vcf \
        > ${sample}heterozygous_SNVs.vcf
    """
}