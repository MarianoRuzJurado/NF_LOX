process ISEC_VARIANTS {
    // Right now the bcftools isec uses Qual >= 200, might be to strict
    tag "$sample"
    conda 'bioconda::bcftools=1.22'

    input:

    tuple val(sample), path(vcf), path(vcf_index)

    path common_snvs
    path common_snvs_index


    output:

    tuple val(sample),
          path("${sample}common_SNVs"),
          emit: isec_output


    script:

    """
    mkdir ${sample}common_SNVs

    bcftools isec \
        -i 'Qual>=200' \
        ${vcf} \
        ${common_snvs} \
        -p ${sample}common_SNVs
    """
}