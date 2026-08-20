process XCISE_IMPLEMENTATION {

    tag "$sample"

    input:
    tuple val(sample), path(bam), path(bai), path(solo), path(vcf)

    output:
    tuple val(sample),
          path("${sample}_chrX_XCISE_LoX_bc2hap.txt"),
          path(solo),
          path(bam),
          path(vcf),
          emit: xcise_output

    script:
    """
    perl ${params.xcise_dir}/xcise_lox.pl \
        -o ${sample} \
        -s ${vcf} \
        -b ${bam}
    """
}