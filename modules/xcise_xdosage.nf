process XCISE_XDOSAGE {

    tag "$sample"

    input:
    tuple val(sample),
          path(haplotype),
          path(solo),
          path(bam),
          path(vcf)

    output:
    tuple val(sample),
          path("${sample}_xdosage.txt"),
          path(haplotype),
          path(solo),
          path(bam),
          path(vcf),
          emit: xdosage_output

    script:
    """
    perl ${params.xcise_dir}/xcise_xdosage.pl \
        -m ${solo}/Gene/filtered/matrix.mtx \
        -f ${solo}/Gene/filtered/features.tsv \
        -c ${solo}/Gene/filtered/barcodes.tsv \
        -l ${haplotype} \
        -o ${sample}
    """
}