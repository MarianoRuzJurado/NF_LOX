process BUILD_XESCAPE_HAPLOTYPE {

    tag "$sample"

    // Add dependency from python script with conda environment
    conda 'conda-forge::python=3.11 bioconda::pysam'

    input:
    tuple val(sample),
          path(bam),
          path(bai),
          path(snv_file)

    output:
    tuple val(sample),
          path("${sample}_Xescape_Haplotype.tsv"),
          emit: xescape_haplotype

    script:
    """
    python ${params.xcise_dir}/build_Xescape_Haplotype.py \
        --bam ${bam} \
        --snv-file ${snv_file} \
        --output ${sample}_Xescape_Haplotype.tsv
    """
}
