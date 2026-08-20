process BUILD_SUMMARY {

    tag "$sample"

    publishDir "${params.outdir}",
    mode: 'copy',
    overwrite: true,
    saveAs: { filename ->
        "${sample}/${filename}"
    }

    conda 'conda-forge::python=3.11 conda-forge::pandas'

    input:
    tuple val(sample),
          path(xescape_haplotype),
          path(xdosage),
          path(lox_haplotype),
          path(xescape_snvs)

    output:
    tuple val(sample),
          path("${sample}_Xescape_final_cell_level.tsv"),
          emit: xescape_final

    script:
    """
    python ${params.xcise_dir}/build_Xescape_final.py \
        --xescape-haplotype ${xescape_haplotype} \
        --xdosage ${xdosage} \
        --lox-haplotype ${lox_haplotype} \
        --xescape-snvs ${xescape_snvs} \
        --output ${sample}_Xescape_final_cell_level.tsv
    """
}