process CREATE_RDS {

    tag "$sample"

    publishDir "${params.outdir}",
        mode: 'copy',
        overwrite: true,
        saveAs: { filename ->
            "${sample}/${filename}"
        }

    conda 'conda-forge::r-base=4.6'

    input:
    tuple val(sample),
          path(solo),
          path(xescape_final)

    output:
    tuple val(sample),
          path("${sample}_Xescape_final.rds"),
          emit: rds

    script:
    """
    Rscript ${params.xcise_dir}/create_Xescape_rds.R \
        --solo ${solo} \
        --metadata ${xescape_final} \
        --output ${sample}_Xescape_final.rds
    """
}
