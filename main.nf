nextflow.enable.dsl=2


include { STAR_ALIGN } from './modules/star_align'
include { SAMTOOLS_INDEX } from './modules/samtools_index'
include { CALL_VARIANTS } from './modules/call_variants'
include { NORMALIZE_VARIANTS } from './modules/normalize_variants'
include { ISEC_VARIANTS } from './modules/isec_variants'
include { WASP_VCF } from './modules/wasp_vcf'
include { STAR_WASP} from './modules/star_wasp'
include { SAMTOOLS_INDEX_WASP } from './modules/samtools_index_wasp'
include { XCISE_IMPLEMENTATION } from './modules/xcise_implementation'
include { XCISE_XDOSAGE } from './modules/xcise_xdosage'

/*
 * Automatic chemistry detection with barcode file
 */

def whitelist_layouts() {
    return [
        '737k-april-2014_rc.txt': [
            name: '10x_3prime_v1',
            cb_start: 1,
            cb_len: 16,
            umi_start: 17,
            umi_len: 10
        ],

        '737k-august-2016.txt': [
            name: '10x_3prime_v2_or_5prime_v1_v2',
            cb_start: 1,
            cb_len: 16,
            umi_start: 17,
            umi_len: 10
        ],

        '3M-february-2018.txt': [
            name: '10x_3prime_v3_v3.1',
            cb_start: 1,
            cb_len: 16,
            umi_start: 17,
            umi_len: 12
        ],

         '3M-february-2018_TRU.txt': [
            name: '10x_3prime_v3_v3.1',
            cb_start: 1,
            cb_len: 16,
            umi_start: 17,
            umi_len: 12
        ],

        '3M-3pgex-may-2023.txt': [
            name: '10x_3prime_v4',
            cb_start: 1,
            cb_len: 16,
            umi_start: 17,
            umi_len: 12
        ],

        '3M-3pgex-may-2023_TRU.txt': [
            name: '10x_3prime_v4',
            cb_start: 1,
            cb_len: 16,
            umi_start: 17,
            umi_len: 12
        ]
    ]
}

workflow {

    println '''
╔══════════════════════════════════════════════════════════════╗
║                                                              ║
║    ███╗   ██╗███████╗       ██╗      ██████╗ ██╗  ██╗        ║
║    ████╗  ██║██╔════╝       ██║     ██╔═══██╗╚██╗██╔╝        ║
║    ██╔██╗ ██║█████╗         ██║     ██║   ██║ ╚███╔╝         ║
║    ██║╚██╗██║██╔══╝         ██║     ██║   ██║ ██╔██╗         ║
║    ██║ ╚████║██║            ███████╗╚██████╔╝██╔╝ ██╗        ║
║    ╚═╝  ╚═══╝╚═╝            ╚══════╝ ╚═════╝ ╚═╝  ╚═╝        ║
║                                                              ║
║══════════════════════════════════════════════════════════════║
║    Nextflow LOX Pipeline - Author: Mariano Ruz Jurado        ║
║    Found problems? Report them under NF_LOX/Issues on github ║
╚══════════════════════════════════════════════════════════════╝
    '''

    /*
     * ---------------------------------------------------------
     * Validate command-line parameters
     * ---------------------------------------------------------
     */

    println "PARAMETERS:"
    println params

    if (!params.sample) {
        error "Missing required parameter: --sample"
    }

    if (!params.r1) {
        error "Missing required parameter: --r1"
    }

    if (!params.r2) {
        error "Missing required parameter: --r2"
    }

    if (!params.reference) {
        error "Missing required parameter: --reference"
    }

    if (!params.whitelist) {
        error "Missing required parameter: --whitelist"
    }



    /*
     * Define reference files
     */

    star_index = file("${params.reference}")

    genome_fasta = file(
        "${params.reference}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    )

    whitelist = file(params.whitelist)
    if (whitelist.getName().endsWith('.gz')) {
        log.error """
        Error: Compressed whitelist detected:

            ${whitelist}

        STARsolo requires an uncompressed barcode whitelist.

        Please decompress it before running the pipeline, for example:

            gunzip ${whitelist.getName()}

        The pipeline will stop.
        """
    }

    common_snvs = file(params.common_snvs)
    common_snvs_index = file("${params.common_snvs}.csi")

    xcise_dir = file(params.xcise_dir)
    vcf4wasp = file("${params.xcise_dir}/vcf4wasp.pl")

    /*
     * Detect 10x barcode and UMI layout
     */

    def layouts = whitelist_layouts()
    def whitelist_name = whitelist.getName()

    if (!layouts.containsKey(whitelist_name)) {
        error """
    Unsupported 10x whitelist:

    ${whitelist_name}

    Supported whitelists:
    ${layouts.keySet().join('\n')}
    """
    }

    def chemistry = layouts[whitelist_name]

    println """
    ---------------------------------------------------------
    Detected 10x library layout
    ---------------------------------------------------------
    Whitelist:  ${whitelist_name}
    Chemistry:  ${chemistry.name}
    CB start:   ${chemistry.cb_start}
    CB length:  ${chemistry.cb_len}
    UMI start:  ${chemistry.umi_start}
    UMI length: ${chemistry.umi_len}
    ---------------------------------------------------------
    """

    /*
     * Check for required files and directories
     */

    if (!star_index.exists()) {
        error """
        STAR index directory does not exist:

        ${star_index}

        Please check --reference.
        """
    }

    if (!genome_fasta.exists()) {
        error """
        Reference FASTA does not exist:

        ${genome_fasta}

        Please check --reference.
        """
    }

    if (!whitelist.exists()) {
        error """
        Barcode whitelist does not exist:

        ${whitelist}

        Please check --whitelist.
        """
    }

    if (!common_snvs.exists()) {
        error """
        Common SNV file does not exist:

        ${common_snvs}
        """
    }

    if (!common_snvs_index.exists()) {
        error """
        Common SNV VCF index does not exist:

        ${common_snvs_index}

        Please check the .csi index for the common SNV VCF.

        """
    }

    if (!vcf4wasp.exists()) {
        error """
        vcf4wasp.pl file does not exist:

        ${vcf4wasp}

        Please check the XCISE resource directory.
        """
    }

    /*
     * ---------------------------------------------------------
     * FASTQ input
     * preserve the R2/R1 order from original XCISE analysis.
     * This step still confuses me, documentation is missing about r2 or r1 first.
     *
     * Original:
     *
     * STAR --readFilesIn R2.fastq.gz R1.fastq.gz
     * has now support for multilane samples, needs testing
     */

    r1_files = params.r1
        .split(',')
        .collect { file(it.trim()) }

    r2_files = params.r2
        .split(',')
        .collect { file(it.trim()) }


    if (r1_files.size() != r2_files.size()) {
        error """
        The number of R1 and R2 FASTQ files does not match.

        R1 files: ${r1_files.size()}
        R2 files: ${r2_files.size()}

        Please provide the same number of R1 and R2 files.
        """
    }

    reads_ch = Channel.of(
        tuple(
            params.sample,
            r2_files,
            r1_files
        )
    )


    /*
     * STARsolo first run without WASP -> normal alignment
     */

    star_ch = STAR_ALIGN(
        reads_ch,
        star_index,
        whitelist,
        chemistry
    )


    /*
     * samtools index
     */

    indexed_bam_ch = SAMTOOLS_INDEX(
        star_ch
    )


    /*
     * Variant calling
     */

    variants_ch = CALL_VARIANTS(
        indexed_bam_ch,
        genome_fasta
    )


    /*
     * Variant normalization
     */

    normalized_ch = NORMALIZE_VARIANTS(
        variants_ch,
        genome_fasta
    )

    /*
     * intersction of variants
     */

    isec_ch = ISEC_VARIANTS(
        normalized_ch,
        common_snvs,
        common_snvs_index
    )

    wasp_vcf_input_ch = isec_ch.map { sample, isec_dir ->
        tuple(sample, isec_dir, vcf4wasp)
    }

    heterozygous_vcf_ch = WASP_VCF(
        wasp_vcf_input_ch
    )

    reads_ch = Channel.of(
        tuple(
            params.sample,
            r2_files,
            r1_files
        )
    )

    wasp_reads_ch = Channel.of(
        tuple(
            params.sample,
            r2_files,
            r1_files
        )
    )

    wasp_bam_ch = STAR_WASP(
        wasp_reads_ch,
        star_index,
        whitelist,
        chemistry,
        heterozygous_vcf_ch
    )

    wasp_indexed_ch = SAMTOOLS_INDEX_WASP(
        wasp_bam_ch
    )

    xcise_input_ch = wasp_indexed_ch
        .join(heterozygous_vcf_ch)
        .map { sample, bam, bai, solo, vcf ->
            tuple(sample, bam, bai, solo, vcf)
    }

    xcise_ch = XCISE_IMPLEMENTATION(
        xcise_input_ch
    )

    xcise_xdosage_ch = XCISE_XDOSAGE(
        xcise_ch
    )

}