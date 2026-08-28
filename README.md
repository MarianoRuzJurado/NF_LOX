# Nextflow LOX Pipeline

Nextflow pipeline for allele-specific X-chromosome analysis and X-escape characterization from single-cell RNA-seq data.

    git clone https://github.com/MarianoRuzJurado/NF_LOX

## Overview

NF_LOX automates the workflow used to generate cell-level X-chromosome haplotype and X-escape information from scRNA-seq data.
With the information of the active haplotype as well as with detection of X-escape genes allele location a loss-of-X can
be derived per cell. Please be aware that sequencing depth is the crucial for the confidence of the analysis.

The pipeline combines:

1. STARsolo alignment of the FASTQ files.
2. BAM indexing with samtools.
3. X-chromosome variant calling with bcftools.
4. Variant normalization.
5. Selection of informative/common X-linked variants.
6. WASP-based filtering of allele-specific reads.
7. A second STARsolo alignment using the WASP-filtered reads.
8. XCISE analysis for X-inactivation state and dosage.
9. Identification of heterozygous SNVs located in X-escape genes.
10. Construction of cell-level X-escape haplotypes.
11. Integration of X-escape haplotypes with XCISE dosage and LoX haplotypes.
12. Generation of the final cell-level X-escape table.

The main final outputs are:
    
    <sample>_Xescape_final.rds
    <sample>_Xescape_final_cell_level.tsv

By default, it is found in:

    NF_LOX_results/

For example:
    
    NF_LOX_results/5kPBMC4_Xescape_final.rds
    NF_LOX_results/5kPBMC4_Xescape_final_cell_level.tsv


## Usage

A typical run is:

    nextflow run /path/to/NF_LOX_Pycharm/main.nf \
        --sample SAMPLENAME \
        --r1 /path/to/L001_R1.fastq.gz,/path/to/L002_R1.fastq.gz \
        --r2 /path/to/L001_R2.fastq.gz,/path/to/L002_R2.fastq.gz \
        --reference /path/to/referencegenome/ \
        --whitelist /path/to/barcodes.txt \

The R1 and R2 arguments are comma-separated lists dependent on the number of lanes used for the experiment. The number of R1 and R2 files must be identical.


## Example

For the 5k PBMC Donor 4 dataset:
https://www.10xgenomics.com/datasets/5k_Human_Donor4_PBMC_3p_gem-x

    nextflow run /mnt/mariano/Temporary_Restore/Mariano/forLukasZ/2026/request_260520_LOX_workflow/NF_LOX_Pycharm/main.nf \
        --sample 5kPBMC4 \
        --r1 /path/to/S4_L001_R1_001.fastq.gz,/path/to/S4_L002_R1_001.fastq.gz,/path/to/S4_L003_R1_001.fastq.gz,/path/to/S4_L004_R1_001.fastq.gz \
        --r2 /path/to/S4_L001_R2_001.fastq.gz,/path/to/S4_L002_R2_001.fastq.gz,/path/to/S4_L003_R2_001.fastq.gz,/path/to/S4_L004_R2_001.fastq.gz \
        --reference /path/to/GRCh38_STAR_ensembl_115/ \
        --whitelist /opt/cellranger-10.1.0/lib/python/cellranger/barcodes/3M-3pgex-may-2023_TRU.txt \

## Output

The final cell-level result and a .rds file are returned, these one can be integrated in downstream single cell analysis for example in Seurat or Scanpy.

Default:

    NF_LOX_results/<sample>_Xescape_final.rds
    NF_LOX_results/<sample>_Xescape_final_cell_level.tsv

Intermediate files remain in the Nextflow `work/` directory of the pipeline and can be inspected for problems.
Failures of the pipeline can be mostly explained by no SNVs found in the .vcf file.
## Requirements

The following software is required to run the pipeline:

- Nextflow 26.04.6 or a compatible recent Nextflow version
- Conda

The followining software will be handled automatically in Conda environments:
- STAR
- samtools
- bcftools
- Python 3.11
- R 4.6
- Python packages:
  - pysam
  - pandas
- R packages:
  - SeuratObject
  - optparse

Dependencies are fixed as follows:

    conda-forge::python=3.11
    conda bioconda::star=2.7.11b
    conda-forge::gzip
    conda-forge::pandas

    bioconda::bcftools=1.22
    bioconda::pysam
    bioconda::samtools=1.22.1
    bioconda::bedtools=2.31.1

## Input files

The pipeline requires:

- Single cell RNA sequencing 10x Genomics FASTQ files
- A STAR genome index
- The corresponding genome FASTA provided in the reference
- A compatible 10x barcode whitelist

Make sure that the reference directory contains at least:

    GRCh38_STAR_ensembl_115/ #or other reference
    ├── STAR index files
    └── Homo_sapiens.GRCh38.dna.primary_assembly.fa

The corresponding FASTA index should also be available:

    Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai

The supplied resources directory contains pipeline-specific files, including:

    XCI_escape_gene_coordinates.bed
    GRCh38_Common_Xlinked_SNVs.repaired.vcf.gz
    GRCh38_Common_Xlinked_SNVs.repaired.vcf.gz.csi
    vcf4wasp.pl
    build_Xescape_Haplotype.py
    build_Xescape_final.py
    create_Xescape_rds.R

`XCI_escape_gene_coordinates.bed` is a stable resource used to identify variants falling within X-escape genes.

## 10x chemistry detection

NF_LOX determines the 10x library layout from the basename of the supplied whitelist.

Currently supported whitelists include:

    737k-april-2014_rc.txt
    737k-august-2016.txt
    3M-february-2018.txt
    3M-february-2018_TRU.txt
    3M-3pgex-may-2023.txt
    3M-3pgex-may-2023_TRU.txt

NOTE: STARsolo requires uncompressed barcode whitelists, please make sure no .gz files are provided.


## Important resource files

These resources are required by the workflow and should remain available:

    resources/
    ├── XCI_escape_gene_coordinates.bed
    ├── GRCh38_Common_Xlinked_SNVs.repaired.vcf.gz
    ├── GRCh38_Common_Xlinked_SNVs.repaired.vcf.gz.csi
    ├── vcf4wasp.pl
    ├── build_Xescape_Haplotype.py
    ├── build_Xescape_final.py
    └── create_Xescape_rds.R

The X-escape BED file defines the genomic coordinates used to determine whether an informative X-linked SNV lies within an X-escape gene.

## Troubleshooting

Please refer to the Issues section on GitHub to report problems.
Alternatively, an email can be sent to ruzjurado@med.uni-frankfurt.de.

## Directory structure

The recommended project structure should follow this:

    NF_LOX_Pycharm/
    ├── main.nf
    ├── modules/
    │   ├── star_align.nf
    │   ├── samtools_index.nf
    │   ├── call_variants.nf
    │   ├── normalize_variants.nf
    │   ├── isec_variants.nf
    │   ├── wasp_vcf.nf
    │   ├── star_wasp.nf
    │   ├── samtools_index_wasp.nf
    │   ├── xcise_implementation.nf
    │   ├── xcise_xdosage.nf
    │   ├── vcf_to_x_snv_bed.nf
    │   ├── intersect_xescape_snvs.nf
    │   ├── build_Xescape_haplotype.nf
    │   ├── build_summary.nf
    │   └── create_rds.nf
    ├── resources/
    │   ├── XCI_escape_gene_coordinates.bed
    │   ├── GRCh38_Common_Xlinked_SNVs.repaired.vcf.gz
    │   ├── GRCh38_Common_Xlinked_SNVs.repaired.vcf.gz.csi
    │   ├── vcf4wasp.pl
    │   ├── build_Xescape_Haplotype.py
    │   ├── build_Xescape_final.py
    │   └── create_Xescape_rds.R
    └── NF_LOX_results/

## Author

Mariano Ruz Jurado
