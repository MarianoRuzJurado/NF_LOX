#!/usr/bin/env python3

import argparse
import pysam
from collections import defaultdict


# ------------------------------------------------------------
# Command-line arguments
# ------------------------------------------------------------

parser = argparse.ArgumentParser(
    description=(
        "Assign cell barcodes to X-escape haplotypes based on "
        "allele observations at informative X-escape SNVs."
    )
)

parser.add_argument(
    "--bam",
    required=True,
    help="Input WASP-filtered BAM file."
)

parser.add_argument(
    "--snv-file",
    required=True,
    help=(
        "X-escape SNV TSV file. Expected columns include: "
        "chrom, start, end, ref, alt, ..., gene."
    )
)

parser.add_argument(
    "--output",
    required=True,
    help="Output TSV file."
)

parser.add_argument(
    "--min-mapq",
    type=int,
    default=20,
    help="Minimum mapping quality. Default: 20."
)

parser.add_argument(
    "--min-bq",
    type=int,
    default=20,
    help="Minimum base quality. Default: 20."
)

parser.add_argument(
    "--min-informative",
    type=int,
    default=1,
    help=(
        "Minimum number of informative allele observations "
        "required for haplotype assignment. Default: 1."
    )
)

args = parser.parse_args()


# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------

BAM_FILE = args.bam
SNV_FILE = args.snv_file
OUT_FILE = args.output

MIN_MAPQ = args.min_mapq
MIN_BQ = args.min_bq
MIN_INFORMATIVE = args.min_informative


# ------------------------------------------------------------
# Read X-escape SNVs
# ------------------------------------------------------------

snvs = []

with open(SNV_FILE) as f:

    for line_number, line in enumerate(f, start=1):

        line = line.strip()

        if not line:
            continue

        # Skip possible header
        if line.startswith("#"):
            continue

        fields = line.split("\t")

        if len(fields) < 5:
            raise ValueError(
                f"SNV file line {line_number} has fewer than 5 columns."
            )

        chrom = fields[0]

        try:
            start = int(fields[1])
            end = int(fields[2])
        except ValueError:
            # Allows a normal TSV header such as:
            # chrom start end ref alt ...
            if line_number == 1:
                continue
            raise ValueError(
                f"Could not parse coordinates on line {line_number}: {line}"
            )

        ref = fields[3].upper()
        alt = fields[4].upper()

        # Current X_escape_SNVs.tsv structure:
        #
        # 0 chrom
        # 1 start
        # 2 end
        # 3 ref
        # 4 alt
        # 5 genotype / additional field
        # ...
        # 8 gene
        #
        # If gene is unavailable, use ".".
        gene = fields[8] if len(fields) > 8 else "."

        # The current TSV uses:
        # start = 0-based
        # end   = 1-based
        #
        # pysam pileup uses 0-based reference positions.
        pos0 = start

        snvs.append(
            {
                "chrom": chrom,
                "pos0": pos0,
                "pos1": pos0 + 1,
                "ref": ref,
                "alt": alt,
                "gene": gene,
            }
        )


print(f"X-escape SNVs: {len(snvs)}")


# ------------------------------------------------------------
# Count alleles per cell
# ------------------------------------------------------------

counts = defaultdict(
    lambda: {
        "ref": 0,
        "alt": 0,
        "other": 0,
        "informative": 0,
        "snvs": set(),
        "genes": set(),
    }
)


# ------------------------------------------------------------
# Open BAM
# ------------------------------------------------------------

print(f"Input BAM: {BAM_FILE}")

bam = pysam.AlignmentFile(BAM_FILE, "rb")


# ------------------------------------------------------------
# Process every informative SNV
# ------------------------------------------------------------

for snv_number, snv in enumerate(snvs, start=1):

    chrom = snv["chrom"]
    pos0 = snv["pos0"]
    ref = snv["ref"]
    alt = snv["alt"]
    gene = snv["gene"]

    for column in bam.pileup(
        chrom,
        pos0,
        pos0 + 1,
        truncate=True,
        stepper="samtools",
        min_base_quality=MIN_BQ,
        min_mapping_quality=MIN_MAPQ,
    ):

        if column.reference_pos != pos0:
            continue

        for pileupread in column.pileups:

            # Skip deletions and splice junctions
            if pileupread.is_del or pileupread.is_refskip:
                continue

            read = pileupread.alignment

            # Require a cell barcode
            if not read.has_tag("CB"):
                continue

            cb = read.get_tag("CB")

            # Ignore unassigned barcode
            if cb == "-":
                continue

            query_pos = pileupread.query_position

            if query_pos is None:
                continue

            base = read.query_sequence[query_pos].upper()

            # ------------------------------------------------
            # Reference allele
            # ------------------------------------------------

            if base == ref:

                counts[cb]["ref"] += 1
                counts[cb]["informative"] += 1
                counts[cb]["snvs"].add(snv["pos1"])
                counts[cb]["genes"].add(gene)

            # ------------------------------------------------
            # Alternative allele
            # ------------------------------------------------

            elif base == alt:

                counts[cb]["alt"] += 1
                counts[cb]["informative"] += 1
                counts[cb]["snvs"].add(snv["pos1"])
                counts[cb]["genes"].add(gene)

            # ------------------------------------------------
            # Other base
            # ------------------------------------------------

            else:

                counts[cb]["other"] += 1

    # Progress indicator
    if snv_number % 10 == 0 or snv_number == len(snvs):
        print(
            f"Processed {snv_number}/{len(snvs)} X-escape SNVs",
            flush=True
        )


bam.close()


# ------------------------------------------------------------
# Classification
# ------------------------------------------------------------

print("Writing output...")


with open(OUT_FILE, "w") as out:

    out.write(
        "Barcode\t"
        "Xescape_Haplotype\t"
        "Ref_Count\t"
        "Alt_Count\t"
        "Other_Count\t"
        "Informative_Count\t"
        "Informative_SNVs\t"
        "Informative_Genes\n"
    )

    for cb in sorted(counts):

        x = counts[cb]

        ref_count = x["ref"]
        alt_count = x["alt"]
        other_count = x["other"]
        informative = x["informative"]

        # ----------------------------------------------------
        # Haplotype assignment
        # ----------------------------------------------------

        if informative < MIN_INFORMATIVE:

            classification = "Low_coverage"

        elif ref_count > 0 and alt_count > 0:

            classification = "Both"

        elif ref_count > 0:

            classification = "Haplotype_A"

        elif alt_count > 0:

            classification = "Haplotype_B"

        else:

            classification = "Low_coverage"

        # ----------------------------------------------------
        # Informative SNVs
        # ----------------------------------------------------

        snv_string = ",".join(
            str(pos)
            for pos in sorted(x["snvs"])
        )

        # ----------------------------------------------------
        # Informative genes
        # ----------------------------------------------------

        gene_string = ",".join(
            sorted(x["genes"])
        )

        # ----------------------------------------------------
        # Write row
        # ----------------------------------------------------

        out.write(
            f"{cb}\t"
            f"{classification}\t"
            f"{ref_count}\t"
            f"{alt_count}\t"
            f"{other_count}\t"
            f"{informative}\t"
            f"{snv_string}\t"
            f"{gene_string}\n"
        )


# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

print()
print("Finished.")
print(f"Output: {OUT_FILE}")
print(f"Cells with informative observations: {len(counts):,}")
print(f"Minimum informative observations: {MIN_INFORMATIVE}")
print(f"Minimum mapping quality: {MIN_MAPQ}")
print(f"Minimum base quality: {MIN_BQ}")
