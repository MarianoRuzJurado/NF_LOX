#!/usr/bin/env python3

import argparse
import pandas as pd
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(
        description="Build an integrated cell-level X-escape/XCI haplotype table."
    )

    p.add_argument(
        "--xescape-haplotype",
        required=True,
        help="Cell-level Xescape haplotype TSV."
    )
    p.add_argument(
        "--xdosage",
        required=True,
        help="Cell-level X-dosage TSV."
    )
    p.add_argument(
        "--lox-haplotype",
        required=True,
        help="Cell-level LoX haplotype TSV."
    )
    p.add_argument(
        "--xescape-snvs",
        required=True,
        help="X-escape SNV annotation TSV."
    )
    p.add_argument(
        "--output",
        required=True,
        help="Output TSV file."
    )

    return p.parse_args()


def main():

    args = parse_args()

    # ------------------------------------------------------------
    # Read input files
    # ------------------------------------------------------------

    print("Reading Xescape haplotype file...")
    hap = pd.read_csv(args.xescape_haplotype, sep="\t")

    print("Reading X-dosage file...")
    xd = pd.read_csv(args.xdosage, sep="\t")

    print("Reading LoX haplotype file...")
    lo = pd.read_csv(args.lox_haplotype, sep="\t")

    print("Reading X-escape SNV annotation...")
    snv = pd.read_csv(
        args.xescape_snvs,
        sep="\t",
        header=None,
        names=[
            "chrom",
            "start",
            "end",
            "ref",
            "alt",
            "gene_chrom",
            "gene_start",
            "gene_end",
            "gene",
            "strand",
        ],
    )

    # ------------------------------------------------------------
    # Standardize barcode column
    # ------------------------------------------------------------

    if "Barcode" not in hap.columns:
        if "CB" in hap.columns:
            hap = hap.rename(columns={"CB": "Barcode"})
        else:
            raise ValueError(
                "Xescape haplotype file has neither 'Barcode' nor 'CB'."
            )

    if "CB" in xd.columns:
        xd = xd.rename(columns={"CB": "Barcode"})

    if "CB" in lo.columns:
        lo = lo.rename(columns={"CB": "Barcode"})

    if "Barcode" not in xd.columns:
        raise ValueError("X-dosage file has no Barcode/CB column.")

    if "Barcode" not in lo.columns:
        raise ValueError("LoX haplotype file has no Barcode/CB column.")

    # ------------------------------------------------------------
    # Build active/inactive haplotypes from LoX classification
    # ------------------------------------------------------------

    lo_map = lo[
        [
            "Barcode",
            "Haplotype",
            "Haplotype_A_UMI",
            "Haplotype_B_UMI",
        ]
    ].copy()

    valid_haplotypes = {"Haplotype_A", "Haplotype_B"}

    lo_map["Active_Haplotype"] = lo_map["Haplotype"].where(
        lo_map["Haplotype"].isin(valid_haplotypes)
    )

    lo_map["Inactive_Haplotype"] = lo_map["Active_Haplotype"].map(
        {
            "Haplotype_A": "Haplotype_B",
            "Haplotype_B": "Haplotype_A",
        }
    )

    # ------------------------------------------------------------
    # Build lookup for X-escape SNVs
    #
    # The Xescape_Haplotype file contains positions as 1-based
    # coordinates in Informative_SNVs.
    #
    # Xescape_SNVs will be represented as:
    #
    #   24055644:C>T;41341173:C>T
    #
    # ------------------------------------------------------------

    snv_position = {}
    snv_gene = {}

    for _, row in snv.iterrows():

        position = str(int(row["end"]))

        variant = (
            f"{position}:"
            f"{row['ref']}>{row['alt']}"
        )

        snv_position[position] = variant
        snv_gene[position] = str(row["gene"])

    # ------------------------------------------------------------
    # Add explicit SNV descriptions to Xescape haplotype table
    # ------------------------------------------------------------

    def expand_snvs(value):

        if pd.isna(value) or str(value).strip() == "":
            return "NA"

        positions = [
            x.strip()
            for x in str(value).split(",")
            if x.strip()
        ]

        return ";".join(
            snv_position.get(pos, pos)
            for pos in positions
        )

    def expand_genes(value):

        if pd.isna(value) or str(value).strip() == "":
            return "NA"

        positions = [
            x.strip()
            for x in str(value).split(",")
            if x.strip()
        ]

        genes = []

        for pos in positions:
            gene = snv_gene.get(pos)

            if gene and gene != "nan":
                genes.append(gene)

        return ";".join(dict.fromkeys(genes)) if genes else "NA"

    hap["Xescape_SNVs"] = hap["Informative_SNVs"].apply(expand_snvs)

    hap["Xescape_SNV_Genes"] = hap["Informative_SNVs"].apply(
        expand_genes
    )

    # ------------------------------------------------------------
    # Select Xescape columns
    # ------------------------------------------------------------

    hap_columns = [
        "Barcode",
        "Xescape_Haplotype",
        "Ref_Count",
        "Alt_Count",
        "Other_Count",
        "Informative_Count",
        "Informative_SNVs",
        "Informative_Genes",
        "Xescape_SNVs",
        "Xescape_SNV_Genes",
    ]

    missing = [
        c for c in hap_columns
        if c not in hap.columns
    ]

    if missing:
        raise ValueError(
            "Missing columns in Xescape haplotype file: "
            + ", ".join(missing)
        )

    hap = hap[hap_columns]

    # ------------------------------------------------------------
    # Select X-dosage columns
    # ------------------------------------------------------------

    dosage_columns = [
        "Barcode",
        "Total_Counts",
        "Escape_X_Counts",
        "Escape_X_Genes",
        "Escape_X_CP10k",
    ]

    missing = [
        c for c in dosage_columns
        if c not in xd.columns
    ]

    if missing:
        raise ValueError(
            "Missing columns in X-dosage file: "
            + ", ".join(missing)
        )

    xd = xd[dosage_columns]

    # ------------------------------------------------------------
    # Merge
    # ------------------------------------------------------------

    print("Merging LoX haplotypes...")
    final = lo_map.merge(
        hap,
        on="Barcode",
        how="outer",
    )

    print("Merging X-dosage...")
    final = final.merge(
        xd,
        on="Barcode",
        how="outer",
    )

    # ------------------------------------------------------------
    # Reorder columns
    # ------------------------------------------------------------

    columns = [
        "Barcode",

        # XCI / LoX classification
        "Active_Haplotype",
        "Inactive_Haplotype",
        "Haplotype",
        "Haplotype_A_UMI",
        "Haplotype_B_UMI",

        # X-escape SNV classification
        "Xescape_Haplotype",
        "Ref_Count",
        "Alt_Count",
        "Other_Count",
        "Informative_Count",
        "Informative_SNVs",
        "Informative_Genes",
        "Xescape_SNVs",
        "Xescape_SNV_Genes",

        # X-escape expression / dosage
        "Total_Counts",
        "Escape_X_Counts",
        "Escape_X_Genes",
        "Escape_X_CP10k",
    ]

    columns = [
        c for c in columns
        if c in final.columns
    ]

    final = final[columns]

    # ------------------------------------------------------------
    # Sort and write
    # ------------------------------------------------------------

    final = final.sort_values("Barcode")

    output = Path(args.output)
    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    print(f"Writing: {output}")

    final.to_csv(
        output,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    # ------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------

    print()
    print("Finished.")
    print(f"Rows:    {len(final):,}")
    print(f"Columns: {len(final.columns)}")
    print()

    if "Haplotype" in final.columns:
        print("LoX haplotype:")
        print(
            final["Haplotype"]
            .value_counts(dropna=False)
            .to_string()
        )
        print()

    if "Xescape_Haplotype" in final.columns:
        print("Xescape haplotype:")
        print(
            final["Xescape_Haplotype"]
            .value_counts(dropna=False)
            .to_string()
        )


if __name__ == "__main__":
    main()
