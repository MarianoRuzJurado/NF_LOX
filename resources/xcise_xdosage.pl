#!/usr/bin/perl -w
use strict;

# XCISE-Xdosage
#
# Second layer of the XCISE-LoX framework.
#
# Calculates per-cell expression of a predefined set of
# XCI-escape genes and normalizes it by total expression.
#
# IMPORTANT:
# This script does NOT assign maternal/paternal identity.
# It also does NOT initially force a 1X/2X biological classification.
# The X-escape expression score should first be inspected empirically.
#
# Input:
#   -m  matrix.mtx
#   -f  features.tsv
#   -c  barcodes.tsv
#   -l  XCISE-LoX bc2hap file
#   -o  output prefix
#
# Output:
#   <output>_xdosage.txt
#   <output>_xdosage_summary.txt


my ($matrix_file, $features_file, $barcodes_file, $lox_file, $sample);

my @escape_genes = qw(
    NLGN4X
    STS
    PUDP
    PNPLA4
    ANOX1
    FAMSC
    OFD1
    GPM6B
    GEMIN8
    CA5B
    ZRSR2
    AP1S2
    S100G
    CTPS2
    EIF1AX
    EIF2S3
    ZFX
    USP9X
    DDX3X
    KDM6A
    UBA1
    CDK16
    KDM5C
    IQSEC2
    SMC1A
    RPS4X
    JPX
);

warn "Reading parameters ...\n";

my $ele = 0;

while ($ele <= $#ARGV) {

    if ($ARGV[$ele] eq '-m') {
        die "Matrix file does not exist: $ARGV[$ele+1]\n"
            unless -e $ARGV[$ele+1];

        $matrix_file = $ARGV[$ele+1];
        $ele += 2;
    }

    elsif ($ARGV[$ele] eq '-f') {
        die "Features file does not exist: $ARGV[$ele+1]\n"
            unless -e $ARGV[$ele+1];

        $features_file = $ARGV[$ele+1];
        $ele += 2;
    }

    elsif ($ARGV[$ele] eq '-c') {
        die "Barcode file does not exist: $ARGV[$ele+1]\n"
            unless -e $ARGV[$ele+1];

        $barcodes_file = $ARGV[$ele+1];
        $ele += 2;
    }

    elsif ($ARGV[$ele] eq '-l') {
        die "LoX file does not exist: $ARGV[$ele+1]\n"
            unless -e $ARGV[$ele+1];

        $lox_file = $ARGV[$ele+1];
        $ele += 2;
    }

    elsif ($ARGV[$ele] eq '-o') {
        $sample = $ARGV[$ele+1];
        $ele += 2;
    }

    else {
        die "Unexpected parameter: $ARGV[$ele]\n";
    }
}


my $usage =
    "Usage: perl $0 " .
    "-m <matrix.mtx> " .
    "-f <features.tsv> " .
    "-c <barcodes.tsv> " .
    "-l <XCISE_LoX_bc2hap.txt> " .
    "-o <output_prefix>\n";


die "Missing matrix file\n$usage"
    unless $matrix_file;

die "Missing features file\n$usage"
    unless $features_file;

die "Missing barcode file\n$usage"
    unless $barcodes_file;

die "Missing LoX file\n$usage"
    unless $lox_file;

die "Missing output prefix\n$usage"
    unless $sample;


warn "Matrix   : $matrix_file\n";
warn "Features : $features_file\n";
warn "Barcodes : $barcodes_file\n";
warn "LoX file : $lox_file\n";
warn "Output   : $sample\n";


#
# ------------------------------------------------------------
# Read features
# ------------------------------------------------------------
#

warn "Reading features ...\n";

my %gene_row = ();
my %found_escape = ();

open my $FF, '<', $features_file
    or die "Cannot open features file: $features_file\n";

my $feature_number = 0;

while (my $line = <$FF>) {

    chomp $line;

    my @a = split /\t/, $line;

    $feature_number++;

    next unless @a >= 2;

    my $gene_id = $a[0];
    my $gene_name = $a[1];

    if (grep { $_ eq $gene_name } @escape_genes) {
        $gene_row{$gene_name} = $feature_number;
        $found_escape{$gene_name} = 1;
    }
}

close $FF;


warn "Total features loaded: $feature_number\n";

warn "XCI-escape genes found in matrix:\n";

foreach my $gene (@escape_genes) {

    if (exists $gene_row{$gene}) {
        warn "    $gene\trow=$gene_row{$gene}\n";
    }
    else {
        warn "    $gene\tNOT FOUND\n";
    }
}


my $n_escape_found = scalar keys %found_escape;

warn "Escape genes found: $n_escape_found / ",
     scalar(@escape_genes), "\n";

die "No XCI-escape genes were found in features.tsv\n"
    unless $n_escape_found;


#
# ------------------------------------------------------------
# Read barcodes
# ------------------------------------------------------------
#

warn "Reading barcodes ...\n";

my @barcodes = ();
my %barcode_index = ();

open my $FB, '<', $barcodes_file
    or die "Cannot open barcode file: $barcodes_file\n";

while (my $line = <$FB>) {

    chomp $line;

    next unless length($line);

    push @barcodes, $line;
    $barcode_index{$line} = $#barcodes;
}

close $FB;

my $n_cells = scalar @barcodes;

warn "Cells loaded: $n_cells\n";


#
# ------------------------------------------------------------
# Read LoX calls
# ------------------------------------------------------------
#

warn "Reading XCISE-LoX haplotype calls ...\n";

my %lox_haplotype = ();
my %lox_hapA = ();
my %lox_hapB = ();

open my $FL, '<', $lox_file
    or die "Cannot open LoX file: $lox_file\n";

my $lox_header = <$FL>;

while (my $line = <$FL>) {

    chomp $line;

    next unless length($line);

    my @a = split /\t/, $line;

    next unless @a >= 4;

    my ($cb, $ha, $hb, $hap) = @a[0..3];

    $lox_haplotype{$cb} = $hap;
    $lox_hapA{$cb} = $ha;
    $lox_hapB{$cb} = $hb;
}

close $FL;

warn "LoX cells loaded: ",
     scalar(keys %lox_haplotype), "\n";


#
# ------------------------------------------------------------
# Read MatrixMarket header
# ------------------------------------------------------------
#

warn "Reading expression matrix ...\n";

open my $FM, '<', $matrix_file
    or die "Cannot open matrix file: $matrix_file\n";


my $line;

while ($line = <$FM>) {

    chomp $line;

    next if $line =~ /^%/;

    last;
}


die "Could not read MatrixMarket dimensions\n"
    unless defined $line;


my ($matrix_rows, $matrix_cols, $matrix_entries) =
    split /\s+/, $line;


warn "Matrix dimensions:\n";
warn "    Rows    : $matrix_rows\n";
warn "    Columns : $matrix_cols\n";
warn "    Entries : $matrix_entries\n";


die "Feature count mismatch: features.tsv=$feature_number matrix=$matrix_rows\n"
    unless $feature_number == $matrix_rows;

die "Barcode count mismatch: barcodes.tsv=$n_cells matrix=$matrix_cols\n"
    unless $n_cells == $matrix_cols;


#
# ------------------------------------------------------------
# Prepare per-cell storage
# ------------------------------------------------------------
#

my @total_counts = ();
my @escape_counts = ();
my @escape_gene_number = ();

for my $i (0 .. $n_cells-1) {
    $total_counts[$i] = 0;
    $escape_counts[$i] = 0;
    $escape_gene_number[$i] = 0;
}


#
# Map row number -> escape gene
#
# This makes matrix parsing fast.
#

my %escape_row = ();

foreach my $gene (@escape_genes) {

    next unless exists $gene_row{$gene};

    $escape_row{$gene_row{$gene}} = $gene;
}


#
# ------------------------------------------------------------
# Parse sparse matrix
# ------------------------------------------------------------
#

warn "Parsing sparse matrix ...\n";

my $entries_read = 0;

while (my $line = <$FM>) {

    next if $line =~ /^\s*$/;
    next if $line =~ /^%/;

    chomp $line;

    my ($row, $col, $value) = split /\s+/, $line;

    next unless defined $value;

    $row--;
    $col--;

    $total_counts[$col] += $value;

    if (exists $escape_row{$row+1}) {

        my $gene = $escape_row{$row+1};

        $escape_counts[$col] += $value;
        $escape_gene_number[$col]++;
    }

    $entries_read++;

    if ($entries_read % 250000 == 0) {
        warn "    Matrix entries processed: $entries_read / $matrix_entries\n";
    }
}

close $FM;

warn "Matrix entries processed: $entries_read\n";


#
# ------------------------------------------------------------
# Calculate normalized X-escape score
# ------------------------------------------------------------
#
# Counts-per-10,000 normalization:
#
# escape CP10k = escape counts / total counts * 10000
#
# This removes the major effect of different sequencing depth.
#

warn "Calculating X-escape expression scores ...\n";

open my $FO, '>',
    $sample.'_xdosage.txt'
    or die "Cannot create output file\n";


print $FO join(
    "\t",
    'CB',
    'Total_Counts',
    'Escape_X_Counts',
    'Escape_X_Genes',
    'Escape_X_CP10k',
    'LoX_Haplotype',
    'LoX_Haplotype_A_UMI',
    'LoX_Haplotype_B_UMI'
), "\n";


my $total_cells = 0;

my $sum_escape_cp10k = 0;

foreach my $i (0 .. $n_cells-1) {

    my $cb = $barcodes[$i];

    my $total = $total_counts[$i];
    my $escape = $escape_counts[$i];
    my $ngenes = $escape_gene_number[$i];

    my $escape_cp10k = 0;

    if ($total > 0) {
        $escape_cp10k = 10000 * $escape / $total;
    }

    $sum_escape_cp10k += $escape_cp10k;

    my $hap =
        exists $lox_haplotype{$cb}
        ? $lox_haplotype{$cb}
        : 'Not_in_LoX';

    my $hapA =
        exists $lox_hapA{$cb}
        ? $lox_hapA{$cb}
        : 0;

    my $hapB =
        exists $lox_hapB{$cb}
        ? $lox_hapB{$cb}
        : 0;

    print $FO join(
        "\t",
        $cb,
        $total,
        $escape,
        $ngenes,
        sprintf("%.6f", $escape_cp10k),
        $hap,
        $hapA,
        $hapB
    ), "\n";

    $total_cells++;
}

close $FO;


#
# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------
#

my $mean_escape_cp10k =
    $total_cells
    ? $sum_escape_cp10k / $total_cells
    : 0;


open my $FS, '>',
    $sample.'_xdosage_summary.txt'
    or die "Cannot create summary file\n";


print $FS "XCISE-Xdosage summary\n";
print $FS "=====================\n";
print $FS "Cells                    : $total_cells\n";
print $FS "Escape genes requested   : ", scalar(@escape_genes), "\n";
print $FS "Escape genes found       : $n_escape_found\n";
print $FS "Mean Escape-X CP10k        : $mean_escape_cp10k\n";
print $FS "\n";
print $FS "Escape genes found:\n";

foreach my $gene (@escape_genes) {

    if (exists $gene_row{$gene}) {
        print $FS "$gene\trow=$gene_row{$gene}\n";
    }
}

print $FS "\n";
print $FS "IMPORTANT:\n";
print $FS "Escape_X_CP10k is a normalized expression score. CP10k = escape counts / total counts * 10000 \n";
print $FS "It is NOT yet a calibrated 1X/2X dosage call.\n";
print $FS "Dosage thresholds should be determined empirically from the cell distribution.\n";

close $FS;


warn "\n";
warn "Finished X-dosage calculation.\n";
warn "Output:\n";
warn "    $sample"."_xdosage.txt\n";
warn "    $sample"."_xdosage_summary.txt\n";
warn "\n";
warn "Mean Escape-X CP10k: $mean_escape_cp10k\n";
