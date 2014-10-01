#!/usr/bin/perl

use Getopt::Long;
use FileHandle;

# Autoflush STDERR

STDERR->autoflush(1);

# Set up universal paths

use path_setup;

# Fetch COSMIC files: CosmicInsMutExport_vxx_xxxxxx, CosmicMutantExport_vxx_xxxxxx, CosmicHGNC_vxx_xxxxxx,
# CosmicCellLineExport_vxx_xxxxxx, cancer_gene_census

print STDERR "--- Fetch COSMIC files, release $cosmic_release\n";

my $cosmic_fn = "ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/CosmicInsMutExport_" . $cosmic_release . ".tsv.gz";

system("wget", $cosmic_fn, "-N", "-P", $mokcatanic_data_root) == 0 || die "+++Cannot fetch $cosmic_fn: $?\n";

$cosmic_fn = "ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/CosmicCompleteExport_" . $cosmic_release . ".tsv.gz";

system("wget", $cosmic_fn, "-N", "-P", $mokcatanic_data_root) == 0 || die "+++Cannot fetch $cosmic_fn: $?\n";

$cosmic_fn = "ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/CosmicMutantExport_" . $cosmic_release . ".tsv.gz";

system("wget", $cosmic_fn, "-N", "-P", $mokcatanic_data_root) == 0 || die "+++Cannot fetch $cosmic_fn: $?\n";

$cosmic_fn = "ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/CosmicHGNC_" . $cosmic_release . ".tsv.gz";

system("wget", $cosmic_fn, "-N", "-P", $mokcatanic_data_root) == 0 || die "+++Cannot fetch $cosmic_fn: $?\n";

# CellLineExport_[version] changed to CellLineProject_[version] in v59_230512.  An error?
# ...and then completely missing in v68, despite being mentioned in the README.  Skipping for now.
#
#$cosmic_fn = "ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/CosmicCellLineProject_" . $cosmic_release . ".tsv.gz"
#
#system("wget", $cosmic_fn, "-N", "-P", $mokcatanic_data_root) == 0 || die "+++Cannot fetch $cosmic_fn: $?\n";

$cosmic_fn = "ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/cancer_gene_census.tsv";

system("wget", $cosmic_fn, "-N", "-P", $mokcatanic_data_root) == 0 || die "+++Cannot fetch $cosmic_fn: $?\n";
