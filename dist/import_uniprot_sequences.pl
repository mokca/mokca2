#!/usr/bin/perl -w

use Getopt::Long;
use FileHandle;
use DBI();

use database_setup;
use path_setup;

# Autoflush STDERR

STDERR->autoflush(1);

# Ticker

my $tick_c = 0;
my $col_wrap = 78;
my $tick_step = 10;
my $tick_char = '.';

# Explain purpose

print STDERR "::: Import Uniprot sequences - Import sequences for all uniprot entries, including isoforms\n";

# Setup universal paths, database details

my $sprot_seq_fn = $mokcatanic_data_root . "uniprot_sprot.fasta.gz";
my $trembl_seq_fn = $mokcatanic_data_root . "uniprot_trembl.fasta.gz";
my $uniprot_iso_fn = $mokcatanic_data_root . "uniprot_sprot_varsplic.fasta.gz";

# Connect to databases

my $dbh_titanic = DBI->connect("DBI:mysql:database=$mokca_db;host=$mokca_host",$mokca_user,$mokca_pw, {'RaiseError' => 1});

# Prepare queries

my $dbq_accs = $dbh_titanic->prepare(q{SELECT DISTINCT accession FROM uniprot});

# Go

my %accessions;
my ( $uniprot_accession, $uniprot_isoform );
my $sequence;

open my $uniprot_fh, "gunzip -c $sprot_seq_fn |" || die "\n+++ Cannot open $sprot_seq_fn to read data: $?\n";
while ($line = <$uniprot_fh>) {
    chomp $line;
    if ($line =~ /^>..\|(.{6})\|/) {
        $uniprot_accession = $1;
        $accessions{$uniprot_accession} = 1;
    }
}
close $uniprot_fh || die "\n+++ Cannot close $sprot_seq_fn to stop reading data: $?\n";

open $uniprot_fh, "gunzip -c $trembl_seq_fn |" || die "\n+++ Cannot open $trembl_seq_fn to read data: $?\n";
while ($line = <$uniprot_fh>) {
    chomp $line;
    if ($line =~ /^>..\|(.{6})\|/) {
        $uniprot_accession = $1;
        $accessions{$uniprot_accession} = 1;
    }
}
close $uniprot_fh || die "\n+++ Cannot close $trembl_seq_fn to stop reading data: $?\n";

$dbq_accs->execute();
while ((my $acc) = $dbq_accs->fetchrow_array) {
    if (!defined $accessions{$acc}) {
        print "--- Missing sequence for $acc\n";
    }
}
