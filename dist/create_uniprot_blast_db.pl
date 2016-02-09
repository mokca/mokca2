#!/usr/bin/perl -w

use strict;

use Getopt::Long;
use FileHandle;
use DBI();

# Autoflush STDERR

STDERR->autoflush(1);

# Ticker

my $tick_c = 0;
my $col_wrap = 78;
my $tick_step = 10;
my $tick_char = '.';

# Explain purpose

print STDERR "::: Create UNIPROT blast database - create database for blasting from accumulated human isoform 1 sequences\n";

# Setup universal paths, database details

use path_setup;
use database_setup;

my $uniprot_fasta_fn = $mokcatanic_local_root . "uniprot_human.seq";

# Connect to databases

my $dbh_titanic = DBI->connect("DBI:mysql:database=$mokca_db;host=$mokca_host",$mokca_user,$mokca_pw, {'RaiseError' => 1});

# Prepare queries

my $dbq_sequences = $dbh_titanic->prepare(q{SELECT sequence.sequence, uniprot.accession FROM sequence, uniprot WHERE sequence.source = 'uniprot' AND sequence.isoform = 1 AND sequence.codes_for = 'protein' AND sequence.source_id = uniprot.uniprot_id});

# Go

my ($uniprot_sequence, $uniprot_accession, $fasta_fh);
open ($fasta_fh, ">$uniprot_fasta_fn") || die "\n+++ Cannot open $uniprot_fasta_fn to write sequences: $?\n";
$dbq_sequences->execute;
while (($uniprot_sequence, $uniprot_accession) = $dbq_sequences->fetchrow_array) {
    
    print $fasta_fh ">$uniprot_accession\n$uniprot_sequence\n";
    
    # Pointless ticker
    
    if ($tick_c % $tick_step == 0) {
        print STDERR $tick_char;
        if ($tick_c % ($tick_step * $col_wrap) == 0 && $tick_c) {
            print STDERR "\n";
        }
    }
    $tick_c++;
}
close ($fasta_fh) || die "\n+++Cannot close $uniprot_fasta_fn: $?\n";

system("makeblastdb","-dbtype","prot","-in",$uniprot_fasta_fn);

# Tick finish

print STDERR "\n";

$dbh_titanic->disconnect;
